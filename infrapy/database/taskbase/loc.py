'''
Created on Oct 31, 2014

@author: omarcillo
'''

from .base import Base

import sys, pdb
import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
from obspy.core import UTCDateTime
from obspy.core.util import AttribDict
from datetime import datetime

from scipy import stats
import numpy as np
import scipy as sc
from IPython import embed
import pisces as ps

from sqlalchemy import func
import cmath
import math
import itertools
import pylab as py
import matplotlib

import matplotlib.pyplot as pl
import matplotlib.mlab as mpy

import time

from ...utils.cart2pol import cart2pol
from ...utils.short_time import short_time

from .. import schema
import matplotlib.dates as mdates

from ...propagation import likelihoods
from ...propagation import infrasound
from ...location import bisl


from multiprocessing import Pool
import multiprocessing

import pkgutil

import warnings

warnings.filterwarnings('error')

import glob


status=0

def run_bislSC_wrapper(args):
    '''
    wrapper for multicore and quality control
    '''
    try:
        global status
        if len(args[0])==0:
            #print 'No detections'
            return []
        #print 'start:',args[-2]

        if bisl_qc_check(args[0], args[1])==False:
            #print 'Assoc may not converge'
            st=[]
        else:
            try:
                # remove multiple phases and keep the earliest
                IDs=[]
                ts=[]
                for ee in args[0]:
                    IDs.append(ee.ID)
                    ts.append(ee.t)
                uIDs=list(set(IDs))
                uID_index=[]
                for ii in range(len(uIDs)):
                    aux1=[]
                    auxt=[]
                    for uu in range(len(IDs)):
                        if IDs[uu]==uIDs[ii]:
                            aux1.append(uu)
                            auxt.append(ts[uu])
                    auxMIN=aux1[np.argmin(auxt)]
                    uID_index.append(auxMIN)
                argsN=[]
                for ii in uID_index:
                    argsN.append(args[0][ii])
                argsP=[]
                for ii in range(len(args)):
                    if ii==0:
                        argsP.append(argsN)
                    else:
                        argsP.append(args[ii])
                st=run_bislSC(*argsP)
            except Warning:
                print('warning')
                st=[]
        status=status+1
        '''
        sys.stdout.write("\rPercent (in core): {0}%".format(round(args[5]*status*100.0/args[-1]))+"  # {0}".format(int(args[-2])))
        sys.stdout.flush()
        '''
        print(["\rPercent (in core): {0}%".format(round(args[5]*status*100.0/args[-1]))+"  # {0}".format(int(args[-2]))])
    except KeyboardInterrupt:
        print('stop at:',args[-2])

        for aa in args[0]:
            print(aa)

        #embed()
        #exit()
        return []

    return st

class LocInfraPy(Base):
    '''
    classdocs
    '''
    # this one keeps only the earliest detection per array
    algorithm='Blom/BISL III'

    def __init__(self, conf_file=[]):
        '''
        Constructor
        '''
        super(LocInfraPy,self).__init__(conf_file,'loc_params')
        self.assocversion=0
        print('Loc version:',self.assocversion)
        self.year=int(self.general_PARAM['year'])
        self.dayofyearini=int(self.general_PARAM['dayofyearini'])
        self.dayofyearend=int(self.general_PARAM['dayofyearend'])
        self.jdayini = int(self.year*1000 + self.dayofyearini)
        self.jdayend = int(self.year*1000 + self.dayofyearend)


        try:
            self.networkName=self.general_PARAM['networktable']
            self.affiliationName=self.general_PARAM['affiliationtable']
        except :
            print('No specific tables')
            self.networkName=None
            self.affiliationName=None
        try:
            self.siteName=self.general_PARAM['sitetable']
            self.wfdiscName=self.general_PARAM['wfdisctable']
        except :
            print('No specific tables')
            self.siteName=None
            self.wfdiscName=None

        self.net=self.task_PARAM['network']
        self.passocid=int(self.task_PARAM['passocid'])
        self.maxlat=float(self.task_PARAM['maxlat'])
        self.minlat=float(self.task_PARAM['minlat'])
        self.maxlon=float(self.task_PARAM['maxlon'])
        self.minlon=float(self.task_PARAM['minlon'])

        try:
            self.priors=int(self.task_PARAM['priors'])
        except :
            print('No specific tables')
            self.priors=0

        self.lims=[[self.minlat,self.maxlat],[self.minlon,self.maxlon]]

        try:
            self.resultstable=self.task_PARAM['resultstable']
        except :
            print('No specific tables')
            self.resultstable=None
        try:
            print('Using multicores')
            self.numcores=int(self.task_PARAM['numcores'])
        except :
            print('Using one core only')
            self.numcores=1

    def database_connecting(self):
        print('connecting')

        session,tables = load_config(self.db_PARAM)
        self.session=session
        self.Site=tables['site']
        self.Wfdisc=tables['wfdisc']


        import pisces.schema.css3 as kba
        class Site(kba.Site):
            __tablename__ = 'site'

        '''
        class Network(kba.Network):
            __tablename__ = 'network'
        '''

        class Affiliation(kba.Affiliation):
            __tablename__ = 'affiliation'

        class Fk_results(schema.fk_results):
            __tablename__ = 'Fk_results'

        class Fk_params(schema.fk_params):
            __tablename__ = 'Fk_params'

        class Fd_params(schema.fd_params):
            __tablename__ = 'Fd_params'

        class AssocParams(schema.assoc_params):
            __tablename__ = 'Assoc_params'

        class AssocResults(schema.assoc_results):
            __tablename__ = 'Assoc_results'

        class LocParams(schema.loc_params):
            __tablename__ = 'loc_params'

        class LocResults(schema.loc_results):
            __tablename__ = 'loc_results'

        print('tables done')
        '''
        for  fdi in self.fdtables:
            exec('class FD'+str(numT)+'(schema.fd_results):\n\t__tablename__ = \''+str(fdi)+'\'')
            exec('self.FD'+str(numT)+'='+'FD'+str(numT))
            numT=numT+1
        '''
        '''
        self.dict_namefk={}

        self.fdtables=[]
        for  fdi in self.fdtables_names:
            self.fdtables.append(type(str(fdi),(schema.fd_results,),{'__tablename__':str(fdi)}))
        '''

        self.Site=Site
        #self.Network=Network
        self.Affiliation=Affiliation
        self.FK_par=Fk_params
        self.FK_results=Fk_results
        self.Fd_par=Fd_params
        #self.Fd_results=Fd_results

        self.AssocPar=AssocParams
        self.AssocResults=AssocResults

        self.LocPar=loc_params
        self.LocResults=loc_results

        self.LocPar.__table__.create(self.session.bind,checkfirst=True)
        self.LocResults.__table__.create(self.session.bind,checkfirst=True)


        self.db_connected=True
        return self.db_connected

        try:
            self.Ploc_Q=self.session.query(self.LocPar). \
                filter(self.LocPar.maxlat==self.maxlat).\
                filter(self.LocPar.minlat==self.minlat). \
                filter(self.LocPar.maxlon==self.maxlon).\
                filter(self.LocPar.minlon==self.minlon).\
                filter(self.LocPar.priors==self.priors).\
                filter(self.LocPar.algorithm==self.algorithm).\
                all()

            if len(self.Ploc_Q)>1:
                print('issue with the database too many parameters entries, there should be just one')
                embed()
            if len(self.Ploc_Q)==1:
                self.Ploc_Q=self.Ploc_Q[0]

        except Exception as x1:
            print("issue with the table or first assoc entered")
            print(x1)
            embed()
            self.Ploc_Q=[]

        if bool(self.Ploc_Q)==False:
            print('New process parameters, write process to INFRA_ASSOC_PARAM table')
            new_row=self.session.query(self.LocPar).count()
            try:
                res=self.LocPar(  maxlat=self.maxlat,\
                                     minlat=self.minlat,\
                                     maxlon=self.maxlon,\
                                     minlon=self.minlon,\
                                     priors=self.priors,\
                                     algorithm=self.algorithm,\
                                     plocid=new_row)
            except Exception as x1:
                print('problem writing to the loc param file')
                print("Unexpected error:", x1)
                embed()
            #embed()
            self.session.add(res)
            self.session.commit()
            self.Ploc_Q=self.session.query(self.LocPar). \
                filter(self.LocPar.maxlat==self.maxlat).\
                filter(self.LocPar.minlat==self.minlat). \
                filter(self.LocPar.maxlon==self.maxlon).\
                filter(self.LocPar.minlon==self.minlon).\
                filter(self.LocPar.priors==self.priors).\
                filter(self.LocPar.algorithm==self.algorithm).\
                one()
            self.plocid=self.Ploc_Q.plocid
            #embed()
        else:
            print('process already in table: Assoc params table')
            self.plocid=self.Ploc_Q.plocid
            print(self.Ploc_Q)


        self.tables_FD={}  # all tables for detection
        self.tables_FK={}  # all tables for fk
        self.db_connected=True

        self.lims_bisl=[[self.Ploc_Q.minlat,self.Ploc_Q.maxlat],[self.Ploc_Q.minlon,self.Ploc_Q.maxlon]]

        return self.db_connected

    def data_retrieving(self):
        '''

        '''

    def data_retrieving(self,dayofyear):
        '''

        '''
        jday = int(self.year*1000 + dayofyear)
        t0 = UTCDateTime( year=self.year, julday=dayofyear).timestamp
        te = UTCDateTime( year=self.year, julday=dayofyear+1).timestamp

        # this should be revised depending on the maximum distance of the limits.
        t0_E=t0
        te_E=te
        self.t0_E=t0
        self.te_E=te
        self.time_initial=t0
        self.time_end=te
        t0_E_jday=UTCDateTime(t0_E).year*1000+UTCDateTime(t0_E).julday
        te_E_jday=UTCDateTime(te_E).year*1000+UTCDateTime(te_E).julday


        self.te_E_jday=te_E_jday
        self.t0_E_jday=t0_E_jday

        try:
            #self.Network_Q=self.session.query(self.Network).filter(self.Network.net==self.net).all()
            self.Affiliation_Q=self.session.query(self.Affiliation).filter(self.Affiliation.net==self.net).all()
        except Exception as ex1:
            print('Error with network retrieving', ex1)
            #embed()
            exit(0)
        #Query_temp=self.session.query(self.AssocResults).filter(self.AssocResults.net==self.net).filter(self.AssocResults.passocid==self.passocid).filter(self.AssocResults.timeini>=t0_E).filter(self.AssocResults.timeini<=te_E).all()

        id_res=self.session.query(self.LocResults)\
                          .filter(self.LocResults.net==self.net)\
                          .filter(self.LocResults.eventid==-1)\
                          .filter(self.LocResults.plocid==self.plocid)\
                          .filter(self.LocResults.timeini==self.time_initial)\
                          .filter(self.LocResults.timeend==self.time_end)\
                          .all()
        if len(id_res)==1:
            return 1
        else:
            print('start getting data for analysis')

        refsta=self.meanlocations()

        Query_temp_evID=self.session.query(self.AssocResults.eventid).filter(self.AssocResults.net==self.net).filter(self.AssocResults.passocid==self.passocid).filter(self.AssocResults.timeini>=t0_E).filter(self.AssocResults.timeini<=te_E).order_by(self.AssocResults.eventid).distinct().all()
        event_ID=[]
        eventDET_ID=[]
        eventDET_TABLE=[]
        eventDET_TABLE_ALL=[]
        eventPRIORS_FILE=[]
        for qq in Query_temp_evID:
            #Query_temp=self.session.query(self.AssocResults).filter(self.AssocResults.eventid==qq[0]).filter(self.AssocResults.qassoc>0)
            Query_temp_D=self.session.query(self.AssocResults).filter(self.AssocResults.eventid==qq[0]).filter(self.AssocResults.qassoc>0).order_by(self.AssocResults.fdid).all()
            if len(Query_temp_D)==0:
                #print 'no results'
                continue

            if not (self.priors==0):
                ss=Query_temp_D[0]
                HH=UTCDateTime((ss.timeini+ss.timeend)/2).hour
                if HH<3 or HH >=21:
                    HH_E=0
                elif HH>=3 and HH<9:
                    HH_E=6
                elif HH>=9 and HH<15:
                    HH_E=12
                elif HH>=15 and HH<21:
                    HH_E=18
                else:
                    print('problem with hour')
                    exit()
                package_directory = os.path.dirname(os.path.abspath(__file__))
                MMS='%02d'%UTCDateTime((ss.timeini+ss.timeend)/2).month
                HH_ES='%02d'%HH_E
                prior_f=glob.glob(package_directory+'/assoc_lib/priors/*_'+MMS+'_'+HH_ES+'00UTC.pri')[0]
                #priors = prop_stats()
                #priors.load(prior_f)
            else:
                prior_f=None

            tempID=[]
            tempTABLE=[]

            for qt in Query_temp_D:
                tempID.append(qt.fdid)
                tempTABLE.append(qt.fdtable)
            if (tempID in eventDET_ID):
                continue

            event_ID.append(qq[0])
            eventDET_ID.append(tempID)
            eventDET_TABLE.append(tempTABLE)
            eventDET_TABLE_ALL.extend(tempTABLE)
            eventPRIORS_FILE.append(prior_f)

        self.det_GEN=[]
        self.det_EVID=[]
        self.det_PRIORSFILE=[]

        for tt in range(len(eventDET_ID)):
            id_res=self.session.query(self.LocResults)\
                                      .filter(self.LocResults.plocid==self.plocid)\
                                      .filter(self.LocResults.net==self.net)\
                                      .filter(self.LocResults.eventid==event_ID[tt])\
                                      .all()

            if bool(id_res)==True:
                sys.stdout.write('\r'+'Event already analyzed and stored:'+str(tt)+'  '+str(len(eventDET_ID))+'  '+str(event_ID[tt]))
                sys.stdout.flush()
                #print 'Event already analyzed and stored:',tt,len(eventDET_ID),event_ID[tt]
                #embed()
                continue


            newtables=set(eventDET_TABLE[tt])

            for tname in newtables:
                if not(tname in list(self.tables_FD.keys())):
                    self.tables_FD[tname]=type(tname.encode('utf-8') ,(schema.fd_results,),{'__tablename__':tname.encode('utf-8')})

            sta_nam=[]
            det_Q=[]
            for ww in range(len(eventDET_ID[tt])):
                fdid=(eventDET_ID[tt][ww])
                ftable=(eventDET_TABLE[tt][ww])
                res_aux=self.session.query(self.tables_FD[ftable]).filter(self.tables_FD[ftable].fdid==fdid).one()
                det_Q.append(res_aux)
                sta_nam.append(res_aux.sta)

            det_tot=[]

            sys.stdout.write('\r'+str(tt)+'  '+str(len(eventDET_ID))+'   '+str(event_ID[tt]))
            sys.stdout.flush()

            if len(set(sta_nam))<2:
                sys.stdout.write('\r'+'Assoc with one detection')
                sys.stdout.flush()
            else:

                for ww in range(len(eventDET_ID[tt])):
                    fdid=(eventDET_ID[tt][ww])
                    ftable=(eventDET_TABLE[tt][ww])
                    #Query_temp_D=self.session.query(self.tables_FD[ftable]).filter(self.tables_FD[ftable].fdid==fdid).one()
                    Query_temp_D=det_Q[ww]
                    fk_tabname=Query_temp_D.fktablename
                    if not( fk_tabname in list(self.tables_FK.keys())):
                        self.tables_FK[fk_tabname]=type(fk_tabname.encode('utf-8') ,(schema.fk_results,),{'__tablename__':fk_tabname.encode('utf-8')})

                    #remember that the locaiton is looking for the point with maximum Fval

                    #Query_temp_FK=self.session.query(self.tables_FK[fk_tabname]).filter(self.tables_FK[fk_tabname].timeini==Query_temp_D.timeini).one()
                    Query_fktempMAX=self.session.query(self.tables_FK[fk_tabname]).filter(self.tables_FK[fk_tabname].timeini==Query_temp_D.maxf_time).one()
                    for refs in refsta:
                        if refs['name']==Query_fktempMAX.sta:
                            aai=refs
                            break
                    det1 = inf_det_global(aai['lat'], aai['lon'], UTCDateTime(Query_fktempMAX.timeini).datetime, Query_fktempMAX.bz,  Query_fktempMAX.fval, aai['numsta'],DetID=aai['name'])
                    det_tot.append(det1)

            self.det_GEN.append(det_tot)
            self.det_EVID.append(event_ID[tt])
            self.det_PRIORSFILE.append(eventPRIORS_FILE[tt])
            #print tt,len(eventDET_ID),event_ID[tt]
        return 2

    def data_processingLOC(self):

        print('data processing Effective',short_time(UTCDateTime(self.time_initial)),short_time(UTCDateTime(self.time_end)))
        if len(self.det_GEN)<1:
            print('this may be a case of a day with no detections')
            return 0
        args_bisl = [0.0] * len(self.det_GEN)
        args_bisl_effect=[]
        index_no_empty=[]
        out= [0.0] * len(self.det_GEN)
        for det_ii in range(len(self.det_GEN)):
            out[det_ii]=[]
            try:
                if not(self.priors==0):
                    priorsL = prop_stats()
                    priorsL.load(self.det_PRIORSFILE[det_ii])
                else:
                    priorsL=None
            except Exception as ex1:
                print('problem loading priors')
                embed()
            #[priorsL] i sjust so it doestno complaing aboptu the len
            args_bisl[det_ii]=(self.det_GEN[det_ii], self.lims,priorsL,50, None, self.numcores,False,False, None,det_ii,len(self.det_GEN))
            valQC=bisl_qc_check(self.det_GEN[det_ii], self.lims)
            if (len(self.det_GEN[det_ii])>0) and (valQC==True):
                index_no_empty.append(det_ii)
                args_bisl_effect.append(args_bisl[det_ii])


        print('tot',len(self.det_GEN),'  effective',len(args_bisl_effect))
        if self.numcores>1:
            print('Available cores :', multiprocessing.cpu_count())
            print('Cores to be used:',self.numcores)
            if self.numcores >= multiprocessing.cpu_count():
                print('There is only:',multiprocessing.cpu_count(), 'available')
                print('infrapy will be using only', multiprocessing.cpu_count()-1)
                self.numcores=multiprocessing.cpu_count()-1
            pool = Pool(processes=self.numcores)
            try:
                out_effec=pool.map(run_bislSC_wrapper,args_bisl_effect)
            except Exception as ex1:
                    pool.terminate()
                    print('Error:',ex1)
                    embed()
                    sys.exit(0)
            except KeyboardInterrupt:
                #embed()
                pool.terminate()
                print('Do not finish processing and discard values')
                sys.exit(0)
            pool.close()

            for ii_noe in range(len(index_no_empty)):
                ii_gen=index_no_empty[ii_noe]
                out[ii_gen]= out_effec[ii_noe]


        else:
            print('multicore capabilities disabled (standard one core calculation)')
            for aar in range(len(args_bisl)):
                if bisl_qc_check(self.det_GEN[aar], self.lims)==False:
                    print('Assoc may not converge')
                    val=[]
                else:
                    val=run_bislSC(*args_bisl[aar])
                out.append(val)

        for det_ii in range(len(self.det_GEN)):
            #print det_ii,' from',len(self.det_GEN)
            stations=[]
            for dd in self.det_GEN[det_ii]:
                stations.append(dd.ID)
            sta_sets=set(stations)

            id_res=self.session.query(self.LocResults)\
                                      .filter(self.LocResults.plocid==self.plocid)\
                                      .filter(self.LocResults.net==self.net)\
                                      .filter(self.LocResults.eventid==self.det_EVID[det_ii])\
                                      .all()

            if bool(id_res)==True:
                print('Event already analyzed and stored',self.det_EVID[det_ii])
                continue
            result=out[det_ii]
            lastEVENTIDQ=self.session.query(func.max(self.LocResults.locid))
            lastEVENTID=lastEVENTIDQ[0][0]
            if lastEVENTID is None:
                lastEVENTID=0


            id_resC=self.session.query(self.LocResults).count()+1

            if bool(id_res)==False:
                if (len(sta_sets)>=2) & (len(result)>0):
                    res=self.LocResults(\
                                          locid=id_resC,\
                                          plocid=self.plocid,\
                                          net=self.net,\
                                          eventid=self.det_EVID[det_ii],\
                                          timeorigmap =UTCDateTime(result['t_MaP']),\
                                          timeorigmean =UTCDateTime(result['t_mean']),\
                                          timeorigmin=UTCDateTime(result['t_min']),\
                                          timeorigmax=UTCDateTime(result['t_max']),\
                                          timeorigvar=UTCDateTime(result['t_var']),\
                                          latorigmean=result['lat_mean'],\
                                          latorigvar=result['lat_var'],\
                                          lonorigmean=result['lon_mean'],\
                                          lonorigvar=result['lon_var'],\
                                          latlonorigcovar=result['covar'],\
                                          latorigmap=result['lat_MaP'],\
                                          lonorigmap=result['lon_MaP'],\
                                          timeini=self.time_initial,\
                                          timeend=self.time_end,\
                                          numstations=len(sta_sets)\
                                          )
                else:
                    res=self.LocResults(\
                                          locid=id_resC,\
                                          plocid=self.plocid,\
                                          net=self.net,\
                                          eventid=self.det_EVID[det_ii],\
                                          numstations=len(sta_sets)\
                                          )
                try:
                    self.session.add(res)
                    self.session.commit()
                except Exception as ex1:
                    print('this event was already analyzed, rolling back and continue to the next event')
                    self.session.rollback()

        print('\rlocations written', len(self.det_GEN))
        id_resC=self.session.query(self.LocResults).count()+1
        lastEVENTIDQ=self.session.query(func.max(self.LocResults.locid))
        lastEVENTID=lastEVENTIDQ[0][0]
        if lastEVENTID is None:
            lastEVENTID=0

        res=self.LocResults(locid=id_resC,\
                eventid=-1,\
                plocid=self.plocid,\
                timeini=self.time_initial,\
                timeend=self.time_end,\
                net=self.net,\
                numstations=-1)
        self.session.add(res)
        self.session.commit()


    def data_processing(self):
        '''
        Constructor
        '''

        try:
            self.Affiliation_Q=self.session.query(self.Affiliation).filter(self.Affiliation.net==self.net).all()
        except Exception as ex1:
            print('Error with network retrieving', ex1)
            exit(0)

        for day in np.arange(self.dayofyearini,self.dayofyearend):
            print('year:',self.year,'day:', day)
            self.dayofyear=day
            retrieve_ret=self.data_retrieving(day)

            if retrieve_ret==2:
                self.data_processingLOC()
            elif retrieve_ret==1 :
                print('day already analyzed')
            else:
                print('there was an error or not data')

        #return self.data_processed


    def meanlocations(self):
        refSTA=[]
        for aai in self.Affiliation_Q:
            try:
                STA_dataM=self.session.query(self.Site).filter(self.Site.sta==aai.sta).filter(((self.Site.offdate==-1)|(self.Site.offdate>self.te_E_jday)) & (self.Site.ondate<self.t0_E_jday)).one()
            except Exception as ex1:
                #print
                print('there is more than just one station:', aai.sta,'  ',ex1)
                embed()
                exit()
            #embed()
            refSTA.append(STA_dataM.refsta)

        refstations_l=list(set(refSTA))
        refsta=[]
        for aai in refstations_l:
            STA_dataM=self.session.query(self.Site).filter(self.Site.refsta==str(aai)).filter(((self.Site.offdate>self.te_E_jday) |(self.Site.offdate==-1))& (self.Site.ondate<self.t0_E_jday)).all()
            array_lo=[]
            array_la=[]
            array_el=[]
            for sta_i in STA_dataM:
                array_la.append(sta_i.lat)
                array_lo.append(sta_i.lon)
                array_el.append(sta_i.elev)
            array_la=np.asarray(array_la)
            array_lo=np.asarray(array_lo)
            array_el=np.asarray(array_el)
            refsta.append({'lon':np.mean(array_lo),'lat':np.mean(array_la),'elev':np.mean(array_el),'name':aai,'numsta':len(array_la)})
        return refsta


if __name__ == '__main__':
    pdetect=AssocInfraPy_LANL('../conf_files/InfraConfig_Assoc')
    pdetect.database_connecting()
    pdetect.data_processing()
