'''
Created on Oct 31, 2014
Updated Jan 2020

@author: omarcillo, fkdd
'''

from .base import Base

import sys, pdb

import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import func
from sqlalchemy import MetaData

import pisces as ps
from pisces.util import load_config
from pisces.io.trace import read_waveform
from obspy.core import UTCDateTime
from obspy.core.util import AttribDict
from datetime import datetime

from scipy import stats
import numpy as np
import scipy as sc
from IPython import embed
import pisces as ps

from sqlalchemy import func
import pylab as py

import matplotlib.pyplot as pl
import matplotlib.mlab as mpy

from ...utils.cart2pol import cart2pol
from ...utils.short_time import short_time
from .. import schema

import matplotlib.dates as mdates
from ...propagation import likelihoods
from ...propagation.likelihoods import InfrasoundDetection
from ...propagation import infrasound
from ...association import hjl
import numpy as np

import pathos.multiprocessing as mp
from multiprocessing import cpu_count

from infrapy.association import hjl
from infrapy.utils import data_io

class AssocInfraPy_LANL(Base):
    '''
    classdocs
    '''

    algorithm='Blom and Euler'

    def __init__(self, conf_file=[]):
        '''
        Constructor
        '''
        super(AssocInfraPy_LANL,self).__init__(conf_file,'AssocLocParams')
        self.assocversion=0
        print('Assoc version:',self.assocversion)
        self.year=int(self.general_PARAM['year'])
        self.dayofyearini=int(self.general_PARAM['dayofyearini'])
        self.dayofyearend=int(self.general_PARAM['dayofyearend'])
        self.jdayini = int(self.year*1000 + self.dayofyearini)
        self.jdayend = int(self.year*1000 + self.dayofyearend)
        self.cpu=self.general_PARAM['cpucnt']
        self.pl = mp.ProcessingPool(cpu_count() - 1)


        self.net=self.task_PARAM['network']
        self.pfdid=self.task_PARAM['pfdetectid']
        self.pfkid=self.task_PARAM['pfkid']
        self.beamwidth=float(self.task_PARAM['beamwidth'])
        self.rangemax=float(self.task_PARAM['rangemax'])
        self.distmax=float(self.task_PARAM['distmax'])
        self.clusterthresh=float(self.task_PARAM['clusterthresh'])
        self.trimthresh=(self.task_PARAM['trimthresh'])
        self.eventdetmin=float(self.task_PARAM['eventdetmin'])
        self.eventarrmin=float(self.task_PARAM['eventarrmin'])
        self.duration=float(self.task_PARAM['duration'])


        try:
            self.resultstable=self.task_PARAM['resultstable']
        except :
            print('No specific tables')
            self.resultstable=None

        listK=list(self.task_PARAM.keys())
        self.fdtables_names=[]
        for li in listK:
            if bool('fdtable_' in li):
                self.fdtables_names.append(self.task_PARAM[li])
        if len( self.fdtables_names)==0:
            print('NO tables with fd results were included, define one fd results tables, or include Fd_results (this is the table where results are written by default, but needs to be specified ) ')
            sys.exit()
        self.num_tables=len(self.fdtables_names)
            #self.fdtables.append('Fd_results')
        #embed()

    def database_connecting(self):
        print('connecting')

        session,tables = load_config(self.db_PARAM)
        self.session=session
        self.Site=tables['site']
        self.Wfdisc=tables['wfdisc']
        self.Affiliation=tables['affiliation']
        import pisces.schema.css3 as kba

        class FK_results(schema.fk_results):
            __tablename__ = 'FK_results'

        class FK_params(schema.fk_params):
            __tablename__ = 'FK_params'

        class FD_params(schema.fd_params):
            __tablename__ = 'FD_params'

        self.dict_namefk={}


        self.fdtables=[]
        for  fdi in self.fdtables_names:
            self.fdtables.append(type(str(fdi),(schema.fd_results,),{'__tablename__':str(fdi)}))


        class ASSOC_params(schema.ASSOC_params):
            __tablename__ = 'ASSOC_params'

        class ASSOC_results(schema.ASSOC_results):
            __tablename__= self.resultstable


        self.FK_par=FK_params
        self.FK_results=FK_results
        self.FD_par=FD_params
        #self.Fd_results=Fd_results

        self.ASSOC_par=ASSOC_params
        self.ASSOC_results=ASSOC_results

        self.ASSOC_par.__table__.create(self.session.bind,checkfirst=True)
        self.ASSOC_results.__table__.create(self.session.bind,checkfirst=True)


        try:
            self.Passoc_Q=self.session.query(self.ASSOC_par). \
                filter(self.ASSOC_par.beamwidth==self.beamwidth).\
                filter(self.ASSOC_par.rangemax==self.rangemax). \
                filter(self.ASSOC_par.clusterthresh==self.clusterthresh).\
                filter(self.ASSOC_par.trimthresh==self.trimthresh).\
                filter(self.ASSOC_par.eventdetmin==self.eventdetmin).\
                filter(self.ASSOC_par.eventarrmin==self.eventarrmin).\
                filter(self.ASSOC_par.duration==self.duration).\
                all()


            if len(self.Passoc_Q)>1:
                print('issue with the database too many parameters entries, there should be just one')
                embed()
            if len(self.Passoc_Q)==1:
                self.Passoc_Q=self.Passoc_Q[0]
        except Exception as x1:
            print("issue with the table or first assoc entered")
            print(x1)
            embed()
            self.Passoc_Q=[]
            print(Passoc_Q)
        if bool(self.Passoc_Q)==False:
            print('New process parameters, write process to INFRA_ASSOC_PARAM table')
            new_row=self.session.query(self.ASSOC_par).count()
            try:
                res=self.ASSOC_par(  beamwidth=self.beamwidth,\
                                     rangemax=self.rangemax,\
                                     clusterthresh=self.clusterthresh,\
                                     trimthresh=self.trimthresh,\
                                     eventdetmin=self.eventdetmin,\
                                     algorithm=self.algorithm,\
                                     eventarrmin=self.eventarrmin,\
                                     duration=self.duration,\
                                     passocid=new_row)
            except Exception as x1:
                print('problem writing to the assoc param file')
                print("Unexpected error:", x1)
                embed()

            self.session.add(res)
            self.session.commit()
            self.Passoc_Q=self.session.query(self.ASSOC_par). \
                filter(self.ASSOC_par.beamwidth==self.beamwidth).\
                filter(self.ASSOC_par.rangemax==self.rangemax). \
                filter(self.ASSOC_par.clusterthresh==self.clusterthresh).\
                filter(self.ASSOC_par.trimthresh==self.trimthresh).\
                filter(self.ASSOC_par.eventdetmin==self.eventdetmin).\
                filter(self.ASSOC_par.eventarrmin==self.eventarrmin).\
                filter(self.ASSOC_par.duration==self.duration).\
                one()
            self.passocid=self.Passoc_Q.passocid
            #embed()
        else:
            print('process already in table: Assoc params table')
            self.passocid=self.Passoc_Q.passocid
            print(self.Passoc_Q)
        self.db_connected=True
        return self.db_connected

    def data_retrievingS(self,win_start, win_end):
        '''

        '''

        # this should be revised depending on the maximum distance of the limits.

        self.time_initial=win_start
        self.time_end=win_end

        id_res=self.session.query(self.ASSOC_results).filter(self.ASSOC_results.net==self.net)\
                                  .filter(self.ASSOC_results.fdid ==-1)\
                                  .filter(self.ASSOC_results.passocid ==self.passocid)\
                                  .filter(self.ASSOC_results.timeini==self.time_initial)\
                                  .filter(self.ASSOC_results.timeend==self.time_end)\
                                  .filter(self.ASSOC_results.qassoc==-1).all()
        if len(id_res)==1:
            return 1
        else:
            print('start getting data for analysis')
        self.Detection=[]
        self.Detection_Q=[]
        #embed()
        refSTA=[]
        for aai in self.Affiliation_Q:

            try:
                STA_dataM=self.session.query(self.Site).filter(self.Site.sta==aai.sta).one()

            except Exception as ex1:
                #print
                print('there is more than just one station:', aai.sta,'  ',ex1)
                embed()
                exit()
            #embed()
            refSTA.append(STA_dataM.refsta)

        refstations_l=list(set(refSTA))
        refsta=[]
        #embed()
        for aai in refstations_l:
            STA_dataM=self.session.query(self.Site).filter(self.Site.refsta==str(aai)).all()
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
        self.det_tot=[]
        self.fdtable_name=[]
        #embed()

        for aai in refsta:
            print('getting data from:',aai['name'])
            # here it looks for any detections in all the tables
            for ti in range(self.num_tables):
                try:

                    fd_res=self.session.query(self.fdtables[ti]).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).filter(self.fdtables[ti].timeini>=self.time_initial).all()
                    #times_ini=self.session.query(self.fdtables[ti].timeini).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).filter(self.fdtables[ti].timeini>=self.time_initial).filter(self.fdtables[ti].timeini<=self.time_end).all()
                    times_ini=self.session.query(self.fdtables[ti].timeini).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).filter(self.fdtables[ti].timeini>=self.time_initial).all()
                    if len(fd_res)>0:
                        print('length results:', len(fd_res))
                        for dqi in range(len(fd_res)):
                            if win_start < fd_res[dqi].timeini < win_end:

                                class Fk_results(schema.fk_results):
                                    __tablename__ = fd_res[0].fktablename

                                fk_res=self.session.query(Fk_results).filter(Fk_results.pfkid==fd_res[dqi].pfkid).filter(Fk_results.sta==aai['name']).filter(Fk_results.timeini.between(fd_res[dqi].timeini,fd_res[dqi].timeend)).all()

                                fval=[]
                                bz=[]

                                for fk_i in fk_res:
                                    fval.append(fk_i.fval)
                                    bz.append(fk_i.bz)

                                fval=np.asarray(fval)
                                try:
                                    fvalm_i=np.argmax(fval)

                        #det1 = inf_det_global(aai['lat'], aai['lon'], UTCDateTime(fd_res.timeini).datetime, fd_res.bz,  fd_res.fval, aai['numsta'])
                                    det1 = (aai['lat'], aai['lon'], fd_res[dqi].timeini, bz[fvalm_i],  fval[fvalm_i], aai['numsta'], fd_res[dqi].fdid,aai['name'],self.fdtables_names[ti])

                                    self.det_tot.append(det1)
                                    self.fdtable_name.append(self.fdtables_names[ti])
                                except Exception as ex1:
                                    pass
                        embed()
                except Exception as x1:
                    print('There is an error',x1)
                    embed()
                    exit()
        print('all data retrieved')
        return 2


    def data_processingASSOC(self):
        '''

        '''
        print('data processing',short_time(UTCDateTime(self.time_initial)),short_time(UTCDateTime(self.time_end)))

        det_list = data_io.db2dets(self.det_tot)
        EVIDs=[]
        embed()
        if len(det_list)>1:
            try:
                #EVIDs,DAQ,CAQ=assoc(det_list, self.lims, float(self.assocthresh), show_result=False,parallel=True,num_cores=self.numcores)
                #labels, dists = hjl.run(det_list, self.clusterthresh, dist_max=self.distmax, bm_width=self.beamwidth, rng_max=self.rangemax, trimming_thresh=self.trimthresh, pool=self.pl)
                labels, dists = hjl.run(det_list, self.clusterthresh, dist_max=self.distmax, bm_width=self.beamwidth, rng_max=self.rangemax,  pool=self.pl)
                clusters, qualities = hjl.summarize_clusters(labels, dists)

                for n in range(len(clusters)):
                    print("Cluster:", clusters[n], '\t', "Cluster Quality:", 10.0**(-qualities[n]))
                    lastEVENTIDQ=self.session.query(func.max(self.ASSOC_results.eventid)).all()
                    lastEVENTID=lastEVENTIDQ[0][0]
                    if lastEVENTID is None:
                        lastEVENTID=int(0)
                    lastEVENTID=lastEVENTID + 1
                    #embed()
                    for nn in range(len(clusters[n])):
                        det_id = clusters[n][nn]
                        id_res=self.session.query(self.ASSOC_results).filter(self.ASSOC_results.net==self.net)\
                                                  .filter(self.ASSOC_results.fdid ==self.det_tot[det_id][6])\
                                                  .filter(self.ASSOC_results.passocid ==self.passocid)\
                                                  .filter(self.ASSOC_results.timeini==self.time_initial)\
                                                  .filter(self.ASSOC_results.timeend==self.time_end)\
                                                  .filter(self.ASSOC_results.qdetcluster==10.0**(-qualities[n]))\
                                                  .filter(self.ASSOC_results.fdtable==self.det_tot[det_id][8])\
                                                  .filter(self.ASSOC_results.sta==self.det_tot[det_id][7]).all()
                        id_resC=self.session.query(self.ASSOC_results).count()+1

                        if bool(id_res)==False:

                            res=self.ASSOC_results(associd=id_resC,\
                                    fdid=self.det_tot[det_id][6],\
                                    eventid=int(lastEVENTID),\
                                    passocid=self.passocid,\
                                    net=self.net,\
                                    timeini=self.time_initial,\
                                    timeend=self.time_end,\
                                    qdetcluster=10.0**(-qualities[n]),\
                                    fdtable=self.det_tot[det_id][8],\
                                    sta=self.det_tot[det_id][7])
                            self.session.add(res)
                            self.session.commit()
                print('associations written', len(clusters))

            except Exception as ex1:
                print('error running assoc:',ex1)
                embed()
                exit()

            '''

            for ii in range(len(EVIDsS)):
                EVIDsN.append({'ID':self.det_tot[sorted_index[ii]].ID,'eventID':EVIDsS[ii],'qdetcluster':DAQS[ii],'qassoc':CAQS[ii],'fdtable':self.fdtable_nameS[ii]})
            lastEVENTIDQ=self.session.query(func.max(self.ASSOC_results.eventid)).all()
            lastEVENTID=lastEVENTIDQ[0][0]
            if lastEVENTID is None:
                lastEVENTID=int(0)

            for ev1 in EVIDsN:
                id_res=self.session.query(self.ASSOC_results).filter(self.ASSOC_results.net==self.net)\
                                          .filter(self.ASSOC_results.fdid ==ev1['ID'])\
                                          .filter(self.ASSOC_results.passocid ==self.passocid)\
                                          .filter(self.ASSOC_results.timeini==self.time_initial)\
                                          .filter(self.ASSOC_results.timeend==self.time_end)\
                                          .filter(self.ASSOC_results.qdetcluster==ev1['qdetcluster'])\
                                          .filter(self.ASSOC_results.fdtable==ev1['fdtable'])\
                                          .filter(self.ASSOC_results.qassoc==ev1['qassoc']).all()
                id_resC=self.session.query(self.ASSOC_results).count()+1

                if bool(id_res)==False:

                    res=self.ASSOC_results(associd=id_resC,\
                            fdid=ev1['ID'],\
                            eventid=int(lastEVENTID+ev1['eventID']),\
                            passocid=self.passocid,\
                            net=self.net,\
                            timeini=self.time_initial,\
                            timeend=self.time_end,\
                            qdetcluster=ev1['qdetcluster'],\
                            fdtable=ev1['fdtable'],\
                            qassoc=ev1['qassoc'])
                    self.session.add(res)
                    self.session.commit()
        print('associations written', len(EVIDs))


        id_resC=self.session.query(self.ASSOC_results).count()+1
        #embed()
        lastEVENTIDQ=self.session.query(func.max(self.ASSOC_results.eventid)).all()
        lastEVENTID=(lastEVENTIDQ[0][0])

        if lastEVENTID is None:
            lastEVENTID=int(0)

        res=self.ASSOC_results(associd=id_resC,\
                passocid=self.passocid,\
                net=self.net,\
                fdtable=self.net,\
                timeini=self.time_initial,\
                timeend=self.time_end,\
                #eventid=lastEVENTID+bytes(1))
                eventid=lastEVENTID)
        self.session.add(res)
        self.session.commit()

'''

    def data_processing(self):
        '''
        Constructor
        '''
        #import pdb; pdb.set_trace()
        try:
            self.Affiliation_Q=self.session.query(self.Affiliation).filter(self.Affiliation.net==self.net).all()
        except Exception as ex1:
            print('Error with network retrieving', ex1)
            exit(0)
            # define  the duration, maximum propagation time assuming 220 m/s minimum
            # celerity, and source window length in minutes
        t_start = UTCDateTime(year=self.year, julday = self.dayofyearini, hour=0, minute=0)
        t_start = np.datetime64(t_start,'s')
        t_end = UTCDateTime(year=self.year, julday = self.dayofyearend, hour=23, minute=59)
        t_end = np.datetime64(t_end,'s')
        days = self.jdayend-self.jdayini
        days = days * 24 * 60
        duration_dd = int((t_end - t_start).astype('m8[s]').astype(float) / 60.0)
        max_prop_tm = int((self.rangemax / 0.22) / self.duration)
        src_win = int(max_prop_tm * 0.5)



        for day in np.arange(self.dayofyearini,self.dayofyearend):

            print('year:',self.year,'day:', day)
            self.dayofyear=day
            for dt in range(0, duration_dd, int(src_win)):
                win_start = t_start + np.timedelta64(dt, 'm')
                win_end = t_start + np.timedelta64(dt + int(src_win + max_prop_tm), 'm')
                win_start = UTCDateTime(win_start.astype(datetime)).timestamp
                win_end = UTCDateTime(win_end.astype(datetime)).timestamp
                print('Window Start',UTCDateTime(win_start))
                print('Window End',UTCDateTime(win_end))

                retrieve_ret=self.data_retrievingS(win_start,win_end)
                for aa in range(len(self.det_tot)):
                    print(UTCDateTime(self.det_tot[aa][2]))

                # retrieve time based on UTC date time string
                if retrieve_ret==2:
                    self.data_processingASSOC()
                elif retrieve_ret==1 :
                    print('time period already analyzed')
                else:
                    print('there was an error or not data')


        #return self.data_processed


if __name__ == '__main__':
    pdetect=AssocInfraPy_LANL('../conf_files/InfraConfig_Assoc')
    pdetect.database_connecting()
    pdetect.data_processing()
