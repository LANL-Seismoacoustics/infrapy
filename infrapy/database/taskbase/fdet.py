'''
Updated May 2020
@fkdd @omarcillo
'''

from .base import Base
from ...utils.cart2pol import cart2pol
from ...utils.short_time import short_time
from .. import schema
from ...detection import beamforming_new

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
from obspy.core import trace
from obspy.core import Stream
from obspy.core.util import AttribDict
from obspy.signal.array_analysis import *
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

import timeit
import time

#matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import matplotlib.mlab as mpy

pl.ioff()

import time

import warnings
warnings.filterwarnings("ignore")



class FDet(Base):

    '''
    classdocs
    '''
    def __init__(self, conf_file=[],ARRAY_NAME=None):
        '''
        Constructor
        '''

        super(FDet,self).__init__(conf_file,'FDetectParams')

        self.year=int(self.general_PARAM['year'])
        self.dayofyearini=int(self.general_PARAM['dayofyearini'])
        self.dayofyearend=int(self.general_PARAM['dayofyearend'])
        if ARRAY_NAME==None:
            self.sta=self.general_PARAM['station']
        else:
            print('Warning')
            print('Station name has been overloaded to:',ARRAY_NAME)
            self.sta=ARRAY_NAME

        try:
            self.chan=self.general_PARAM['channel']
        except:
            print('channel set by default to BDF')
            self.chan='BDF'
        try:
            self.cpucnt=self.general_PARAM['cpucnt']
        except:
            print('cpu count set by default to 8')
            self.cpucnt=8
        try:
            self.domain=self.general_PARAM['domain']
        except:
            print('domain set by default to frequency')
            self.domain='freq'
        self.jdayini = int(self.year*1000 + self.dayofyearini)
        self.jdayend = int(self.year*1000 + self.dayofyearend)
        #self.database=self.general_PARAM['database']

        try:
            self.siteName=self.general_PARAM['sitetable']
            self.wfdiscName=self.general_PARAM['wfdisctable']
        except :
            print('No specific tables')
            self.siteName=None
            self.wfdiscName=None

        self.pthreshold=float(self.task_PARAM['pthreshold'])
        self.detwinlen=self.task_PARAM['detwinlen']
	self.adapwlen=self.task_PARAM['adaptivewlen']
        #self.pthr=float(self.task_PARAM['pthreshold'])
        self.cthr=float(self.task_PARAM['corrthreshold'])
        self.minlen=int(self.task_PARAM['mineventlength'])
        self.dsegmin=float(self.task_PARAM['dsegmin'])
        self.detmethod=(self.task_PARAM['detmethod'])
        self.tb_prod=float(self.task_PARAM['tb_prod'])
        self.backazlim=float(self.task_PARAM['backazlim'])

        try:
            self.fdresults=self.task_PARAM['fdresults']
        except KeyError:
            self.fdresults='FD_RESULTS'

        if self.fdresults=='None':
            print('Using default table Fd_results')
            self.fdresults='FD_RESULTS'

        try:
            self.fkresults=self.task_PARAM['fkresults']
        except KeyError:
            self.fkresults='FK_RESULTS'

        if self.fkresults=='None':
            print('Using default table Fk_results')
            self.fkresults='FK_RESULTS'

        self.pfkid=self.task_PARAM['pfkid']

        try:
            self.numsources=int(self.task_PARAM['numsources'])
        except KeyError:
            self.numsources=1


    def database_connecting(self):

        session,tables = load_config(self.db_PARAM)
        self.session=session
        self.Site=tables['site']
        self.Wfdisc=tables['wfdisc']

        try:
            class Fk_results(schema.fk_results):
                __tablename__ = self.fkresults
        except Exception as Ex1:
            print('fk_res table already defined')
        try:
            class Fk_params(schema.fk_params):
                __tablename__ = 'FK_PARAMS'
        except Exception as Ex1:
            print('FK_PARAMS table already defined')
        try:
            class Fd_results(schema.fd_results):
                __tablename__ = self.fdresults
        except Exception as Ex1:
            print('fd_res table already defined')
        try:
            class Fd_params(schema.fd_params):
                __tablename__ = 'FD_PARAMS'
        except Exception as Ex1:
            print('FD_PARAMS table already defined')
        try:
            self.FK_par=Fk_params
            self.FK_results=Fk_results

            self.FD_par=Fd_params
            self.FD_par.__table__.create(self.session.bind,checkfirst=True)

            self.FD_results=Fd_results
            self.FD_results.__table__.create(self.session.bind,checkfirst=True)

        except Exception as Ex1:
            print(Ex1)
            import warnings
            warnings.filterwarnings("ignore")
            print("table creation issue line 653")
            embed()
            self.db_connected=False
            return self.db_connected
        try:
            self.PFDetect_Q=self.session.query(self.FD_par). \
                filter(self.FD_par.detwinlen==self.detwinlen).\
                filter(self.FD_par.pthreshold==self.pthreshold). \
                filter(self.FD_par.cthr==self.cthr).\
                filter(self.FD_par.minlen==self.minlen).\
                all()

            if len(self.PFDetect_Q)>1:
                print('issue with the database too many parameters entries, there should be just one')
                embed()
            if len(self.PFDetect_Q)==1:
                self.PFDetect_Q=self.PFDetect_Q[0]
        except Exception as ex1:
            print("issue with the table",ex1)
            self.PFDetect_Q=[]
        if bool(self.PFDetect_Q)==False:
            print('New process parameters, write process to INFRA_DETECT_PARAM table')
            new_row=self.session.query(self.FD_par).count()
            res=self.FD_par(detwinlen=self.detwinlen,\
                                pthreshold=self.pthreshold,\
                                cthr=self.cthr,\
                                minlen=self.minlen,\
                                pfdid=new_row)
            self.session.add(res)
            self.session.commit()
            self.PFDetect_Q=self.session.query(self.FD_par).filter(self.FD_par.pfdid==new_row).one()
        else:
            print('process already in table: INFRA_DETECT_PARAM table')
            #embed()
            print(self.PFDetect_Q[:])
        self.db_connected=True
        return self.db_connected

    def data_retrieving(self,dayofyear):
        '''


        '''
        jday = int(self.year*1000 + dayofyear)
        t0 = UTCDateTime( year=self.year, julday=dayofyear).timestamp
        te = t0 + 86400
        self.time_initial=t0
        self.time_end=te
        # the window previous should be divided by 2 as the adaptive win is centered
        self.prewindowlen=int(float(self.detwinlen)/2)
        self.fkpar=self.session.query(self.FK_par).filter(self.FK_par.pfkid==self.pfkid).one()
        # this checks if data was processed
        id_res=self.session.query(self.FD_results).filter(self.FD_results.pfdid==self.PFDetect_Q.pfdid)\
                                  .filter(self.FD_results.pfkid==self.fkpar.pfkid)\
                                  .filter(self.FD_results.timeini==self.time_initial)\
                                  .filter(self.FD_results.timeend==self.time_end)\
                                  .filter(self.FD_results.maxfc==-1)\
                                  .filter(self.FD_results.fktablename==self.fkresults)\
                                  .filter(self.FD_results.sta==self.sta).all()
        if len(id_res)==1:
            return 1

        self.nchan=self.session.query(self.FK_results.nchan).filter(self.FK_results.timeini>=(t0-self.prewindowlen))\
                                      .filter(self.FK_results.timeini<=(te+self.prewindowlen))\
                                      .filter(self.FK_results.pfkid==self.pfkid)\
                                      .filter(self.FK_results.sta==self.sta).distinct().all()
        if len(self.nchan)==0:
            print('no fk results')
            return 1


        if len(self.nchan)>1:
            print('the number of channels changed sometime this day')
            print('this should be addressed in a later implementation')
            print('we will be using the min value')
        try:
            self.nchan=min(self.nchan)
        except Exception as ex1:
            print('problem with data')
            #embed()
            exit()

        self.nchan=self.nchan[0]

        #self.fkdayres=self.session.query(self.FK_results).filter(self.FK_results.timeini>=(t0-self.prewindowlen))\
                                     # .filter(self.FK_results.timeini<=(te+self.prewindowlen))\
                                      #.filter(self.FK_results.pfkid==self.pfkid)\
                                      #.filter(self.FK_results.sourcenum<self.numsources)\
                                      #.filter(self.FK_results.sta==self.sta).order_by(self.FK_results.timeini).all()

        self.fkdayres=self.session.query(self.FK_results).filter(self.FK_results.timeini>=(t0-self.prewindowlen))\
                                      .filter(self.FK_results.timeini<=(te+self.prewindowlen))\
                                      .filter(self.FK_results.pfkid==self.pfkid)\
                                      .filter(self.FK_results.sta==self.sta).order_by(self.FK_results.timeini).all()

        band_w=float(self.fkpar.freqmax)-float(self.fkpar.freqmin)
        win=float(self.fkpar.beamwinlen)
        # we need to retrieve the number of channels....

        num_ch=self.nchan
        self.dofnum = 0.5*(2*win*band_w);
        self.dofden = 0.5*(2*win*band_w*(num_ch-1))


        self.time_lim=self.session.query(func.min(self.FK_results.timeini),func.max(self.FK_results.timeini)).filter(self.FK_results.timeini>(t0-self.prewindowlen))\
                                              .filter(self.FK_results.timeini<(te+self.prewindowlen))\
                                              .filter(self.FK_results.pfkid==self.pfkid)\
                                              .filter(self.FK_results.sta==self.sta)



        if len(self.fkdayres) ==0 :
            return 0
        return 2

    def event_detection(self):
        '''

        '''
        #parameters for fstatistics
        self.x = np.arange(0.01,50.01,0.05)
        self.c_values = np.arange(1.0,100.0,0.5)
        timeW=self.time_lim[0][0]+self.prewindowlen
        time_te=self.time_lim[0][1]
        j=0
        # remember the number of windows is the number of steps wlen
        num_prewindowlen=int(self.prewindowlen/(int(self.fkpar.beamwinlen)-int(self.fkpar.beamwinstep)))

        F2=[]
        pf=[]
        time_pf=[]
        corr=[]
        bz=[]
        sl=[]
        c_found=[]
        times_fk=[]
        fval_fk=[]
        cum_c=[]
        WL=self.prewindowlen
        NWL=num_prewindowlen
        for fki in self.fkdayres:
            times_fk.append(fki.timeini)
            fval_fk.append(fki.fval)
        times_fk=np.asarray(times_fk)
        fval_fk=np.asarray(fval_fk)
        times_fk_effective=[]
        fval_fk_effective=[]


        for t_i in range(len(times_fk)):
            ti=times_fk[t_i]
            if ((ti<self.time_initial) or (ti > self.time_end)) :
                #print 'pre-window or post-window'
                continue

            short_f=fval_fk[(times_fk>(ti-WL)) & (times_fk<(ti+WL))]

            times_fk_effective.append(ti)
            fval_fk_effective.append(fval_fk[t_i])

            if len(short_f) < (2*(NWL)-1):
                #print 'not enough data'
                short_f=np.empty(int(2*(NWL)-1))
                short_f[:]=np.NaN

            cum_c.append(short_f)
        c_val=[]
        #print 'dist fit'
        cum_c_corr=[]

        f_theor = stats.f.pdf(self.x,self.dofnum,self.dofden)
        f_theor=f_theor[:-1]
        #corr = np.corrcoef(f_theor,f_data)[0,1]

        count=0
        for cc in cum_c:
            if np.all(np.isnan(cc))==True:
                c_val.append(np.nan)
                corr.append(0)
            else:
                valsH,bins=np.histogram(cc,self.x)
                ind_valM=np.argmax(valsH)
                valsH,bins=np.histogram(cc/bins[ind_valM],self.x)
                try:
                    corr.append(np.corrcoef(f_theor,valsH)[0,1])
                    c_val.append(bins[ind_valM])
                except Exception as ex1:
                    corr.append(0)
                    c_val.append(np.nan)
                    #print 'lin 394'
            sys.stdout.write("\rPercent: {0}%".format(int(100*float(count)/len(cum_c))))
            sys.stdout.flush()
            count=count+1
        print(' ')

        fval_fk_effective=np.asarray(fval_fk_effective)
        c_val=np.asarray(c_val)
        fval_fk_effective_corr=fval_fk_effective/c_val

        for fvc in fval_fk_effective_corr:
            if np.isnan(fvc)==True:
                pf.append(0)
            else:
                pf.append(stats.f.cdf(fvc ,self.dofnum,self.dofden))
        time_pf=times_fk_effective
        fval_corr=fval_fk_effective_corr

        detections = self.getdet(time_pf,fval_corr,pf,1-self.pthreshold,corr,self.cthr,int(self.minlen/int(self.fkpar.beamwinlen-self.fkpar.beamwinstep)),fval=fval_fk_effective,c=c_val)
        print('dofnum:',self.dofnum,'  dofnum:',self.dofden)
        print('Number of detections in the day:',len(detections))
        if len(detections)>0:
            for idet in range(len(detections)):
                id_resC=self.session.query(self.FD_results).count()+1
                id_res=self.session.query(self.FD_results).filter(self.FD_results.pfdid==self.PFDetect_Q.pfdid)\
                                                  .filter(self.FD_results.pfkid==self.fkpar.pfkid)\
                                                  .filter(self.FD_results.timeini==detections[idet]['ti'])\
                                                  .filter(self.FD_results.timeend==detections[idet]['te'])\
                                                  .filter(self.FD_results.fktablename==self.fkresults)\
                                                  .filter(self.FD_results.sta==self.sta).all()

                if bool(id_res)==False:
                    '''
                    auxD={'ti':time[start_pt[i]],'te':time[stop_pt[i]],/
                                        'tmaxpfc':tmax,'maxpfc':maxpfc,/
                                        'maxfc':maxfc,'maxfo':maxfo,'c':cval}
                    '''
                    res=self.FD_results(fdid=id_resC,\
                            pfdid=self.PFDetect_Q.pfdid,\
                            pfkid=self.fkpar.pfkid,\
                            timeini=detections[idet]['ti'],\
                            timeend=detections[idet]['te'],\
                            maxfc=detections[idet]['maxfc'],\
                            maxfo=detections[idet]['maxfo'],\
                            c=detections[idet]['c'],\
                            sta=self.sta,\
                            fktablename=self.fkresults)
                    self.session.add(res)
                    self.session.commit()
                    #print('det: ', short_time(UTCDateTime(detections[idet]['ti'])),' len(s):',float(detections[idet]['te'])-float(detections[idet]['ti']),' max fval:',detections[idet]['maxfo'])
                else:
                    print('detection already in table')
        else:
            print('No detections in this day')


        id_resC=self.session.query(self.FD_results).count()+1
        id_res=self.session.query(self.FD_results).filter(self.FD_results.pfdid==self.PFDetect_Q.pfdid)\
                                          .filter(self.FD_results.pfkid==self.fkpar.pfkid)\
                                          .filter(self.FD_results.timeini==self.time_lim[0][0])\
                                          .filter(self.FD_results.timeend==self.time_lim[0][1])\
                                          .filter(self.FD_results.maxfc==-1)\
                                          .filter(self.FD_results.fktablename==self.fkresults)\
                                          .filter(self.FD_results.sta==self.sta).all()
        if bool(id_res)==False:
            res=self.FD_results(fdid=id_resC,\
                    pfdid=self.PFDetect_Q.pfdid,\
                    pfkid=self.fkpar.pfkid,\
                    timeini=self.time_initial,\
                    timeend=self.time_end,\
                    fktablename=self.fkresults,\
                    maxfc=-1,\
                    sta=self.sta)
            self.session.add(res)
            self.session.commit()
    def event_detection2(self):

        beam_results = []
        times = []
        for fki in self.fkdayres:
            beam_results = beam_results + [[fki.bz, fki.trvel, fki.fval, fki.nchan]]
            times = times + [[fki.timeini]]

        beam_results=np.array(beam_results)
        times = np.array(times)[:, 0]
        tt = []
        for t in range(len(times)):
            ttt = int(times[t])
            tttt = np.datetime64(ttt,'s')
            tt.append(tttt)
        times = tt
        M = self.nchan
        #detections = beamforming_new.detect_signals(times, beam_results, float(self.detwinlen), det_thresh=float(self.pthresholdeshold), min_seq=int(self.dsegmin), method=self.detmethod, TB_prod=self.tb_prod,channel_cnt=M)

        #detections = beamforming_new.detect_signals(times, beam_results, float(self.detwinlen), det_thresh=float(1-self.pthreshold), min_seq=int(self.dsegmin), back_az_lim=self.backazlim, method=self.detmethod, TB_prod=self.tb_prod,channel_cnt=M, use_det_mask=False)
	detections=beamforming_new.detect_signals(times,beam_results,float(self.adapwlen), float(self.detwinlen)*band_w,M,det_thresh=float(1-self.pthreshold),min_seq=int(self.dsegmin),back_az_lim=self.backazlim,fixed_thresh=det_thresh)
        if len(detections)>0:
            for idet in range(len(detections)):
                ti = UTCDateTime(detections[idet][0].astype(datetime)).timestamp
                te = ti + float(self.detwinlen)
                id_resC=self.session.query(self.FD_results).count()+1
                id_res=self.session.query(self.FD_results).filter(self.FD_results.pfdid==self.PFDetect_Q.pfdid)\
                                                  .filter(self.FD_results.pfkid==self.fkpar.pfkid)\
                                                  .filter(self.FD_results.timeini==ti)\
                                                  .filter(self.FD_results.timeend==te)\
                                                  .filter(self.FD_results.fktablename==self.fkresults)\
                                                  .filter(self.FD_results.sta==self.sta).all()

                if bool(id_res)==False:

                    res=self.FD_results(fdid=id_resC,\
                            pfdid=self.PFDetect_Q.pfdid,\
                            pfkid=self.fkpar.pfkid,\
                            timeini=ti,\
                            timeend=te,\
                            tonset=detections[idet][1],\
                            tend=detections[idet][2],\
                            trvel=detections[idet][4],\
                            c=detections[idet][5],\
                            sta=self.sta,\
                            fktablename=self.fkresults)

                    self.session.add(res)
                    self.session.commit()
                    #print('det: ', short_time(UTCDateTime(detections[idet]['ti'])),' len(s):',float(detections[idet]['te'])-float(detections[idet]['ti']),' max fval:',detections[idet]['maxfo'])
                else:
                    print('detection already in table')
        else:
            print('No detections in this day')


    def distfit(self,x,F,dofnum,dofden,c):
        '''
        Returns the correlation coefficient between the theoretical F-distribution and the
        distribution of the input F-statistics, scaled by a c-value

        Inputs:
        - x is the set of F-statistics for binning
        - F are the measured F-statistics
        - dofnum is the number of degrees of freedom on the numerator
        - dofden is the number of degrees of freedom on the denominator
        - c is the c value

        Outputs:
        - corr is the correlation coefficient
        '''
        try:

            Fmod = F/c
        except:
            return 0
        f_theor = stats.f.pdf(x,dofnum,dofden)
        f_data = np.histogram(Fmod,x)[0]
        x = x[0:len(x)-1]
        f_theor = f_theor[0:len(f_theor)-1]
        corr = np.corrcoef(f_theor,f_data)[0,1]
        return corr

    def getdet(self,time,fval_corr,pf,p,corr,corr_thres,min_no_samples,fval=None,c=None):
        '''
        Obtains the start and end points for significant detections given a p-value.
        Inputs:
        - pf is the array of p-values as a function of time
        - p is the threshold p-value
        Outputs:
        - start_pt is the list of indices corresponding to the starts of detections
        - stop_pt is the list of indices corresponding to the ends of detections
        '''
        in_det = 0        # Flag that indicates if detector is within a detection
        # Finding start/stop points for detections:
        start_pt = []; stop_pt = []; i = -1
        for pf_i in pf:
            i = i + 1
            corr_i = corr[i]
            if (in_det == 0):
                if (pf_i >= p and corr_i >= corr_thres):
                    in_det = 1
                    start_pt.append(i)
                else:
                    continue
            else:
                if (pf_i >= p and corr_i >= corr_thres):
                    continue
                else:
                    in_det = 0
                    stop_pt.append(i)

        # Removing detection at the end:
        if (len(start_pt) > len(stop_pt)):
            stop_pt.append(len(corr)-1)

        detections=[]
        for i in range(0,len(start_pt)):
            # in (stop_pt[i] - start_pt[i]) >= min_no_samples the >=  is required! if min_no_samples=1
            if ((stop_pt[i] - start_pt[i]) >= min_no_samples):
                pfmax_i=np.argmax(pf[start_pt[i]:stop_pt[i]])+start_pt[i]
                tmax=time[pfmax_i]
                maxpfc=pf[pfmax_i]
                maxfc=fval_corr[pfmax_i]
                if (fval is not None):
                    maxfo=fval[pfmax_i]
                else:
                    maxfo=-1
                if (c is not None):
                    cval=c[pfmax_i]
                else:
                    cval=-1
                auxD={'ti':time[start_pt[i]],'te':time[stop_pt[i]],\
                                    'tmaxfc':tmax,'maxpfc':maxpfc,\
                                    'maxfc':maxfc,'maxfo':maxfo,'c':cval}
                detections.append(auxD)
        return detections


    def data_processing(self):
        '''
        Constructor
        '''
        for day in np.arange(self.dayofyearini,self.dayofyearend):
            print('year:',self.year,'day:', day)
            self.dayofyear=day
            retrieve_ret=self.data_retrieving(day)
            if retrieve_ret==2:
                if self.domain == 'time':
                    print('running event detection in time domain')
                    self.event_detection()
                elif self.domain == 'freq':
                    print('running event detection in freq domain')
                    self.event_detection2()
                else:
                    print('invalid domain for data processing')
            elif retrieve_ret==1 :
                print('day already analyzed')
            elif retrieve_ret==0:
                print('No FK analysis for this day ')
            else:
                print('there was an error or not data in fk table')

        #return self.data_processed



if __name__ == '__main__':
    pdetect=FDetectInfraPy_LANL('../conf_files/IMS_CONF_FILES/InfraConfigI57US')
    pdetect.database_connecting()
    pdetect.data_processing()
    #embed()
