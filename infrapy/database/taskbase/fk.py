 #!/usr/bin/env python -W ignore::DeprecationWarning
'''
Updated Fall 2019
@fkdd @omarcillo
'''

from .base import Base

#from ...utils.cart2pol import cart2pol

from ...utils.get_arraywaveforms import get_arraywaveforms
from .. import schema
from ...detection import beamforming
from ...detection import beamforming_new

import sys, pdb

import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import func
from sqlalchemy import MetaData

import pisces as ps
from pisces.util import load_config

from obspy.core import UTCDateTime
from obspy.core import trace
from obspy.core import Stream
from obspy.signal.array_analysis import *
import obspy.signal.invsim as inv
import numpy as np
from scipy import signal
from IPython import embed
import cmath
import math
import itertools
import matplotlib
import re as re
import time
from numpy import linalg as LA
from multiprocessing import Pool
import multiprocessing
import itertools
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

import timeit

status=0
inte=0


def stream2array(str):
    """ turn obspy stream into 2D array

    the data from stream is extracted and copied to a 2D numpy array
     - ** parameters**, **types**, **return**, and **return types**::
         :param arg1: stream to be converted
        :type arg1: obspy stream
        :return: 2D array with time series for all channels
        :rtype: numpy array
    """

    tI=[]
    tE=[]
    for sti in str:
        tI.append(float(sti.stats.starttime))
        tE.append(float(sti.stats.endtime))
    try:
        tIMAX=np.max(np.asarray(tI))
        tEMIN=np.min(np.asarray(tE))
    except Exception as ex1:
        print('error at stream2array, FK:',ex1)
        embed()
        exit()
    str.trim(starttime=UTCDateTime(tIMAX),endtime=UTCDateTime(tEMIN))
    num_samples=[]
    for sti in str:
        num_samples.append(sti.stats.npts)
    if len(np.unique(num_samples))>1:
        print('there is an issue with the number of samples,')
        print('maybe different sampling rates')
        return []
    else:
        for st_i in range(len(str)):
            if str[st_i].stats.calib==0:
                print("there is an error with calib number")
                print('setting calib to 1')
                str[st_i].stats.calib=1
            if st_i==0:
                mat=str[st_i].data*str[st_i].stats.calib
            else:
                mat=np.vstack((mat,str[st_i].data*str[st_i].stats.calib))
    return mat

def stream2mult_array(str,wlen,overlap):
    """ turn obspy stream into array of short 2D arrays

    the data from a stream is extracted, sliced, and copied to 2D numpy arrays
     - ** parameters**, **types**, **return**, and **return types**::
         :param str: stream to be converted
        :param wlen: length of short arrays
        :param overlap: length of overlap between short arrays

        :type str: obspy stream
        :type wlen: length of windows in seconds
        :type overlap: length of overlap in seconds

        :return: 2D array with time series for all channels
        :rtype: numpy array

        :return: array of times of the beginning of the windows
        :rtype: list of UTCDateTime

        :return: array of window lengths
        :rtype: list of seconds
    """
    arr=stream2array(str)
    sps=str[0].stats.sampling_rate
    npts=str[0].stats.npts
    timeI=str[0].stats.starttime
    wlen_sam=sps*wlen
    overlap_sam=sps*overlap
    curr=0
    mat=[]
    time=[]
    lenA=[]
    try:
        cosVal=inv.cosTaper(int(wlen_sam),p=0.1)
    except:
        cosVal=inv.cosine_taper(int(wlen_sam),p=0.1)
    while curr+wlen_sam<=npts:
        mat.append(arr[:,int(curr):int(curr+wlen_sam)]*cosVal)
        time.append(UTCDateTime(float(timeI)+float(curr)/sps))
        lenA.append(len(arr[0,int(curr):int(curr+wlen_sam)]))
        curr=curr+wlen_sam-overlap_sam
    mat=np.asarray(mat)
    return mat,time,lenA

def extARRAYPROC_ZIP(argu):
    """ allow processing the data using multicore capabilities

    This is a temporary method, as python cannot process using multicores within classes.
    So we have to call a method from within the class to run multicores. you provide a
     - ** parameters**, **types**, **return**, and **return types**::
         :param argu: list of lists with differeent elelemtn
        :type argu: list of lists
        :return: list of results of each of the element of the analysis
        :rtype: numpy array
    """

    global inte
    if inte==1:
        return []
    try:
        global status
        try:
            if argu[14]=='bartlett':
                method='bartlett'
            elif argu[14]=='music':
                method='music'
            elif argu[14]=='music32':
                method='music32'
            elif argu[14]=='music16':
                method='music16'
            elif argu[14]=='music32D':
                method='music32D'
            elif argu[14]=='music16D':
                method='music16D'
            else:
                print('the method requested for processing is not available')
                exit()
                #method='bartlett'
            #ret=beamforming.fkBARLET(argu[0],argu[1],argu[2],argu[3],argu[4],argu[5],argu[6],argu[8],argu[10],argu[11],argu[12],argu[13])
            ret=beamforming.fkPROC(argu[14],argu[0],argu[1],argu[2],argu[3],argu[4],argu[5],argu[6],argu[8],argu[10],argu[11],argu[12],argu[13])
        except Exception as ex1:
            print('problem processing the data():',ex1)
            ret=[0,0,-1,0,0,0,0,0,0]
        status=status+1
        if status%10==0:
            sys.stdout.write("\rPercent: {0}%".format(int(argu[9]*round(status*100.0/argu[7]))))
            sys.stdout.flush()
    except KeyboardInterrupt:
        print(inte)
        inte=1
        ret=[0,0,-1,0,0,0,0,0,0]
    return ret

class Fk(Base):
    ''' FK processing main class
    This class provides all the features to process array data (infrasound and seismic)
    using a database. In case a fully featured database is not available, this tool can interact with
    a serverless database.
    '''
    algorithm='None'
    def __init__(self, conf_file=[],ARRAY_NAME=None):
        """ initialize the processing parameters

        The configuration file is read and parameters are extracted for processing
         - ** parameters**, **types**, **return**, and **return types**::
             :param conf_file: configuration file
            :param ARRAY_NAME (optional): overrides the array selected in the configuration file
            :type conf_file: text file
            :type ARRAY_NAME: string name of the reference station to be processed

            :return: no return
            :rtype: no return
        """
        super(Fk,self).__init__(conf_file,'FKParams')

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


        self.fkversion=0
        self.name=self.task_PARAM['name']
        self.freqmin=float(self.task_PARAM['freqmin'])
        self.freqmax=float(self.task_PARAM['freqmax'])
        self.beamwinlen=float(self.task_PARAM['beamwinlen'])
        self.beamwinstep=float(self.task_PARAM['beamwinstep'])

        try:
            self.backazmin=float(self.task_PARAM['backazmin'])
        except:
            print('bzmin set by default to -180')
            self.backazmin=-180
        try:
            self.backazmax=float(self.task_PARAM['backazmax'])
        except:
            print('bzmax set by default to 180')
            self.backazmax=180
        try:
            self.backazstep=float(self.task_PARAM['backazstep'])
        except:
            print('bz step set by default to 1.5')
            self.backazstep=1.5
        try:
            self.trvelmin=float(self.task_PARAM['trvelmin'])
        except:
            print('trace vel min set by default to 300')
            self.trvelmin=300
        try:
            self.trvelmax=float(self.task_PARAM['trvelmax'])
        except:
            print('trace vel max set by default to 600')
            self.trvelmax=600
        try:
            self.trvelstep=float(self.task_PARAM['trvelstep'])
        except:
            print('trace vel step set by default to 2.5')
            self.trvelstep=2.5
        try:
            self.minslowness=float(self.task_PARAM['minslowness'])
        except:
            print('minimum slowness set by default to -3.6')
            self.minslowness=-3.6
        try:
            self.maxslowness=float(self.task_PARAM['maxslowness'])
        except:
            print('maximum slowness set by default to 3.6')
            self.maxslowness=3.6
        try:
            self.stepslowness=float(self.task_PARAM['stepslowness'])
        except:
            print('slowness step set by default to 0.1')
            self.maxslowness=0.1


        self.slow=np.arange(self.minslowness,self.maxslowness+self.stepslowness,self.stepslowness)

        try:
            self.resultstable=self.task_PARAM['fkresults']
        except KeyError:
            self.resultstable='Fk_results'


        try:
            self.decimate=int(self.task_PARAM['decimate'])
        except KeyError:
            self.decimate=1

        try:
            self.numsources=int(self.task_PARAM['numsources'])
        except KeyError:
            self.numsources=1

        try:
            self.func_fk=self.task_PARAM['func_fk']
        except KeyError:
            self.func_fk=None


        self.jdayini = int(self.year*1000 + self.dayofyearini)
        self.jdayend = int(self.year*1000 + self.dayofyearend)


        try:
            self.algorithm=self.task_PARAM['algorithm']
        except:
            print('method set by default to Bartlett')
            self.algorithm='bartlett'

        if self.algorithm=='bartlett':
            print('Bartlett method was selected')
        elif self.algorithm=='capon':
            print('Capon method was selected')
        elif self.algorithm=='music':
            print('Music method was selected')
        elif self.algorithm=='music16':
            print('Music method was selected, div 16')
        elif self.algorithm=='music32':
            print('Music method was selected, div 32')
        elif self.algorithm=='music32D':
            print('Music method was selected, div 32 with dynamic number of sources')
        elif self.algorithm=='music16D':
            print('Music method was selected, div 16 with dynamic number of sources')
        else:
            print('This method is not implemented here')
            print('select either bartlett or capon method and run the script again')
            sys.exit()

    def database_connecting(self):
        """ connect to database and write parameters for analysis

        The system connects to the database and establishes the health of the connection, and writes the parameters for processing
         - ** parameters**, **types**, **return**, and **return types**::
            :return: state of the connection True (successful), FALSE (problem), for connection or
            :rtype: boolean

        Update: Feb 2018
        DB connection established with a user-defined configuration file to specify which standard database tables a user is targeting.
        Pisces provides a function that reads the format and returns the necessary tables. The configuration file also specifies infrapy-specific database tables for fk processing
        """

        session,tables = load_config(self.db_PARAM)
        self.session=session
        self.Site=tables['site']
        self.Wfdisc=tables['wfdisc']


        class Fk_results(schema.fk_results):
                __tablename__ = self.resultstable

        class Fk_params(schema.fk_params):
                __tablename__ = 'FK_PARAMS'
        self.FK_par=Fk_params
        self.FK_results=Fk_results
        self.FK_par.__table__.create(self.session.bind,checkfirst=True)
        self.FK_results.__table__.create(self.session.bind,checkfirst=True)


        try:
            self.FK_parQuery=self.session.query(self.FK_par). \
                filter(self.FK_par.freqmax==self.freqmax).filter(self.FK_par.freqmin==self.freqmin). \
                filter(self.FK_par.beamwinlen==self.beamwinlen).filter(self.FK_par.beamwinstep==self.beamwinstep).\
                filter(self.FK_par.backazmin==self.backazmin).\
                filter(self.FK_par.backazmax==self.backazmax).\
                filter(self.FK_par.backazstep==self.backazstep).\
                filter(self.FK_par.trvelmin==self.trvelmin).\
                filter(self.FK_par.trvelmax==self.trvelmax).\
                filter(self.FK_par.trvelstep==self.trvelstep).\
                filter(self.FK_par.name==self.name).\
                filter(self.FK_par.minslowness==self.minslowness).\
                filter(self.FK_par.maxslowness==self.maxslowness).\
                filter(self.FK_par.stepslowness==self.stepslowness).\
                filter(self.FK_par.name==self.name).\
                filter(self.FK_par.numsources==self.numsources).\
                filter(self.FK_par.domain==self.domain).\
                filter(self.FK_par.algorithm==self.algorithm).\
                all()

            if len(self.FK_parQuery)>1:
                print('issue with the database too many parameters entries, there should be just one')
                embed()
            if len(self.FK_parQuery)==1:
                self.FK_parQuery=self.FK_parQuery[0]
        except Exception as ex1:
            print("issue with the table,",ex1)
            embed()
            self.FK_parQuery=[]
        if bool(self.FK_parQuery)==False:
            new_row=self.session.query(self.FK_par).count()
            print('New process parameters, write process to fk_params table pfkid=',new_row)
            res=self.FK_par(freqmax=self.freqmax,\
                                freqmin=self.freqmin,\
                                beamwinlen=self.beamwinlen,\
                                beamwinstep=self.beamwinstep,\
                                backazmin=self.backazmin,\
                                backazmax=self.backazmax,\
                                backazstep=self.backazstep,\
                                trvelmin=self.trvelmin,\
                                trvelmax=self.trvelmax,\
                                trvelstep=self.trvelstep,\
                                minslowness=self.minslowness,\
                                maxslowness=self.maxslowness,\
                                stepslowness=self.stepslowness,\
                                name=self.name, \
                                numsources=self.numsources, \
                                domain = self.domain, \
                                algorithm=self.algorithm, \
                                pfkid=new_row)
            self.session.add(res)
            try:
             self.session.commit()
            except Exception as ex1:
             print(ex1,", there is a problem with these parameters")
             embed()
            self.FK_parQuery=self.session.query(self.FK_par). \
                filter(self.FK_par.freqmax==self.freqmax).filter(self.FK_par.freqmin==self.freqmin). \
                filter(self.FK_par.beamwinlen==self.beamwinlen).filter(self.FK_par.beamwinstep==self.beamwinstep).\
                filter(self.FK_par.backazmin==self.backazmin).\
                filter(self.FK_par.backazmax==self.backazmax).\
                filter(self.FK_par.backazstep==self.backazstep).\
                filter(self.FK_par.trvelmin==self.trvelmin).\
                filter(self.FK_par.trvelmax==self.trvelmax).\
                filter(self.FK_par.trvelstep==self.trvelstep).\
                filter(self.FK_par.minslowness==self.minslowness).filter(self.FK_par.maxslowness==self.maxslowness).\
                filter(self.FK_par.stepslowness==self.stepslowness).\
                filter(self.FK_par.name==self.name).\
                filter(self.FK_par.numsources==self.numsources).\
                filter(self.FK_par.domain==self.domain).\
                filter(self.FK_par.algorithm==self.algorithm)\
                .one()
            #embed()
        else:
            print('process already in  fk_params table, pfkid:', self.FK_parQuery.pfkid)

        self.db_connected=True
        return self.db_connected

    def data_retrieving(self,dayofyear):
        """ retrieve one day of data from the database and write back results

         - ** parameters**, **types**, **return**, and **return types**::
             :param dayofyear: day of year that will be retrieved
            :type dayofyear: integer
            :type ARRAY_NAME: string name of the reference station to be processed

            :return: integer with state of health of the stream, 0 problem with retrieving the data; 1 this data was already analyzed; 2 successful retrieving the
                      data stream
            :rtype: integer

        """
        t0 = UTCDateTime( year=self.year, julday=dayofyear).timestamp
        te = t0 + 86400
        self.time_initial=t0
        self.time_end=te
        #embed()
        id_res=self.session.query(self.FK_results).filter(self.FK_results.timeini>t0)\
                                              .filter(self.FK_results.timeini<te)\
                                              .filter(self.FK_results.chan==self.chan)\
                                              .filter(self.FK_results.pfkid==self.FK_parQuery.pfkid)\
                                              .filter(self.FK_results.sta==self.sta).all()

        if len(id_res)>0:
            id_res_max=self.session.query(func.max(self.FK_results.timeend)).filter(self.FK_results.timeini>t0)\
                                                  .filter(self.FK_results.timeini<te)\
                                                  .filter(self.FK_results.chan==self.chan)\
                                                  .filter(self.FK_results.pfkid==self.FK_parQuery.pfkid)\
                                                  .filter(self.FK_results.sta==self.sta).all()
            if id_res_max[0][0] >= (te-self.beamwinlen):
                return 1
        else:
            print('no previous results in the database')
        try:
            self.aa=get_arraywaveforms(self.session,self.Site,self.Wfdisc,self.sta,t0=t0,te=te,channel=self.chan)
            if len(self.aa[0])<60:
                return 0
        except Exception as ex1:
            print('problem with get_arraywaveforms',ex1)
            return 0
        if len(self.aa.traces)==0:
            print('maybe a channel different than:', self.chan)
            print('no data')
            return 0
        num_s=[]
        for aai in self.aa:
            num_s.append(aai.stats.npts)
        num_s=np.asarray(num_s)
        if not(len(np.unique(num_s)) == 1):
            print('there is a problem with one channel')
            return 0
        return 2

    def array_processing(self):
        t0=time.time()
        try:
            self.aa.detrend()
        except Exception as ex1:
            print("Problem with detrending(353):",ex1)
            for aai in self.aa:
                if (type(aai.data) is np.ma.core.MaskedArray) == True:
                    aai.data=aai.data-aai.data.mean()
                else:
                    aai.data=aai.data-aai.data.mean()
        if not( self.decimate==1):
            self.aa.decimate(self.decimate)

        sps=self.aa[0].stats.sampling_rate
        [b,a]=signal.butter(3,[self.freqmin/(sps/2),self.freqmax/(sps/2)],'bandpass')
        for aai in self.aa:
            try:
                aai.data=signal.filtfilt(b,a,aai.data)
            except Exception as ex1:
                print('problem filtering 572')
                embed()
                exit()
        num_traces=len(self.aa.traces)
        timeI=[]
        timeE=[]
        for ti in self.aa:
            timeI.append(ti.stats.starttime)
            timeE.append(ti.stats.endtime)
        timeIT=np.max(timeI)
        timeET=np.min(timeE)

        window_l=float(self.FK_parQuery.beamwinlen)
        overwindow_l=float(self.FK_parQuery.beamwinstep)
        # info required for the FK

        X,Y=beamforming.getXY_array(self.aa)
        self.lat2MeanKM=np.asarray(Y,dtype='float')
        self.lon2MeanKM=np.asarray(X,dtype='float')
        sps=self.aa[0].stats.sampling_rate
        win_l=self.FK_parQuery.beamwinlen
        x=self.lon2MeanKM
        y=self.lat2MeanKM
        nchan = float(len(self.lat2MeanKM))
        slow=self.slow
        self.fN,self.freqsN,self.mult_vectors=beamforming._mult_vectors(win_l,sps,slow,x,y,nchan,fmin=self.freqmin,fmax=self.freqmax)

        id_resC=self.session.query(self.FK_results).count()+1
        count_min=0
        id_resCHECK=self.session.query(self.FK_results).filter(self.FK_results.timeini>=self.time_initial)\
                                                  .filter(self.FK_results.timeini<=self.time_end)\
                                                  .filter(self.FK_results.pfkid==self.FK_parQuery.pfkid)\
                                                  .filter(self.FK_results.sta==self.sta)\
                                                  .order_by(self.FK_results.timeini)
        if id_resCHECK.count()==0:
            print('No previous results with this ID, all day will be processed')
            #timeRES_ini=0
            timeRES_end=0
            timeW_0=UTCDateTime(timeIT)
        else:
            timeRES_ini=id_resCHECK[-1].timeini
            timeRES_end=id_resCHECK[-1].timeend
            timeW_0=UTCDateTime(timeRES_ini+window_l-overwindow_l)
            print('Previous results with this ID, starting at:',UTCDateTime(timeW_0))
            if (timeW_0 > self.aa[0].stats.endtime) or ((timeW_0 +window_l)> self.aa[0].stats.endtime):
                print('there is not enough data to process this stream')
                return self.data_processed
            else:
                self.aa=self.aa.slice(starttime=timeW_0,endtime=self.aa[0].stats.endtime)

        timeW=timeW_0
        count_min=0
        count_ST=0
        arraySTREAM_time=[]
        counS=0
        count_min
        timeW=timeW_0
        arraySTREAM,arraySTREAM_time,lenARRAY=stream2mult_array(self.aa,window_l,overwindow_l)
        out=[]
        #argumentSTREAM=[]
        argumentSTREAMS=[]
        global status
        status=0
        self.cpucnt = int(self.cpucnt)
        for tre_i in range(len(arraySTREAM)):
            argumentSTREAMS.append([arraySTREAM[tre_i],sps,self.slow,self.mult_vectors,self.fN,X,Y,len(arraySTREAM),arraySTREAM_time[tre_i],self.cpucnt,self.func_fk,self.freqsN,self.sta,self.numsources,self.algorithm])
        if self.cpucnt>1:
            print('Available cores :', multiprocessing.cpu_count())
            print('Cores to be used:',self.cpucnt)
            if self.cpucnt >= multiprocessing.cpu_count():
                print('There is only:',multiprocessing.cpu_count(), 'available')
                print('infrapy will be using only', multiprocessing.cpu_count()-1)
                self.cpucnt=multiprocessing.cpu_count()-1
            pool = Pool(processes=self.cpucnt)
            try:
                out=pool.map(extARRAYPROC_ZIP,argumentSTREAMS)
            except KeyboardInterrupt:
                pool.terminate()
                print('Do not finish processing and discard values')
                sys.exit(0)
            pool.close()
        else:
            print('multicore capabilities disabled (standard one core calculation)')
            for aar in range(len(argumentSTREAMS)):
                val=extARRAYPROC_ZIP(argumentSTREAMS[aar])
                out.append(val)

        id_resC=self.session.query(self.FK_results).count()+1

        for r_T in range(len(out)):
            for r_i in range(len(out[r_T])):
                try:
                    res=self.FK_results(pfkid=self.FK_parQuery.pfkid,\
                                            timeini=float(arraySTREAM_time[r_T]),timeend=float(arraySTREAM_time[r_T]+window_l),\
                                            fkid=id_resC,sta=self.sta,\
                                            bz=out[r_T][r_i][0],slofk=out[r_T][r_i][1],\
                                            fval=out[r_T][r_i][2],xcorrvalmin=out[r_T][r_i][3],xcorrvalmax=out[r_T][r_i][4],\
                                            xcorrvalmean=out[r_T][r_i][5],\
                                            rmsval=out[r_T][r_i][6],\
                                            sx=out[r_T][r_i][7],\
                                            sy=out[r_T][r_i][8],\
                                            chan=self.chan,\
                                            sourcenum=r_i,\
                                            nchan=len(self.aa))
                    id_resC=id_resC+1
                    self.session.add(res)
                except Exception as ex1:
                    print('There was a problem writing the results to DB')
                    print('The error is :', ex1)
                    embed()
                    exit(0)
            if int(len(out)/10)>0:
                if r_T%int(len(out)/10)==0 or  r_T==(len(out)-1):
                    print(UTCDateTime(arraySTREAM_time[r_T]))
                    self.session.commit()
        try:
            #print 'not saving 578'
            self.session.commit()
        except:
            print('Nothing to commit')

        print('processing time:',time.time()-t0,' s')
        return self.data_processed

    def array_processing_2(self):

        t0=time.time()

        lo=[]
        la=[]
        for ii in range(len(self.aa)):
            try:
                self.aa[ii].stats.coordinates
                coord=0
            except Exception as ex1:
                #print('No coordinates')
                coord=1
            if coord==0:
                lo.append(self.aa[ii].stats.coordinates.longitude)
                la.append(self.aa[ii].stats.coordinates.latitude)
            else:
                lo.append(self.aa[ii].stats.sac['stlo'])
                la.append(self.aa[ii].stats.sac['stla'])
        loc = np.array([la,lo])
        latlon = np.transpose(loc)
        x, t, t1, geom = beamforming_new.stream_to_array_data(self.aa,latlon)
        back_az_vals = np.arange(self.backazmin, self.backazmax, self.backazstep)
        trc_vel_vals = np.arange(self.trvelmin, self.trvelmax, self.trvelstep)
        slowness = beamforming_new.build_slowness(back_az_vals, trc_vel_vals)
        delays = beamforming_new.compute_delays(geom, slowness)

        pl = Pool(int(self.cpucnt))
        num_traces=len(self.aa.traces)
        timeI=[]
        timeE=[]
        for ti in self.aa:
            timeI.append(ti.stats.starttime)
            timeE.append(ti.stats.endtime)
        timeIT=np.max(timeI)
        timeET=np.min(timeE)


        window_l=float(self.FK_parQuery.beamwinlen)
        overwindow_l=float(self.FK_parQuery.beamwinstep)

        id_resC=self.session.query(self.FK_results).count()+1
        count_min=0
        id_resCHECK=self.session.query(self.FK_results).filter(self.FK_results.timeini>=self.time_initial)\
                                                  .filter(self.FK_results.timeini<=self.time_end)\
                                                  .filter(self.FK_results.pfkid==self.FK_parQuery.pfkid)\
                                                  .filter(self.FK_results.sta==self.sta)\
                                                  .order_by(self.FK_results.timeini)
        if id_resCHECK.count()==0:
            print('No previous results with this ID, all day will be processed')
            #timeRES_ini=0
            timeRES_end=0
            timeW_0=UTCDateTime(timeIT)
        else:
            timeRES_ini=id_resCHECK[-1].timeini
            timeRES_end=id_resCHECK[-1].timeend
            timeW_0=UTCDateTime(timeRES_ini+window_l-overwindow_l)
            print('Previous results with this ID, starting at:',UTCDateTime(timeW_0))
            if (timeW_0 > self.aa[0].stats.endtime) or ((timeW_0 +window_l)> self.aa[0].stats.endtime):
                print('there is not enough data to process this stream')
                return self.data_processed
            else:
                self.aa=self.aa.slice(starttime=timeW_0,endtime=self.aa[0].stats.endtime)
        times, beam_results = [], []
        for win_start in np.arange(t[0], t[-1], self.beamwinstep):
            if win_start + self.beamwinlen > t[-1]:
                break

            times = times + [[t1 + np.timedelta64(int(win_start), 's')]]

            X, S, f = beamforming_new.fft_array_data(x, t, window=[win_start, win_start + self.beamwinlen])
            beam_power = beamforming_new.run(X, S, f, geom, delays, [self.freqmin, self.freqmax], method=self.algorithm, pool=pl, normalize_beam=True)
            peaks = beamforming_new.find_peaks(beam_power, back_az_vals, trc_vel_vals, signal_cnt=1)
            beam_results = beam_results + [[peaks[0][0], peaks[0][1], peaks[0][2], peaks[0][2] / (1.0 - peaks[0][2]) * (x.shape[0] - 1)]]
            #embed()

        times = np.array(times)[:,0]
        tt = []

        for t in range(len(times)):
            ttt = UTCDateTime(times[t].astype(datetime))
            tt.append(ttt.timestamp)
        beam_results = np.array(beam_results)
        pl.close()
        id_resC=self.session.query(self.FK_results).count()+1

        for r_T in range(len(beam_results)):
            try:
                res=self.FK_results(pfkid=self.FK_parQuery.pfkid,\
                                        timeini=tt[r_T],timeend=tt[r_T] + (int(self.beamwinstep)),\
                                        fkid=id_resC,sta=self.sta,\
                                        bz=beam_results[r_T][0],trvel=beam_results[r_T][1],\
                                        coher=beam_results[r_T][2],\
                                        fval=beam_results[r_T][3],\
                                        chan=self.chan,\
                                        sourcenum=self.numsources,\
                                        nchan=len(self.aa))
                id_resC=id_resC+1
                self.session.add(res)
            except Exception as ex1:
                print('There was a problem writing the results to DB')
                print('The error is :', ex1)
                embed()
                exit(0)

        try:
            #print 'not saving 578'
            self.session.commit()
        except:
            print('Nothing to commit')

        print('processing time:',time.time()-t0,' s')
        return self.data_processed


    def data_processing(self):
        """ process the data in self.aa

         - ** parameters**, **types**, **return**, and **return types**::
            :return: state of processing True (succesful), FALSE (problem)
            :rtype: boolean

        """
        for day in np.arange(self.dayofyearini,self.dayofyearend):
            print('year:',self.year,'day:', day)
            retrieve_ret=self.data_retrieving(day)
            if retrieve_ret==2:
                if self.domain == 'time':
                    self.array_processing()
                elif self.domain == 'freq':
                    self.array_processing_2()
                else:
                    print('invalid domain for data processing')

            elif retrieve_ret==1 :
                print('day already analyzed')
            else:
                print('there was an error or not data')
        return self.data_processed

if __name__ == '__main__':
        """ example for using the tool

        """

        print('FK step for InfraPy')
        try:
            AR=sys.argv[2]
            print('run FD, for array:',AR)
        except Exception as ex1:
            print('Array to be processed from configuration file')
            AR=None

        print('run FK, with configuration param:',DB)
        pdetect=InfrapyFk(DB,ARRAY_NAME=AR)
        pdetect.database_connecting()
        pdetect.data_processing()
