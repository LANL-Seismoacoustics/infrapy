

#!/usr/bin/env python
import sys, pdb
import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
#from pisces.io.trace import read_waveform
from obspy.core import UTCDateTime
from obspy.core import trace
from obspy.core import Stream
from obspy.core.util import AttribDict
from obspy.signal.array_analysis import *
from datetime import datetime

import numpy as np
import scipy as sc
from IPython import embed
import pisces as ps

from pisces import request
from sqlalchemy import func


import matplotlib


#matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import matplotlib.mlab as mpy
import pylab as py

import matplotlib.dates as mdates


import infrapy.database.schema as schema
from infrapy.utils.short_time import short_time

import time

import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print fd results for specific array and FK/FD combination")
    parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
    parser.add_argument('-a', dest='array',required=True,help="array name, e.g.: -a I37NO")
    parser.add_argument('-f','--pfkid', dest='pfkid',required=True,help='fk parameter id, e.g.: -f 0')
    parser.add_argument('-j','--pfdid', dest='pfdid',required=True,help='fd parameter id, e.g.: -j 0')
    parser.add_argument('-t', dest='fdresults',required=False,help="specific table with detection results, e.g.: -t fd_I37")
    parser.add_argument('-s', dest='tS',required=False,help="starttime limit, e.g.: -s /'2015-03-02T00:00:00/'")
    parser.add_argument('-e', dest='tE',required=False,help="endtime limit, e.g.: -s /'2015-03-03T00:00:00/'")
    parser.add_argument('-o',dest='outtext',required=False,help='outtext file name, e.g.: -o res_FILE')

    args = parser.parse_args()
    if args.sq:
        sq=args.sq
    if args.array:
        array=args.array

    if args.fdresults:
        fdresultsT=args.fdresults
    else:
        print("default table FD_RESULTS")
        fdresultsT='FD_RESULTS'

    if args.pfkid:
        pfkid=args.pfkid

    if args.pfdid:
        pfdid=args.pfdid

    if args.tS:
        t_S=UTCDateTime(args.tS)
        print('setting ini time:', t_S)
    else:
        t_S=None

    if args.tE:
        t_E=UTCDateTime(args.tE)
        print('setting end time:', t_E)
    else:
        t_E=None


    if args.outtext:
        outtext=args.outtext
        fid=open(outtext,'w')
    else:
        print("default table FD_RESULTS")
        outtext=None
        fid=None


    class Fd_params(schema.fd_params):
        __tablename__ = 'FD_PARAMS'

    class Fd_results(schema.fd_results):
        __tablename__ = fdresultsT

    try:
        if sq.count('oracle')>0:
            session=ps.db_connect(sq)
            session_type='oracle'
            from global_ import Site, Origin, Wfdisc_raw

        elif sq.count('sqlite')>0:
            print('SQLITE database')
            session=ps.db_connect(sq)
            session_type='sqlite'
            from pisces.tables.kbcore import Site, Origin, Wfdisc
            '''
            class Site(kba.Site):
                __tablename__ = 'site'

            class Wfdisc(kba.Wfdisc):
                __tablename__ = 'wfdisc'
            '''
        else:
            print('No standard database, try anyway')
            session=ps.db_connect(self.database)
            session_type='unknown'
    except Exception as ex1:
        print('Connection failed:', ex1)
        embed()
        sys.exit()


    #fd_par=session.query(Fd_params).filter(Fd_results.pfdid==pfdid).all()
    fd_res=session.query(Fd_results).filter(Fd_results.sta==array).filter(Fd_results.pfkid==pfkid).filter(Fd_results.pfdid==pfdid).all()
    if not (t_S==None):
        fd_res=fd_res.filter(Fd_results.timeini>=float(t_S)).filter(Fd_results.timeini<=float(t_E)).order_by(Fd_results.timeini).all()
    else:
        fd_res=fd_res
    #fd_res=session.query(Fd_results).filter(Fd_results.sta==array).filter(Fd_results.pfkid==pfkid).filter(Fd_results.pfdid==pfdid).filter(Fd_results.maxfo>0).all()



    if  len(fd_res)==0:
        print('No results with pfkid:', pfkid, ' and pfdid:',pfdid)
        exit()
    #STRH=    'fdid\t'+' pfdid\t'+'pfkid\t'+'timeini\t'+'timeini(machine)\t'+'length\t'+'timeMAX\t'+'Cval\t'+'fvalMAX\t'+'xcorrMAX\t'+'bzMAX\t'+'trveMAX\n'
    STRH=    'fdid\t'+' pfdid\t'+'pfkid\t'+'timeini\t'+'timeini(machine)\t'+'length\t'+'Cval\t'+'fval\t'+'xcorr\t'+'bz\t'+'trvel\n'
    if not (fid==None):
        fid.write(STRH)

    print(STRH, end=' ')
    TT0=time.time()

    try:
        class Fk_results(schema.fk_results):
            __tablename__ = fd_res[0].fktablename
    except Exception as ex1:
        pass
    all_res=session.query(Fk_results).all()
    print(1000*(time.time()-TT0))
    timeI = np.array([row.timeini for row in all_res])
    print(1000*(time.time()-TT0))
    #pdb.set_trace()
    k=0
    for fd_i in fd_res:
        timeINI=time.time()
        t_INI=UTCDateTime(fd_i.timeini)
        t_length=float(UTCDateTime(fd_i.timeend))-float(UTCDateTime(fd_i.timeini))
        #t_delay=float(UTCDateTime(fd_i.maxfc_time))-float(UTCDateTime(fd_i.timeini))
        #pdb.set_trace()
        '''
        try:
            class Fk_results(ab.fk_results):
                __tablename__ = fd_i.fktablename
        except Exception, ex1:
            pass
        '''

        TT1=time.time()
        if fd_i.maxfc == -1:
            continue
        else:
            fk_res=session.query(Fk_results).filter(Fk_results.pfkid==fd_i.pfkid).filter(Fk_results.sta==array).filter(Fk_results.timeini.between(fd_i.timeini,fd_i.timeend)).all()
        #fk_res=all_res[np.argwhere((timeI>=fd_i.timeini)&(timeI<=fd_i.timeend))]
        #all_res.filter(Fk_results.pfkid==fd_i.pfkid).filter(Fk_results.sta==array).filter(Fk_results.timeini.between(fd_i.timeini,fd_i.timeend)).all()
            ind=np.argwhere((timeI>=fd_i.timeini)&(timeI<=fd_i.timeend))
        #pdb.set_trace()
            TT2=time.time()
            fval=[]
            xcorr=[]
            bz=[]
            tr=[]

            #fk_res=[]
            #for iir in ind:
            #    fk_res.append(all_res[iir])

            for fk_i in fk_res:
                fval.append(fk_i.fval)

                bz.append(fk_i.bz)
                if fk_i.slofk != -1.0:
                    tr.append(1./np.asanyarray(fk_i.slofk))
                else:
                    tr.append(fk_i.trvel)
                if fk_i.xcorrvalmax != -1.0:
                    xcorr.append(fk_i.xcorrvalmax)
                else:
                    xcorr.append(fk_i.coher)

            fval=np.asarray(fval)
            fvalm_i=np.argmax(fval)
            #embed()
            STRD=str(fd_i.fdid)+'\t'+str(fd_i.pfdid)+'\t'+str(fd_i.pfkid)+'\t'+str(short_time(t_INI))+'\t'+str(float(t_INI))+'\t'+str(t_length)+'\t'+str(fd_i.c)+'\t'+str(fd_i.maxfo)+'\t'
            STRD=STRD+'{:.2f}'.format(xcorr[fvalm_i])+'\t'+'{:.1f}'.format(bz[fvalm_i])+'\t'+'{:.3f}'.format(tr[fvalm_i])+'\n'
            if not (fid==None):
                fid.write(STRD)
            #print STRD
            #print 1000*(TT1-timeINI),1000*(TT2-timeINI),1000*(time.time()-timeINI)
            if (k%100)==0:
                print(UTCDateTime(fd_i.timeini))
            k=k+1

    sys.exit(0)
