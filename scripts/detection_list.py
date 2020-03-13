## Fransiska K Dannemann Dugick
## fransiska at lanl dot gov
## 08.06.2019

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
from infrapy.utils.get_header_table import get_header_table

import time

import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print list of detections to be used for assocation/location processing")
    parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
    parser.add_argument('-a', dest='array',required=True,help="array name, e.g.: -a I37NO")
    parser.add_argument('-f','--pfkid', dest='pfkid',required=True,help='fk parameter id, e.g.: -f 0')
    parser.add_argument('-j','--pfdid', dest='pfdid',required=True,help='fd parameter id, e.g.: -j 0')
    parser.add_argument('-t', dest='fdresults',required=False,help="specific table with results, e.g.: -t fd_I37")
    parser.add_argument('-s', dest='tS',required=False,help="starttime plot, e.g.: -s /'2014-03-02T00:00:00/'")
    parser.add_argument('-e', dest='tE',required=False,help="endtime plot, e.g.: -s /'2014-03-03T00:00:00/'")
    parser.add_argument('-F', dest='fval',required=False,help="limit Fvalue, e.g.: -F 0")
    parser.add_argument('-bzmin', dest='bzmin',required=False,help="limit min bz, e.g.: -bzmin 0")
    parser.add_argument('-bzmax', dest='bzmax',required=False,help="limit max bz, e.g.: -bzmax 360")
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

    if args.fval:
        fval=args.fval
        print('setting fval:', fval)
    else:
        fval=0

    if args.bzmin:
        bzmin=args.bzmin
        print('setting bzmin:', bzmin)
    else:
        bzmin=None


    if args.bzmax:
        bzmax=args.bzmax
        print('setting bzmax:', bzmax)
    else:
        bzmax=None


    if args.outtext:
        outtext=args.outtext
        fid=open(outtext,'w')
    else:
        print("no outfile specified")
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



    if  len(fd_res)==0:
        print('No results with pfkid:', pfkid, ' and pfdid:',pfdid)
        exit()
    #STRH=    'fdid\t'+' pfdid\t'+'pfkid\t'+'timeini\t'+'timeini(machine)\t'+'length\t'+'timeMAX\t'+'Cval\t'+'fvalMAX\t'+'xcorrMAX\t'+'bzMAX\t'+'trveMAX\n'
    STRH=    '#Lat\t'+' Lon\t'+'Time\t'+'Back Az\t'+'F-stat\t'+'Array Dim\n'
    if not (fid==None):
        fid.write(STRH)

    siteQ=session.query(Site).filter(Site.refsta==array).all() ## need to query here for only array named

    try:
        class Fk_results(schema.fk_results):
            __tablename__ = fd_res[0].fktablename
    except Exception as ex1:
        pass
    all_res=session.query(Fk_results).all()
    timeI = np.array([row.timeini for row in all_res])

    k=0
    for fd_i in fd_res:
        timeINI=time.time()
        t_INI=UTCDateTime(fd_i.timeini)
        t_length=float(UTCDateTime(fd_i.timeend))-float(UTCDateTime(fd_i.timeini))

        fk_res=session.query(Fk_results).filter(Fk_results.pfkid==fd_i.pfkid).filter(Fk_results.sta==array).filter(Fk_results.timeini.between(fd_i.timeini,fd_i.timeend)).all()
        #fk_res=all_res[np.argwhere((timeI>=fd_i.timeini)&(timeI<=fd_i.timeend))]
        #all_res.filter(Fk_results.pfkid==fd_i.pfkid).filter(Fk_results.sta==array).filter(Fk_results.timeini.between(fd_i.timeini,fd_i.timeend)).all()
        ind=np.argwhere((timeI>=fd_i.timeini)&(timeI<=fd_i.timeend))
        fval=[]
        bz=[]
        chan=[]

        for fk_i in fk_res:
            fval.append(fk_i.fval)
            bz.append(fk_i.bz)
            chan.append(fk_i.nchan)


        fval=np.asarray(fval)
        try:
            fvalm_i=np.argmax(fval)
            STRD=str(siteQ[0].lat)+'\t'+str(siteQ[0].lon)+'\t'+str(UTCDateTime(t_INI))+'\t'
            STRD=STRD+'{:.1f}'.format(bz[fvalm_i])+'\t'+'{:.3f}'.format(fval[fvalm_i])+'\t'+'{:.3f}'.format(chan[fvalm_i])+'\n'
            if not (fid==None):
                fid.write(STRD)
        except Exception as ex1:
            #embed()
            pass


            # Lat 	Lon			Time					Back Az	F-stat	Array Dim
        #STRD=str(siteQ[0].lat)+'\t'+str(siteQ[0].lon)+'\t'+str(short_time(t_INI))+'\t'
        #STRD=STRD+'{:.1f}'.format(bz[fvalm_i])+'\t'+'{:.3f}'.format(fval[fvalm_i])+'\t'+'{:.3f}'.format(chan[fvalm_i])+'\n'
        #if not (fid==None):
            fid.write(STRD)


    fid.close()
