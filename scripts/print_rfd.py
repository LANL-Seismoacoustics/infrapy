

#!/usr/bin/env python
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pl

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

import matplotlib.mlab as mpy
import pylab as py

import matplotlib.dates as mdates
import infrapy.database.schema as schema
from infrapy.utils.get_header_table import get_header_table

class Fk_params(schema.fk_params):
    __tablename__ = 'FK_PARAMS'

from infrapy.utils.get_arraywaveforms import get_arraywaveforms
import argparse


if __name__ == '__main__':
    #
    parser = argparse.ArgumentParser(description="Print contents of detection table list for a specific FD/FK processing combination")
    parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
    parser.add_argument('-a', dest='array',required=True,help="array name, e.g.: -a I37NO")
    parser.add_argument('-f','--pfkid', dest='pfk_id',required=True,help="FK parameter ID, e.g.: -f 3")
    parser.add_argument('-j','--pfdid', dest='pfd_id',required=True,help="FD parameter ID, e.g.: -j 3")
    parser.add_argument('-t', dest='fdresults',required=False,help="specific table with detection results, e.g.: -t fd_I37")
    parser.add_argument('-s', dest='tS',required=False,help="starttime limit, e.g.: -s /'2015-03-02T00:00:00/'")
    parser.add_argument('-e', dest='tE',required=False,help="endtime limit, e.g.: -s /'2015-03-03T00:00:00/'")
    parser.add_argument('-o',dest='outtext',required=False,help='file to save output, e.g.: -o res_FILE')

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

    if args.pfd_id:
        pfd_id=args.pfd_id

    if args.pfk_id:
        pfk_id=args.pfk_id

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

    if not (outtext==None):
        fid.write('\n'+str(get_header_table(Fd_results)))


    #fk_par=session.query(Fk_params).filter(Fk_params.pfkid==pfk_id)
    res=session.query(Fd_results).filter(Fd_results.pfkid==pfk_id).filter(Fd_results.pfdid==pfd_id).filter(Fd_results.sta==array)
    if not (t_S==None):
        res=res.filter(Fd_results.timeini>=float(t_S)).filter(Fd_results.timeini<=float(t_E)).order_by(Fd_results.timeini).all()
    else:
        res=res.order_by(Fd_results.timeini).all()

    for res_i in res:
        fid.write('\n'+str(res_i))
        print(UTCDateTime(res_i.timeini))

    fid.close()
