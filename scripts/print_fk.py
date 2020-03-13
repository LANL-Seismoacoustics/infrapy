## Fransiska K Dannemann Dugick
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


import infrapy.database.schema as schema

class Fk_params(schema.fk_params):
    __tablename__ = 'FK_PARAMS'

import getopt
import argparse

from infrapy.utils.get_header_table import get_header_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read fk results for specific array and FK parameter ID')
    parser.add_argument('-d', dest='sq',required=True,help='name of the database connection, e.g.: -d sqlite:///mydb.sqlite')
    parser.add_argument('-a', dest='array',required=True,help='array name, e.g.: -a I37NO')
    parser.add_argument('-t', dest='fkresults',required=False,help='specific table with results, e.g.: -t fk_I37')
    parser.add_argument('-f', dest='pfk_id',required=True,help='FK parameter ID to be printed, e.g.: -f 3')
    parser.add_argument('-s', dest='tS',required=False,help="starttime, e.g.: -s /'2015-03-02T00:00:00/'")
    parser.add_argument('-e', dest='tE',required=False,help="endtime, e.g.: -s /'2015-03-03T00:00:00/'")
    parser.add_argument('-o',dest='outtext',required=False,help='outfile, e.g.: -o res_FILE.txt')


    args = parser.parse_args()
    if args.sq:
        sq=args.sq
    if args.array:
        array=args.array

    if args.fkresults:
        fkresultsT=args.fkresults
    else:
        print('default table FK_RESULTS')
        fkresultsT='FK_RESULTS'

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
        print('no outfile specified')
        outtext=None

    class Fk_results(schema.fk_results):
        __tablename__ = fkresultsT

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

        else:
            print('No standard database, try anyway')
            session=ps.db_connect(self.database)
            session_type='unknown'
    except Exception as ex1:
        print('Connection failed:', ex1)
        embed()
        sys.exit()

    if not (outtext==None):
        fid.write('\n'+str(get_header_table(Fk_results)))

    #fk_par=session.query(Fk_params).filter(Fk_params.pfkid==pfk_id)
    res=session.query(Fk_results).filter(Fk_results.pfkid==pfk_id).filter(Fk_results.sta==array)
    if not (t_S==None):
        res=res.filter(Fk_results.timeini>=float(t_S)).filter(Fk_results.timeini<=float(t_E)).order_by(Fk_results.timeini).all()
    else:
        res=res.order_by(Fk_results.timeini).all()

    for res_i in res:
        fid.write('\n'+str(res_i))
        print(UTCDateTime(res_i.timeini))

    fid.close()
