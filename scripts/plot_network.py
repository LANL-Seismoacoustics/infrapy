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

pl.ion()

from infrapy.utils.cart2pol import cart2pol
#sys.path.insert(0, '../tables')
import infrapy.database.schema as schema
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Plot locations of arrays within network")
    parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
    parser.add_argument('-a', dest='array',required=False,help="array name, e.g.: -a I37NO")
    args = parser.parse_args()
    if args.sq:
        sq=args.sq
    if args.array:
        array=args.array
    else:
        print("plotting everything")
        arrays=[]

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

    if len(array)==0:
        siteQ=session.query(Site).all()
    else:
        siteQ=session.query(Site).all() ## need to query here for only array named
    for site_i in siteQ:
        pl.plot(site_i.lon,site_i.lat,'o')
        pl.text(site_i.lon,site_i.lat,site_i.refsta)
    pl.axis('equal')
    embed()
    pl.show()
    #embed()
