#!/usr/bin/env python
import sys, pdb
import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import func
import sqlalchemy.exc as exc
import sqlalchemy.orm.exc as oexc
import configparser

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
from pisces.util import load_config
from pisces.io.trace import read_waveform


import pylab as py
import matplotlib
import matplotlib.pyplot as pl
import matplotlib.mlab as mpy
import matplotlib.dates as mdates
import calendar as ca

pl.ioff()


# db2sac.py I37NO 2015-11-30T00:00 2015-12-01

if __name__ == '__main__':

    STA=sys.argv[1]
    timeini=sys.argv[2]
    timeend=sys.argv[3]
    try:
        chan=sys.argv[4]
    except Exception as ex1:
        chan='BDF'
    ttI=UTCDateTime(timeini)
    ttE=UTCDateTime(timeend)
    jdayI = int(ttI.year*1000 + ttI.day)
    jdayE = int(ttE.year*1000 + ttE.day)

    config = configparser.ConfigParser()
    config.read('gnem.cfg')
    session, tables = load_config(config['database'])
    session=session
    Site=tables['site']
    Wfdisc_raw=tables['wfdisc']

    print('from ',UTCDateTime(timeini), 'to ', UTCDateTime(timeend))


    si_res = session.query(Site).filter(Site.refsta == STA).filter(Site.ondate <= jdayI).filter(Site.offdate >= jdayE).all()
    for si in si_res:
        wfres=session.query(Wfdisc_raw).filter(Wfdisc_raw.sta==si.sta).filter(Wfdisc_raw.time>float(ttI)).filter(Wfdisc_raw.time<float(ttE)).filter(Wfdisc_raw.chan==chan).all()
        print('working on:', si.sta)
        for wf in wfres:
               trace = wf.to_trace()
               trace.stats.sac=AttribDict({'stla':si.lat,'stlo': si.lon,'stel': si.elev})
               trace.write(STA+str(wf.wfid) + '.sac', format='SAC')
               print('writing ',STA+str(wf.wfid) + '.sac')
