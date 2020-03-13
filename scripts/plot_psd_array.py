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

#from cart2pol import cart2pol
#sys.path.insert(0, '../tables')
#import infrapy.database.schema as ab
import matplotlib.patches as patches


from infrapy.utils.get_arraywaveforms import get_arraywaveforms
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="calculate PSD for an specific array")
    parser.add_argument('-d', dest='sq',required=True,help='name of the database connection, e.g.: -d sqlite:///mydb.sqlite')
    parser.add_argument('-a', dest='array',required=True,help='array name, e.g.: -a I37NO')
    parser.add_argument('-S', dest='tS',required=True,help="starttime plot, e.g.: -S /'2015-03-02T00:00:00/'")
    parser.add_argument('-E', dest='tE',required=True,help="endtime plot, e.g.: -E /'2015-03-03T00:00:00/'")
    parser.add_argument('-l', dest='wlen',required=False,help="subwindow length in seconds, e.g.: -l 100")
    

    args = parser.parse_args()
    if args.sq:
        sq=args.sq
    if args.array:
        array=args.array
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

    if args.wlen:
        wlen=float(args.wlen)
        print('setting end time:', t_E)
    else:
        wlen=100

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


    aaTT=get_arraywaveforms(session,Site,Wfdisc,array,t_S,t_E)

    print('printing all waveform')
    #aaT=get_arraywaveforms(session,Site,Wfdisc,array,fd_res[ii].timeini-0.5*dt,fd_res[ii].timeend+0.5*dt)
    aaT=aaTT.slice(starttime=t_S,endtime=t_E)
    aaT.filter('highpass',freq=0.1)
    fig=pl.figure()
    ax1 = pl.subplot(211)
    ax2 = pl.subplot(212,sharex=ax1)

    timeD=np.linspace(float(aaT[0].stats.starttime), float(aaT[0].stats.endtime), aaT[0].stats.npts)
    ax1.plot(timeD-timeD[0],aaT[0].data*aaT[0].stats.calib)
    samps=aaT[0].stats.npts
    nfft=int(aaT[0].stats.sampling_rate*wlen)
    x=aaT[0].data*aaT[0].stats.calib
    Pxx, freqs,t=mpy.specgram(x, NFFT=nfft, Fs=aaT[0].stats.sampling_rate, noverlap=nfft/2)
    ref=20E-6

    aa=ax2.pcolormesh(t,freqs,10*np.log10(Pxx/(ref*ref)))
        #aa=ax2.pcolormesh(t,freqs,10*np.log10(Pxx/(ref*ref)))
    ref=20E-6
    box3=ax2.get_position()
    #ax2.set_yscale('log')
    cbar1=pl.colorbar(aa,ax=ax2)
    ax2.set_position(box3)
    ax2.set_ylabel('Frequency (Hz)')
    cbar1.ax.set_aspect(40)
    cbar1.set_label('Power (dB/Hz rel. 20 uPa)',fontsize=12)
    cbar1.ax.set_position([box3.x0+1.02*box3.width,box3.y0,box3.width,box3.height])
    ax2.set_yticks([0.1,1,2,5,10])
    ax2.set_yticklabels([0.1,1,2,5,10])
    ax2.set_ylim([0.1,10])
    ax2.set_xlabel('Time (s) after '+aaT[0].stats.starttime.strftime('%y-%m-%d %H:%M:%S'))
    ax1.set_ylabel('Pressure (Pa) ')
    ax1.set_title(array)
    t_i=timeD[0]-timeD[0]
    t_e=timeD[-1]-timeD[0]
    ax1.set_xlim(t_i,t_e)
    pl.figure()
    pl.semilogy(freqs,np.mean(Pxx,axis=1))
    ax=pl.gca()
    ax.set_xlim(0.1,10)
    pl.show()
