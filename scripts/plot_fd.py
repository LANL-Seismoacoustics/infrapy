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
from numpy import append
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
import infrapy.database.schema as schema
import matplotlib.patches as patches


class Fk_params(schema.fk_params):
    __tablename__ = 'fk_params'

class Fd_params(schema.fd_params):
    __tablename__ = 'fd_params'

from infrapy.utils.get_arraywaveforms import get_arraywaveforms
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="read and plot all the f detetions for an specific array")
    parser.add_argument('-d', dest='sq',required=True,help='name of the database connection, e.g.: -d sqlite:///mydb.sqlite')
    parser.add_argument('-a', dest='array',required=True,help='array name, e.g.: -a I37NO')
    parser.add_argument('-f','--pfkid', dest='pfkid',required=True,help='FK parameter ID to be plot, e.g.: -f 3')
    parser.add_argument('-j','--pfdid', dest='pfdid',required=True,help='FD parameter ID to be plot, e.g.: -j 0')
    parser.add_argument('-t', dest='fkresults',required=False,help="specific table with results, e.g.: -t fk_I37")
    parser.add_argument('-T', dest='fdresults',required=False,help="specific table with results, e.g.: -T fd_I37")
    parser.add_argument('-s', dest='tS',required=False,help="starttime plot, e.g.: -s /'2015-03-02T00:00:00/'")
    parser.add_argument('-e', dest='tE',required=False,help="endtime plot, e.g.: -s /'2015-03-03T00:00:00/'")

    args = parser.parse_args()
    if args.sq:
        sq=args.sq
    if args.array:
        array=args.array
    if args.pfkid:
        pfkid=int(args.pfkid)
    if args.pfdid:
        pfdid=int(args.pfdid)

    if args.fkresults:
        fkresultsT=args.fkresults
    else:
        print("default table FK_RESULTS")
        fkresultsT='FK_RESULTS'
    if args.fdresults:
        fdresultsT=args.fdresults
    else:
        print("default table FD_RESULTS")
        fdresultsT='FD_RESULTS'

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


    class Fk_results(schema.fk_results):
        __tablename__ = fkresultsT


    class Fd_results(schema.fd_results):
        __tablename__ = fdresultsT


    print('Using pfdid:',pfdid,' and pfkid:',pfkid)

    if t_E is None:
        fd_res=session.query(Fd_results).filter(Fd_results.pfdid==pfdid).filter(Fd_results.pfkid==pfkid).filter(Fd_results.sta==array).filter(Fd_results.maxfc>0).all()
    else:
        fd_res=session.query(Fd_results).filter(Fd_results.pfdid==pfdid).filter(Fd_results.pfkid==pfkid).filter(Fd_results.timeini>=float(t_S)).filter(Fd_results.timeini<=float(t_E)).filter(Fd_results.sta==array).filter(Fd_results.maxfc>0).all()
    fk_par=session.query(Fk_params).filter(Fk_params.pfkid==pfkid).all()


    aaTT=get_arraywaveforms(session,Site,Wfdisc,array,t_S,t_E)
    aa0=aaTT.copy()
    aaTT.filter('bandpass',freqmin=fk_par[0].freqmin,freqmax=fk_par[0].freqmax)
    aaTT.taper(0.1)
    h1=aaTT.plot(show=False,handle=True)
    ax1=h1.get_axes()

    for ii in range(len(fd_res)):
        aaT=get_arraywaveforms(session,Site,Wfdisc,array,fd_res[ii].timeini,fd_res[ii].timeend)
        fk_res=session.query(Fk_results).filter(Fd_results.pfkid==pfkid).filter(Fk_results.sta==array).filter(Fk_results.timeini>=fd_res[ii].timeini).filter(Fk_results.timeini<=fd_res[ii].timeend).all()
        x0=mdates.date2num(UTCDateTime(fd_res[ii].timeini).datetime)
        xF=mdates.date2num(UTCDateTime(fd_res[ii].timeend).datetime)
        for axi in ax1:
             yl=axi.get_ylim()
             axi.add_patch(patches.Rectangle((x0,yl[0]),(xF-x0),yl[1]-yl[0],edgecolor='yellow',facecolor='yellow'))


    print('printing all waveform')
    for ii in range(len(fd_res)):
        t0=fd_res[ii].timeini
        te=fd_res[ii].timeend
        dt=te-t0
        aaTTA=aa0.copy()
        aaTTA.filter('bandpass',freqmin=fk_par[0].freqmin,freqmax=fk_par[0].freqmax)
        aaTTA.taper(0.1)
        aaT=aaTTA.trim(starttime=UTCDateTime(fd_res[ii].timeini)-0.5*dt,endtime=UTCDateTime(fd_res[ii].timeend)+0.5*dt)
        #aaT=get_arraywaveforms(session,Site,Wfdisc,array,fd_res[ii].timeini-0.5*dt,fd_res[ii].timeend+0.5*dt)
        fk_res=session.query(Fk_results).filter(Fk_results.sta==array).filter(Fk_results.pfkid==fd_res[ii].pfkid).filter(Fk_results.timeini>=fd_res[ii].timeini).filter(Fk_results.timeini<=fd_res[ii].timeend).all()
        fig=pl.figure()
        ax1 = pl.subplot(211)
        ax2 = pl.subplot(212)

        fig2=pl.figure()
        bx1 = pl.subplot(211)
        bx2 = pl.subplot(212)


        timeD=np.linspace(float(aaT[0].stats.starttime), float(aaT[0].stats.endtime), aaT[0].stats.npts)

        ax1.plot(timeD-timeD[0],aaT[0].data*aaT[0].stats.calib)
        bx1.plot(timeD-timeD[0],aaT[0].data*aaT[0].stats.calib)

        nfft=int(aaT[0].stats.npts/32)
        x=aaT[0].data*aaT[0].stats.calib
        Pxx, freqs,t=mpy.specgram(x, NFFT=nfft, Fs=aaT[0].stats.sampling_rate, noverlap=nfft/2)
        ref=20E-6
        #aa=ax2.pcolormesh(t,freqs,10*np.log10(Pxx/(ref*ref)),vmin=-20,vmax=60)
        #embed()
        PxxS=Pxx[(freqs>=fk_par[0].freqmin)&(freqs<=fk_par[0].freqmax),:]
        freqsS=freqs[(freqs>=fk_par[0].freqmin)&(freqs<=fk_par[0].freqmax)]

        aa=ax2.pcolormesh(t,freqsS,10*np.log10(PxxS/(ref*ref)))
        #aa=ax2.pcolormesh(t,freqs,10*np.log10(Pxx/(ref*ref)))
        ref=20E-6
        box3=ax2.get_position()
        ax2.set_yscale('log')
        cbar1=pl.colorbar(aa,ax=ax2)
        ax2.set_position(box3)
        ax2.set_ylabel('Frequency (Hz)')
        cbar1.ax.set_aspect(40)
        cbar1.set_label('Power (dB/Hz rel. 20 uPa)',fontsize=12)
        cbar1.ax.set_position([box3.x0+1.02*box3.width,box3.y0,box3.width,box3.height])
        ax2.set_yticks([0.1,1,2,5,10])
        ax2.set_yticklabels([0.1,1,2,5,10])
        ax2.set_ylim([fk_par[0].freqmin,fk_par[0].freqmax])
        ax2.set_xlabel('Time (s) after '+aaT[0].stats.starttime.strftime('%y-%m-%d %H:%M:%S'))
        ax1.set_ylabel('Pressure (Pa) ')
        ax1.set_title(array)

        tt=[]
        ff=[]
        for fk_i in fk_res:
            tt.append((fk_i.timeini+fk_i.timeend)/2)
            ff.append(fk_i.fval)
        tt=np.asarray(tt)
        ff=np.asarray(ff)
        bx2.plot(tt-float(aaT[0].stats.starttime),ff ,'.')
    pl.show()
