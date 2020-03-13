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



from infrapy.utils.get_arraywaveforms import get_arraywaveforms
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read fk results for specific array and FK parameter ID")
    parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
    parser.add_argument('-a', dest='array',required=True,help="array name, e.g.: -a I37NO")
    parser.add_argument('-f','--pfkid', dest='pfk_id',required=True,help="FK parameter ID to be plot, e.g.: -f 3")
    parser.add_argument('-t', dest='fkresults',required=False,help="specific table with results, e.g.: -t fk_I37")
    parser.add_argument('-domain', dest='domain',required=True,help="domain for processing, e.g.: -domain freq")
    parser.add_argument('-s', dest='tS',required=False,help="starttime plot, e.g.: -s /'2014-03-02T00:00:00/'")
    parser.add_argument('-e', dest='tE',required=False,help="endtime plot, e.g.: -s /'2014-03-03T00:00:00/'")
    parser.add_argument('-F', dest='fval',required=False,help="limit Fvalue, e.g.: -F 0")
    parser.add_argument('-bzmin', dest='bzmin',required=False,help="limit min bz, e.g.: -bzmin 0")
    parser.add_argument('-bzmax', dest='bzmax',required=False,help="limit max bz, e.g.: -bzmax 360")
    parser.add_argument('-fname', dest='fname',required=False,help="output file name e.g.: -fname test_run1")

    args = parser.parse_args()
    if args.sq:
        sq=args.sq
    if args.array:
        array=args.array

    if args.domain:
        domain=args.domain

    if args.fkresults:
        fkresultsT=args.fkresults
    else:
        print("default table FK_RESULTS")
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

    if args.fval:
        fval=args.fval
        print('setting fval:', fval)
    else:
        fval=0

    if args.bzmin:
        bzmin=args.bzmin
        print('setting bzmin:', bzmin)
    else:
        if domain == 'time':
            bzmin=0
        else:
            bzmin = -180

    if args.bzmax:
        bzmax=args.bzmax
        print('setting bzmax:', bzmax)
    else:
        if domain == 'time':
            bzmax=360
        else:
            bzmax=180


    if args.fname:
        fname=args.fname
        print('setting outputfile:', fname)
        name=1
    else:
        fname=None
        name=0

    class Fk_results(schema.fk_results):
        __tablename__ = fkresultsT
        
    class Fk_params(schema.fk_params):
        __tablename__ = 'FK_PARAMS'


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


    fk_par=session.query(Fk_params).filter(Fk_params.pfkid==pfk_id).all()
    for fd_i in fk_par:
        res=session.query(Fk_results).filter(Fk_results.pfkid==fd_i.pfkid).filter(Fk_results.sta==array).filter(Fk_results.fval>=fval).filter(Fk_results.bz>=bzmin,Fk_results.bz<=bzmax).all()
        #res=session.query(Fk_results)
        t_ini=session.query(func.min(Fk_results.timeini)).filter(Fk_results.pfkid==fd_i.pfkid).filter(Fk_results.sta==array)
        t_end=session.query(func.max(Fk_results.timeini)).filter(Fk_results.pfkid==fd_i.pfkid).filter(Fk_results.sta==array)


        print('start plotting')
        fig1=pl.figure(figsize=(24,6))
        ax1=pl.subplot(1,1,1)

        fig2=pl.figure(figsize=(24,6))
        ax2=pl.subplot(1,1,1,sharex=ax1)

        fig3=pl.figure(figsize=(24,6))
        ax3=pl.subplot(1,1,1,sharex=ax1)



        tM_t=[]
        bz_t=[]
        fval_t=[]

        slofk_t=[]

        for res_i in res:
            tt1=UTCDateTime(res_i.timeend)
            tt2=UTCDateTime(res_i.timeini)
            mean_time=(float(tt1)+float(tt2))/2
            tM_t.append(UTCDateTime(mean_time).datetime)
            #tM_t.append(UTCDateTime(res_i.timeini))
            bz_t.append(res_i.bz)
            fval_t.append(res_i.fval)
            #embed()
            if res_i.slofk != -1.0:
                slofk_t.append(1./np.asanyarray(res_i.slofk))
            else:
                slofk_t.append(res_i.trvel)

        formatDATE= mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
        ax1.xaxis.set_major_formatter(formatDATE)

        bz_t=np.asarray(bz_t)
        if domain == 'time':
            im = ax1.scatter(tM_t,bz_t,c=slofk_t, cmap=pl.cm.jet, edgecolors='none',vmin=0.2,vmax=0.4)
        else:
            im = ax1.scatter(tM_t,bz_t,c=slofk_t, cmap=pl.cm.jet, edgecolors='none',vmin=200,vmax=400)

        colorbar = pl.colorbar(im, ax=ax1)
        colorbar.set_label('Trace Velocity',fontsize=12)
        ax2.plot_date(tM_t,fval_t,'.k')
        ax3.plot_date(tM_t,slofk_t,'.k')
        #cbar=pl.colorbar(slofk_t, cmap='jet', ax=ax1)
        #cbar.ax.set_aspect(40)

        #cbar.ax.set_position([box3.x0+1.02*box3.width,box3.y0,box3.width,box3.height])
        if domain == 'time':
            ax1.set_ylim(0,360)
        else:
            ax1.set_ylim(-180,180)



        ax1.set_ylabel('Back-Azimuth (degrees)',fontsize=14,fontweight='bold')
        ax2.set_ylabel('F Value',fontsize=14,fontweight='bold')
        #ax2.set_yscale('log')
        ax3.set_ylabel('Trace Velocity (km/s)',fontsize=14,fontweight='bold')


        if not t_S==None:
            ax1.set_xlim(mdates.date2num(t_S),mdates.date2num(t_E))
        else:
            #embed()
            ax1.set_xlim(mdates.date2num(tM_t[0]),mdates.date2num(tM_t[-1]))

        ax1.grid('on')
        ax2.grid('on')
        ax3.grid('on')
        #str_t=str(UTCDateTime.now())
        if name==1:
            fig1.savefig(fname+'_backaz_trvel.png')
            fig2.savefig(fname+'_fvals.png')
            fig3.savefig(fname+'_trvels.png')


        pl.show()
