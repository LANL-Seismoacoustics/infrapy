#!/usr/bin/env python

## test


import sys, pdb
import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
import sqlalchemy.exc as exc
import sqlalchemy.orm.exc as oexc
from sqlalchemy import func

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


import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pl
import matplotlib.mlab as mpy
import pylab as py

import matplotlib.dates as mdates
import infrapy.database.schema as schema

import argparse

import pdb


if __name__ == '__main__':

    '''
    db2db.py -d oracle://:username:port/account -a I18DK -f 1 -t resultsI18DK -o resultsN
    '''

    parser = argparse.ArgumentParser(description="Export FK table from one DB to another")
    parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
    parser.add_argument('-a', dest='array',required=True,help="array name, e.g.: -a I37NO")
    parser.add_argument('-f','--pfkid', dest='pfk_id',required=True,help="FK parameter ID to be extracted, e.g.: -f 3 or -f *")
    parser.add_argument('-t', dest='fkresults',required=False,help="specific table with results, e.g.: -t fk_I37")
    parser.add_argument('-o', dest='outputtable',required=True,help="specific sqlite name with results, e.g.: -o fk_I37")

    parser.add_argument('-s', dest='tS',required=False,help="starttime plot, e.g.: -s /'2015-03-02T00:00:00/'")
    parser.add_argument('-e', dest='tE',required=False,help="endtime plot, e.g.: -s /'2015-03-03T00:00:00/'")


    args = parser.parse_args()
    if args.sq:
        sq=args.sq
    if args.array:
        array=args.array

    if args.fkresults:
        fkresultsT=args.fkresults
    else:
        print("default table FK_RESULTS")
        fkresultsT='FK_RESULTS'

    if args.pfk_id:
        if args.pfk_id=='*':
            print('Extract all pfk_id')
            pfk_id=None
        else:
            pfk_id=args.pfk_id
    else:
        pfk_id=None

    if args.outputtable:
        outputtable=args.outputtable

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

    if args.wf:
        if args.wf==1:
            wf=1
        else:
            wf=0
    else:
        wf=0

    class Fk_params(schema.fk_params):
        __tablename__ = 'FK_PARAMS'

    class Fk_results(schema.fk_results):
        __tablename__ = fkresultsT

    #import pdb; pdb.set_trace()


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


    fk_par=session.query(Fk_params).all()

    if pfk_id==None:
        res=session.query(Fk_results).filter(Fk_results.sta==array).all()
        t_ini=session.query(func.min(Fk_results.timeini)).filter(Fk_results.sta==array).one()
        t_end=session.query(func.max(Fk_results.timeini)).filter(Fk_results.sta==array).one()
    else:
        res=session.query(Fk_results).filter(Fk_results.pfkid==pfk_id).filter(Fk_results.sta==array).all()
        t_ini=session.query(func.min(Fk_results.timeini)).filter(Fk_results.pfkid==pfk_id).filter(Fk_results.sta==array).one()
        t_end=session.query(func.max(Fk_results.timeini)).filter(Fk_results.pfkid==pfk_id).filter(Fk_results.sta==array).one()

    site_res=session.query(Site).filter(Site.refsta==array).all()


    sessionOUT=ps.db_connect('sqlite:///'+outputtable)
    sessionOUT.commit()
    Fk_params.__table__.create(sessionOUT.bind,checkfirst=True)
    Fk_results.__table__.create(sessionOUT.bind,checkfirst=True)
    try:
        Wfdisc.__table__.create(sessionOUT.bind,checkfirst=True)
        Site.__table__.create(sessionOUT.bind,checkfirst=True)
    except Exception as ex1:
        import pisces.tables.kbcore as kb
        kb.Wfdisc.__table__.create(sessionOUT.bind,checkfirst=True)
        kb.Site.__table__.create(sessionOUT.bind,checkfirst=True)
    resS=[]


    print('Loading the entire fk params table')

    for ii in fk_par:
        try:
            sessionOUT.add(Fk_params(*ii))
            sessionOUT.commit()
        except exc.IntegrityError as e:
            sessionOUT.rollback()
            print('this FK param ID already exists')

    print('Loading fk results')
    numT=5*int(len(res)/100.0)
    cnt=0
    NT=len(res)
    print('number of fk results:',NT)
    for ii in res:
        try:
            sessionOUT.add(Fk_results(*ii))
            sessionOUT.commit()
        except exc.IntegrityError as e:
            print('error',e)
            sessionOUT.rollback()
        if (cnt%numT==0):
            print(100*float(cnt)/NT)
            sessionOUT.commit()
        cnt=cnt+1
    sessionOUT.commit()
    print('Loading site and wfdisc info')

    for ssi in site_res:
        print(ssi.sta)
        try:
            sessionOUT.add(kb.Site(*ssi))
            sessionOUT.commit()
        except Exception as ex1:
            sessionOUT.rollback()
            print('this site ID already exists')
        if wf==1:
            resWF=request.get_wfdisc_rows(session,Wfdisc,sta=ssi.sta,t1=t_ini[0],t2=t_end[0],chan=res[0].chan)
            for ri in resWF:
                try:
                    sessionOUT.add(kb.Wfdisc(*ri))
                    sessionOUT.commit()
                except exc.IntegrityError as e:
                    #print 'this WF row already exists'
                    sessionOUT.rollback()

    sessionOUT.commit()

    session.close()
    sessionOUT.close()
