

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

import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read Association results for specific Network")
    parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
    parser.add_argument('-n', dest='net',required=True,help="network name, e.g.: -n SMU3S")
    parser.add_argument('-t', dest='assocresults',required=False,help="specific table with results, e.g.: -t assoc_r")
    parser.add_argument('-i','--passocid', dest='passocid',required=True,help='assoc parameter id, e.g.: -i 0')
    parser.add_argument('-o',dest='outtext',required=False,help='outtext file name, e.g.: -o res_FILE')
    #parser.add_argument('-e','--pfdid', dest='pfdid',required=True,help='fd parameter id, e.g.: -e 0')
    args = parser.parse_args()

    if args.sq:
        sq=args.sq
    if args.net:
        net=args.net

    if args.assocresults:
        assocresultsT=args.assocresults
    else:
        print("default table FD_RESULTS")
        assocresultsT='assoc_results'

    if args.passocid:
        passocid=int(args.passocid)


    if args.outtext:
        outtext=args.outtext
        fid=open(outtext,'w')
    else:
        print("default table FD_RESULTS")
        outtext=None
        fid=None

    class Assoc_params(schema.ASSOC_params):
        __tablename__ = 'assoc_params'

    AssocT=type(assocresultsT,(schema.ASSOC_results,),{'__tablename__':assocresultsT})
    '''
    class Assoc_results(schema.assoc_results):
        __tablename__ = assocresultsT
    '''

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

    '''

        refSTA=[]
        for aai in self.Affiliation_Q:

            try:
                STA_dataM=self.session.query(self.Site).filter(self.Site.sta==aai.sta).one()

            except Exception as ex1:
                #print
                print('there is more than just one station:', aai.sta,'  ',ex1)
                embed()
                exit()
            #embed()
            refSTA.append(STA_dataM.refsta)

        refstations_l=list(set(refSTA))
        refsta=[]
        #embed()
        for aai in refstations_l:
            STA_dataM=self.session.query(self.Site).filter(self.Site.refsta==str(aai)).all()
            array_lo=[]
            array_la=[]
            array_el=[]
            for sta_i in STA_dataM:
                array_la.append(sta_i.lat)
                array_lo.append(sta_i.lon)
                array_el.append(sta_i.elev)
            array_la=np.asarray(array_la)
            array_lo=np.asarray(array_lo)
            array_el=np.asarray(array_el)
            refsta.append({'lon':np.mean(array_lo),'lat':np.mean(array_la),'elev':np.mean(array_el),'name':aai,'numsta':len(array_la)})
        self.det_tot=[]
        self.fdtable_name=[]
        self.fktables_names=[]
        for aai in refsta:
            #print('getting data from:',aai['name'])
            # here it looks for any detections in all the tables
            for ti in range(self.num_tables):
                try:

                    fd_res=self.session.query(self.fdtables[ti]).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).all()
                    times_ini=self.session.query(self.fdtables[ti].timeini).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).all()
                    fk_table_names=self.session.query(self.fdtables[ti].fktablename).filter(self.fdtables[ti].sta==aai['name']).filter(self.fdtables[ti].pfdid==self.pfdid).filter(self.fdtables[ti].pfkid==self.pfkid).all()
                    #embed()
                    if len(fk_table_names)>0:
                        for tnfk in fk_table_names:
                            self.fktables_names.append(tnfk[0])

                            #embed()

                except Exception as x1:
                    print('There is an error',x1)
                    embed()
                    exit()

        tt = np.unique(self.fktables_names)
        self.fktables=[]
        for  fki in tt:
            self.fktables.append(type(str(fki),(schema.fk_results,),{'__tablename__':str(fki)}))

        return self.db_connected
    '''

    assoc_res=session.query(AssocT).filter(AssocT.passocid==passocid).filter(AssocT.net==net).order_by(AssocT.associd).all()
    assoc_res_table_names=session.query(AssocT.fdtable).filter(AssocT.passocid==passocid).filter(AssocT.net==net).distinct().all()
    assoc_res_eventid=session.query(AssocT.eventid).filter(AssocT.passocid==passocid).filter(AssocT.net==net).order_by(AssocT.eventid).distinct().all()

    tt = np.unique(assoc_res_table_names)
    fdtables=[]
    fd_table_names=[]
    fk_table_names=[]
    fktables=[]
    for fdi in tt:
        fdtables.append(type(str(fdi),(schema.fd_results,),{'__tablename__':str(fdi)}))
        fd_table_names.append(fdi)

    if  len(assoc_res)==0:
        print('No results with passocid:', passocid)
        exit()
    print(' ')
    for evid in assoc_res_eventid:

        #embed()
        assoc_aux=session.query(AssocT).filter(AssocT.passocid==passocid).filter(AssocT.net==net).filter(AssocT.eventid==evid[0]).order_by(AssocT.associd).all()

        for ai in assoc_aux:
            try:
                for aa in range(len(fd_table_names)):
                    if fd_table_names[aa]==ai.fdtable:
                        fd_res=session.query(fdtables[aa]).filter(fdtables[aa].fdid==ai.fdid).one()
            except Exception as ex1:
                print('problem with table')
                embed()
            fk_table_names.append(fd_res.fktablename)

    tt = np.unique(fk_table_names)
    fktables=[]
    for  fki in tt:
        fktables.append(type(str(fki),(schema.fk_results,),{'__tablename__':str(fki)}))

    for evid in assoc_res_eventid:

        #embed()
        assoc_aux=session.query(AssocT).filter(AssocT.passocid==passocid).filter(AssocT.net==net).filter(AssocT.eventid==evid[0]).order_by(AssocT.associd).all()
        for ai in assoc_aux:
            try:
                for aa in range(len(fd_table_names)):
                    if fd_table_names[aa]==ai.fdtable:
                        fd_res=session.query(fdtables[aa]).filter(fdtables[aa].fdid==ai.fdid).one()
            except Exception as ex1:
                print('problem with detection table')

            print(short_time(UTCDateTime(fd_res.timeini)), int(fd_res.timeend-fd_res.timeini),fd_res.sta, ai.eventid, end=' ')
            try:
                for bb in range(len(tt)):
                    if tt[bb]==fd_res.fktablename:

                        fk_res=session.query(fktables[bb]).filter(fktables[bb].pfkid==fd_res.pfkid).filter(fktables[bb].timeini==float(fd_res.timeini)).filter(fktables[bb].sta==fd_res.sta).all()
            except Exception as ex1:
                print('problem with fk table')
                embed()
            if not (fid==None):
                STRD=str(ai.eventid)+' '+fd_res.sta+' '+str(fd_res.fdid)+' '+short_time(UTCDateTime(fd_res.timeini))+' '+str(int(fd_res.timeend-fd_res.timeini))+' '+str(int(fk_res[0].bz))+' '+str(ai.qdetcluster)+'\n'
                fid.write(STRD)
            '''
            Query_fktempMAX=session.query(temp).filter(temp.pfkid==fd_res.pfkid).filter(temp.timeini==float(fd_res.maxf_time)).filter(temp.sta==fd_res.sta)
            for qq in Query_fktempMAX:
            '''
            #embed()
    fid.close()
            #exit()
        #exit()
    sys.exit(0)
