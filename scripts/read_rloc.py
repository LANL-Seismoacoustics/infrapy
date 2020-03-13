

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


import schema.infrapy_schema_abs as ab
from utils.short_time import short_time

import argparse


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Read Association results for specific Network")
	parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
	parser.add_argument('-n', dest='net',required=True,help="network name, e.g.: -n SMU3S")
	parser.add_argument('-t', dest='assocresults',required=False,help="specific table with results, e.g.: -t assoc_r")
	parser.add_argument('-i','--passocid', dest='passocid',required=True,help='assoc parameter id, e.g.: -i 0')
	parser.add_argument('-o',dest='outtext',required=False,help='fd parameter id, e.g.: -o res_FILE')
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
		fid=open(outtext,'wb')
	else:
		print("default table FD_RESULTS")
		outtext=None
		fid=None

	
	class Assoc_params(ab.assoc_params):
		__tablename__ = 'assoc_params'
		
	AssocT=type(assocresultsT,(ab.assoc_results,),{'__tablename__':assocresultsT})
	'''
	class Assoc_results(ab.assoc_results):
		__tablename__ = assocresultsT
	'''	

	try:
		session=ps.db_connect(sq)
		if sq.count('oracle')>0:
			session_type='oracle'

			else:
				import pisces.schema.kbcore as kba
				'''
				class Site(kba.Site):
					__tablename__ = 'site'

				class Wfdisc(kba.Wfdisc):
					__tablename__ = 'wfdisc'
				'''
		elif sq.count('sqlite')>0:
			print('SQLITE database')
			session_type='sqlite'
			import pisces.schema.kbcore as kba
			'''
			class Site(kba.Site):
				__tablename__ = 'site'

			class Wfdisc(kba.Wfdisc):
				__tablename__ = 'wfdisc'
			'''
		else:
			print('No standard database, exiting')
			exit()

			session_type='unknown'
	except Exception as ex1:
		print('Connection failed:', ex1)
		embed()
		sys.exit()


	assoc_res=session.query(AssocT).filter(AssocT.passocid==passocid).filter(AssocT.net==net).filter(AssocT.qassoc>0).order_by(AssocT.associd).all()
	assoc_res_table_names=session.query(AssocT.fdtable).filter(AssocT.passocid==passocid).filter(AssocT.net==net).filter(AssocT.qassoc>0).distinct().all()
	assoc_res_eventid=session.query(AssocT.eventid).filter(AssocT.passocid==passocid).filter(AssocT.net==net).filter(AssocT.qassoc>0).order_by(AssocT.eventid).distinct().all()
	
	dict_namefd={}
	dict_namefk={}
	for tn in assoc_res_table_names:
		try:
			dict_namefd[tn[0]]=type(tn[0].encode('utf-8') ,(ab.fd_results,),{'__tablename__':tn[0].encode('utf-8')})
		except Exception as ex1:
			embed()
		temp_name=[]
		temp_name=dict_namefd[tn[0]]
		fk_table_names=session.query(temp_name.fktablename).distinct().all()
		#print fk_table_names
		for tnfk in fk_table_names:
			try:
				dict_namefk[tnfk[0]]=type(tnfk[0].encode('utf-8') ,(ab.fk_results,),{'__tablename__':tnfk[0].encode('utf-8')})
			except Exception as ex1:
				print(ex1,'130')

				
	if  len(assoc_res)==0:
		print('No results with passocid:', passocid)
		exit()
	print(' ')	
	for evid in assoc_res_eventid:
			
		#embed()
		assoc_aux=session.query(AssocT).filter(AssocT.passocid==passocid).filter(AssocT.net==net).filter(AssocT.eventid==evid[0]).filter(AssocT.qassoc>0).order_by(AssocT.associd).all()			
		stas=[]
		
		if assoc_aux[0].qassoc<1E-4:
			print('low assoc')
			continue
		
		
		
		for ai in assoc_aux:
			try:
				fd_res=session.query(dict_namefd[ai.fdtable]).filter(dict_namefd[ai.fdtable].fdid==ai.fdid).one()
			except Exception as ex1:
				print('problem with table')
				embed()
			stas.append(fd_res.sta)	

		if len(list(set(stas)))<3:
			continue
		
		for ai in assoc_aux:
			try:
				fd_res=session.query(dict_namefd[ai.fdtable]).filter(dict_namefd[ai.fdtable].fdid==ai.fdid).one()
			except Exception as ex1:
				print('problem with table')
			print(short_time(UTCDateTime(fd_res.timeini)), int(fd_res.timeend-fd_res.timeini),fd_res.sta, ai.eventid, end=' ') 
			temp=[]
			temp=dict_namefk[fd_res.fktablename]
			#Query_fktemp=session.query(temp).filter(temp.pfkid==fd_res.pfkid).filter(temp.timeini>=float(fd_res.timeini)).filter(temp.timeend<=float(fd_res.timeend)).filter(temp.sta==fd_res.sta)
			Query_fktempM=session.query(temp).filter(temp.pfkid==fd_res.pfkid).filter(temp.timeini==float(fd_res.maxf_time)).filter(temp.sta==fd_res.sta)

			for qq1 in Query_fktempM:
					print(int(qq1.bz),ai.qassoc)
					break
			if not (fid==None):
				STRD=short_time(UTCDateTime(fd_res.timeini))+' '+str(int(fd_res.timeend-fd_res.timeini))+' '+fd_res.sta+' '+str(ai.eventid)+' '+str(int(qq1.bz))+' '+str(ai.qassoc)+'\n'
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

