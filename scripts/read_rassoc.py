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
		assocresultsT='ASSOC_results'

	if args.passocid:
		passocid=int(args.passocid)


	class Assoc_params(schema.ASSOC_params):
		__tablename__ = 'ASSOC_params'

	AssocT=type(assocresultsT,(schema.ASSOC_results,),{'__tablename__':assocresultsT})
	'''
	class Assoc_results(ab.assoc_results):
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


	assoc_res=session.query(AssocT).filter(AssocT.passocid==passocid).filter(AssocT.net==net).filter(AssocT.qassoc>0).order_by(AssocT.associd).all()
	assoc_res_table_names=session.query(AssocT.fdtable).filter(AssocT.passocid==passocid).filter(AssocT.net==net).filter(AssocT.qassoc>0).distinct().all()
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
		print(fk_table_names)
		for tnfk in fk_table_names:
			try:
				dict_namefk[tnfk[0]]=type(tnfk[0].encode('utf-8') ,(ab.fk_results,),{'__tablename__':tnfk[0].encode('utf-8')})
			except Exception as ex1:
				print(ex1,'130')

		'''

		Query_fktemp=self.session.query(fktable).filter(fktable.pfkid==self.pfkid).filter(fktable.timeini>=float(dqi.timeini)).filter(fktable.timeend<=float(dqi.timeend)).filter(fktable.sta==aai['name'])

		Query_fktempMAX=self.session.query(fktable).filter(fktable.pfkid==self.pfkid).filter(fktable.timeini==float(dqi.maxf_time)).filter(fktable.sta==aai['name']).one()
		'''


	if  len(assoc_res)==0:
		print('No results with passocid:', passocid)
		exit()


	for ai in assoc_res:
		#print ai
		try:
			fd_res=session.query(dict_namefd[ai.fdtable]).filter(dict_namefd[ai.fdtable].fdid==ai.fdid).one()
		except Exception as ex1:
			embed()
		print(short_time(UTCDateTime(fd_res.timeini)), int(fd_res.timeend-fd_res.timeini),fd_res.sta, ai.eventid, end=' ')
		temp=[]
		temp=dict_namefk[fd_res.fktablename]
		Query_fktemp=session.query(temp).filter(temp.pfkid==fd_res.pfkid).filter(temp.timeini>=float(fd_res.timeini)).filter(temp.timeend<=float(fd_res.timeend)).filter(temp.sta==fd_res.sta)
		Query_fktempMAX=session.query(temp).filter(temp.pfkid==fd_res.pfkid).filter(temp.timeini==float(fd_res.maxf_time)).filter(temp.sta==fd_res.sta)
		for qq in Query_fktempMAX:
			print(qq.bz)

	sys.exit(0)
