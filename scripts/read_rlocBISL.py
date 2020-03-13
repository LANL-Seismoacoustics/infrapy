

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
from code.location.assoc_lib.bisl import calc_conf_ellipse
from sqlalchemy import func

from code.location.assoc_lib.likelihood import *
from code.location.assoc_lib.propagation import *
from code.location.assoc_lib.assoc import *
from code.location.assoc_lib.bisl import *

from utils.short_time import short_time


print('change')
from utils.get_mean_locations import get_mean_locations

from code.location.assoc_lib.plotting import *
import argparse

def get_header_table(table_name):
	inf_table=[(c.name,c.info) for c in table_name.__table__.columns]
	hea=''
	for jj in inf_table:
		w=int(jj[1]['width'])
		st_header=str('%-'+str(jj[1]['width'])+'s')
		str_header_f=str(st_header%jj[0])
		if len(str_header_f)>w:
			hea=hea+' '+str_header_f[0:w-1]
		else:
			hea=hea+' '+str_header_f
	return hea

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Read Location results for specific Network")
	parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
	parser.add_argument('-n', dest='net',required=True,help="network name, e.g.: -n SMU3S")
	#parser.add_argument('-t', dest='assocresults',required=False,help="specific table with results, e.g.: -t assoc_r")
	parser.add_argument('-l','--plocid', dest='plocid',required=True,help='loc parameter id, e.g.: -l 0')
	#parser.add_argument('-i','--passocid', dest='passocid',required=True,help='loc parameter id, e.g.: -i 0')
	parser.add_argument('-o',dest='outtext',required=False,help='fd parameter id, e.g.: -o res_FILE')
	#parser.add_argument('-e','--pfdid', dest='pfdid',required=True,help='fd parameter id, e.g.: -e 0')
	parser.add_argument('-s', dest='tS',required=False,help="starttime plot, e.g.: -s /'2015-03-02T00:00:00/'")
	parser.add_argument('-e', dest='tE',required=False,help="endtime plot, e.g.: -s /'2015-03-03T00:00:00/'")
	args = parser.parse_args()

	if args.sq:
		sq=args.sq
	if args.net:
		net=args.net
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


	if args.plocid:
		plocid=int(args.plocid)

	if args.outtext:
		outtext=args.outtext
		fid=open(outtext,'wb')
	else:
		print("default table FD_RESULTS")
		outtext=None
		fid=None

	class LocParams(ab.LocParams):
		__tablename__ = 'LocParams'

	class LocResults(ab.LocResults):
		__tablename__ = 'LocResults'


	class Assoc_params(ab.assoc_params):
		__tablename__ = 'assoc_params'

	class Assoc_results(ab.assoc_results):
		__tablename__ = 'assoc_results'



	try:
		session=ps.db_connect(sq)
		if sq.count('oracle')>0:
			session_type='oracle'


			else:
				import pisces.schema.kbcore as kba

				class Site(kba.Site):
					__tablename__ = 'site'

				class Wfdisc(kba.Wfdisc):
					__tablename__ = 'wfdisc'

				class Affiliation(kba.Affiliation):
					__tablename__ = 'affiliation'

		elif sq.count('sqlite')>0:
			print('SQLITE database')
			session_type='sqlite'
			import pisces.schema.kbcore as kba
			class Site(kba.Site):
				__tablename__ = 'site'

			class Wfdisc(kba.Wfdisc):
				__tablename__ = 'wfdisc'
			class Affiliation(kba.Affiliation):
				__tablename__ = 'affiliation'

		else:
			print('No standard database, exiting')
			exit()

			session_type='unknown'
	except Exception as ex1:
		print('Connection failed:', ex1)
		embed()
		sys.exit()
	if t_S==None:
		loc_res=session.query(LocResults).filter(LocResults.eventid>=0).filter(LocResults.plocid==plocid).order_by(LocResults.timeorigmean).all()
		t_S=session.query(func.min(LocResults.timeini)).filter(LocResults.eventid>=0).filter(LocResults.plocid==plocid).all()

	else:
		loc_res=session.query(LocResults).filter(LocResults.eventid>=0).filter(LocResults.plocid==plocid).order_by(LocResults.timeorigmean)\
		.filter(LocResults.timeini>=float(t_S)).filter(LocResults.timeini<=float(t_E)).all()
	#loc_res=session.query(LocResults).filter(LocResults.eventid>=0).filter(LocResults.plocid==plocid).order_by(LocResults.timeorigmean).all()



	fig1=pl.figure()
	ax1=pl.subplot(1,1,1)
	loc_res_clean=[]
	for ll in loc_res:
		if ll.eventid==-1 or ll.numstations<2 or ll.latorigmean< -180:
			#print ll
			continue
		try:
			conflat,conflon=calc_conf_ellipse([ll.latorigmean,ll.lonorigmean],[ll.latorigvar,ll.lonorigvar,ll.latlonorigcovar],95.0)
		except Exception as ex1:
			print('error boundary calc:',ex1)
			#embed()

		if ll.numstations==2:
			ax1.plot(ll.lonorigmean,ll.latorigmean,'bo',markeredgecolor='b')
			ax1.plot(conflon,conflat,'b')
		else:
			ax1.plot(ll.lonorigmean,ll.latorigmean,'mo',markeredgecolor='m')
			ax1.plot(conflon,conflat,'m')

		loc_res_clean.append(ll)

	ax1.set_xlim([-120,-108])
	ax1.set_ylim([36,46])
	ax1.grid('on')
	ax1.set_title('fig_ploc'+str(plocid))
	fig1.savefig('fig_ploc'+str(plocid)+'.pdf')
	pl.show()

	exit()
	dict_namefd={}
	dict_namefk={}

	for ll in loc_res_clean:
		f = open('EE_'+str(ll.eventid)+'.txt', 'w')
		strAU='loc num:'+str(ll.eventid)
		print(strAU)
		f.write(strAU+'\n')

		strAU='***** Loc results   *****'
		print(strAU)
		f.write(str(strAU+'\n'))

		strAU=get_header_table(LocResults)
		print(strAU)
		f.write(strAU+'\n')

		strAU=str(ll)
		print(strAU)
		f.write(strAU+'\n')

		assoc_res=session.query(Assoc_results).filter(Assoc_results.eventid==ll.eventid).order_by(Assoc_results.eventid).all()
		det_tot=[]
		for tn in assoc_res:

			strAU='***** Assoc  *****'
			print(strAU)
			f.write(str(strAU)+'\n')

			strAU=get_header_table(Assoc_results)
			print(strAU)
			f.write(strAU+'\n')

			if tn.fdid<0:
				continue

			strAU=tn
			print(strAU)
			f.write(str(strAU)+'\n')

			if not(tn.fdtable in list(dict_namefd.keys())):
				dict_namefd[tn.fdtable]=type(tn.fdtable.encode('utf-8') ,(ab.fd_results,),{'__tablename__':tn.fdtable.encode('utf-8')})
			det_tab=dict_namefd[tn.fdtable]
			#embed()
			fd_res=session.query(det_tab).filter(det_tab.fdid==tn.fdid).one()
			strAU='***** Detection  *****'
			print(strAU)
			f.write(str(strAU)+'\n')

			strAU=get_header_table(det_tab)
			print(strAU)
			f.write(strAU+'\n')

			strAU=str(fd_res)
			print(strAU)
			f.write(strAU+'\n')

			if not(fd_res.fktablename in list(dict_namefk.keys())):
				dict_namefk[fd_res.fktablename]=type(fd_res.fktablename.encode('utf-8') ,(ab.fk_results,),{'__tablename__':fd_res.fktablename.encode('utf-8')})

			fk_tab=dict_namefk[fd_res.fktablename]

			Query_fk=session.query(fk_tab).filter(fk_tab.timeini>=fd_res.timeini).filter(fk_tab.timeini<=fd_res.timeend).all()
			strAU='***** FK  *****'
			print(strAU)
			f.write(str(strAU)+'\n')

			strAU=get_header_table(fk_tab)
			print(strAU)
			f.write(strAU+'\n')

			for ff in Query_fk:
				strAU=str(ff)
				print(strAU)
				f.write(strAU+'\n')

			staFK_all=session.query(Site).filter(Site.refsta==Query_fk[0].sta).all()
			staFK=staFK_all[0]
			det1 = inf_det_global(staFK.lat, staFK.lon, UTCDateTime(Query_fk[0].timeini).datetime, Query_fk[0].bz,  Query_fk[0].fval, 1,DetID=staFK.sta)
			det_tot.append(det1)
		#embed()
		f.close()
		#exit()

