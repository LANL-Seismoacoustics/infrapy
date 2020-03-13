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

if __name__ == '__main__':

	sq=sys.argv[1]
	array=sys.argv[2]

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
		sys.exit()
	#embed()
	resSite=session.query(Site).filter(Site.sta.like('%'+array+'%')).all()
	#embed()

	for rS_i in resSite:
		rS_i.refsta=array
		session.commit()

	resWfdisc=session.query(Wfdisc).filter(Wfdisc.sta.like('%'+array+'%')).all()
	#embed()

	for rS_i in resWfdisc:
		rS_i.datatype='sc'
		session.commit()
