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
import argparse
from IPython import embed

class Fd_params(schema.fd_params):
	__tablename__ = 'FD_PARAMS'


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description="Print combination of parameters used for FD processing")
	parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
	args = parser.parse_args()
	if args.sq:
		sq=args.sq
	session=ps.db_connect(sq)
	fd_par=session.query(Fd_params).all()
	#print(Fd_params.__doc__)
	for fd_i in fd_par:
		print('pfdid:',str(fd_i.pfdid),'p-thresh:',str(fd_i.pthreshold),'min event len',str(fd_i.minlen),'detection win len:',str(fd_i.detwinlen))
	if len(fd_par)==0:
		print('No F Detection parameters in the database')
