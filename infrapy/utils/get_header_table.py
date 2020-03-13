

import sys, pdb
import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
from obspy.core import trace
from obspy.core import Stream
from obspy.core.util import AttribDict
from datetime import datetime

import numpy as np
import pisces as ps
from pisces import request

import time


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

