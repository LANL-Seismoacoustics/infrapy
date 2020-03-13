#!/usr/bin/env python

from obspy.core import read
from obspy.core import Stream
import matplotlib

matplotlib.use('TkAgg')
from pylab import *
from scipy import stats
import numpy as np
from pisces.io.trace import read_waveform
from obspy.core import UTCDateTime
from obspy.core import trace
from obspy.core import Stream
from obspy.core.util import AttribDict
from obspy.signal.array_analysis import *
from obspy.core import read, Trace, Stream, UTCDateTime
from obspy.core import Trace
import matplotlib.pyplot as pl
import matplotlib.mlab as mpy
import pisces as ps


from scipy import signal

from datetime import datetime
import glob

#from pisces.schema.css3 import Base
import pisces.schema.kbcore as kba
import os
from IPython import embed
'''
BUK2 =[68.8803,15.5963,]
BUK3 =[68.8804,15.6283,]
BUK1 =[68.8887,15.5822,]

SKD2 =[69.2288,16.0271]
SKD3 =[69.2213, 16.0371]
SKD4 =[69.2307,16.0403]
'''


pl.ioff()

if __name__ == '__main__':

	ref=sys.argv[1]
	ARRAY_NAME=sys.argv[2]
	loc=sys.argv[3]
	try:
		CHAN=sys.argv[4]
	except Exception as e:
		print('no chan supplied using channel 1')
		CHAN=str(1)
		
	try:
		DAS=sys.argv[5]
	except Exception as e:
		print('no chan supplied using channel 1')
		DAS=None

	base=os.getcwd()

	try:
		os.mkdir(ARRAY_NAME+'_SAC')
	except:
		print('already in the folder')
	os.chdir(ARRAY_NAME+'_SAC')
	print('converting reftek to ms')

	if DAS is None:
		os.system("/opt/passcal/bin/rt2ms"+" -f ../" + ref+" -s "+ARRAY_NAME+" --channelskeep "+CHAN+" -R -Y ")
	else:
		print('converting das:'+DAS)
		os.system("/opt/passcal/bin/rt2ms"+" -f ../" + ref+" -s "+ARRAY_NAME+" --channelskeep "+CHAN+" -R -Y  --daskeep "+DAS)

	os.chdir(base)
	mfolders = glob.glob(ARRAY_NAME+'_SAC/R*')
	if len(mfolders)==0:
		mfolders = glob.glob(ARRAY_NAME+'_SAC/Y*/R*')
		#embed()

	print('converting ms to sac')
	for foli in mfolders:
		mfolders = glob.glob(base+'/'+foli+'/*.'+CHAN+'.m')
		if len(mfolders)==0:
			mfolders = glob.glob(base+'/'+foli+'/*.'+CHAN)

		aT=Stream()
		for mi in mfolders:

			os.system("/opt/passcal/bin/ms2sac -c "+loc+" "+mi+" "+mi+".sac" )
			fl=mi.split('/')
			fi=fl[-1]
			fii=fi.split('.')
			try:
				aa=read(mi+".sac",format="SAC")
			except Exception as ex1:
				print('error reading', ex1)
				embed()
				sys.exit()

			aa.merge()
			aT.append(aa[0])
			os.remove(mi)
			os.remove(mi+".sac")
		aT.merge()
		try:
			aT[0].stats.station=ARRAY_NAME
		except Exception as ex1:
			print('error',ex1)
			embed()
			sys.exit()


		try:
			#aT.write(base+'/'+ARRAY_NAME+'_SAC/'+ARRAY_NAME+'_'+fii[1]+'_'+fii[0]+'.sac', format='SAC')
			aT.write(base+'/'+ARRAY_NAME+'_'+fii[1]+'_'+fii[0]+'.sac', format='SAC')
		except Exception as ex1:
			print('issue:',ex1)
			embed()
			for tr in aT:
				if isinstance(tr.data, np.ma.masked_array):
					tr.data = tr.data.filled(fill_value=0)
			aT.write(base+'/'+ARRAY_NAME+'_SAC/'+ARRAY_NAME+'_'+fii[1]+'_'+fii[0]+'.sac', format='SAC')

			#continue
	pass
