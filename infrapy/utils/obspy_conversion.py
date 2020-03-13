import numpy as np
from obspy.core import UTCDateTime

from obspy.core import trace
from obspy.core import Stream
from obspy.core.util import AttribDict



def stream2array(str):
	""" turn obspy stream into 2D array

	the data from stream is extracted and copied to a 2D numpy array
	 - ** parameters**, **types**, **return**, and **return types**::
	 	:param arg1: stream to be converted
		:type arg1: obspy stream
		:return: 2D array with time series for all channels
		:rtype: numpy array
	"""

	tI=[]
	tE=[]
	for sti in str:
		tI.append(float(sti.stats.starttime))
		tE.append(float(sti.stats.endtime))
	try:
		tIMAX=np.max(np.asarray(tI))
		tEMIN=np.min(np.asarray(tE))
	except Exception as ex1:
		print('error at stream2array, FK:',ex1)
		pdb.set_trace()
		exit()
	str.trim(starttime=UTCDateTime(tIMAX),endtime=UTCDateTime(tEMIN))
	num_samples=[]
	for sti in str:
		num_samples.append(sti.stats.npts)
	if len(np.unique(num_samples))>1:
		print('there is an issue with the number of samples,')
		print('maybe different sampling rates')
		return []
	else:
		for st_i in range(len(str)):
			if str[st_i].stats.calib==0:
				print("there is an error with calib number")
				print('setting calib to 1')
				str[st_i].stats.calib=1
			if st_i==0:
				mat=str[st_i].data*str[st_i].stats.calib
			else:
				mat=np.vstack((mat,str[st_i].data*str[st_i].stats.calib))
	return mat
