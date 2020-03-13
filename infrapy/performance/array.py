from obspy.core import UTCDateTime
from obspy.core import trace
from obspy.core import Stream
from obspy.core.util import AttribDict
from obspy.signal.array_analysis import *
#import obspy.signal.util.utlGeoKm

import matplotlib.pyplot as pl

from obspy.signal.util import utlGeoKm, nextpow2
import numpy as np
from utils.cart2pol import cart2pol

import cmath
import math
import itertools

from scipy import interpolate
from scipy import signal

import warnings
warnings.filterwarnings("error")

def freqBANDS_lower(fc,order):
	return fc/np.power(2,1./order)

def freqBANDS_higher(fc,order):
	return fc*np.power(2,1./order)


def freqmodBANDS_lower(fc,order):
	G=10**0.3
	return fc/np.power(G,1./order)



def freqmodBANDS_higher(fc,order):
	return fc*np.power(G,1./order)



def freqBANDS(central,order,num):
	fc=1.0
	fca=fc
	freq=[]
	for i in range(num):
		fca=freqBANDS_higher(fca,order)
		freq.append(fca)
	fca=fc
	for i in range(num):
		fca=freqBANDS_lower(fca,order)
		freq.append(fca)
	freq=np.asarray(freq)
	return np.sort(freq)

if __name__ == "__main__":
    	print(freqBANDS(1,3,10))
