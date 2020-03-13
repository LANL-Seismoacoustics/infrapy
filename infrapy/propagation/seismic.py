# infrapy.propagation.seismic.py
#
# Seismic propagation models for association and localization.
#
# Author            Philip Blom (pblom@lanl.gov)

import sys
import pickle
import imp
import time
import itertools

import numpy as np

from pyproj import Geod

from scipy.interpolate import interp1d

np.seterr(over='ignore')
sph_proj = Geod(ellps='sphere')

# ######################### #
#     ak135 Travel Time     #
#           Tables          #
# ######################### #

# phase specific variances
ak135_p_sigma = 3.0
ak135_s_sigma = 6.0

# load ak135 travel time tabl and interpolate
# file format: arc [deg] : p_time [s] : s_time [s]
ak135_tbls = np.loadtxt(imp.find_module('infrapy')[1] + '/resources/travelTimeTables/ak135_1st_arrivals.dat', skiprows=1)

ak135_p_tr_time_interp = interp1d(np.radians(ak135_tbls[:, 0]) * 6370.997, ak135_tbls[:, 1], kind='cubic')
ak135_s_tr_time_interp = interp1d(np.radians(ak135_tbls[:, 0]) * 6370.997, ak135_tbls[:, 2], kind='cubic')

def ak135_p_tr_time(rng):
    return ak135_p_tr_time_interp(rng)

def ak135_s_tr_time(rng):
    return ak135_s_tr_time_interp(rng)
