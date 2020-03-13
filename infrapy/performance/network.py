# infrapy.performance.network.py
#
# Methods to estimate network performance using
# the various propagation models included in the
# infrapy toolkit.
#
# Author            Philip Blom (pblom@lanl.gov)

import sys
import datetime
import time
import itertools
import random
import copyreg

import numpy as np

from scipy.integrate import quad
from scipy.interpolate import interp1d

from pyproj import Geod

from ..propagation import likelihoods as lklhds
from ..utils import prog_bar
from ..utils import latlon as ll

# ################################ #
#    Set Integration Parameters    #
#       and WGS84 Ellipsoid        #
# ################################ #
int_opts = {'limit': 100, 'epsrel': 1.0e-3}
wgs84 = Geod(ellps='WGS84')

def tloss_ccdf(tloss_distribution, threshold):
    return quad(tloss_distribution, threshold, 0.0, limit=100, epsrel= 1.0e-3)[0]


def detection_probability(lat, lon, nodes, model, det_cnt_min=2):
    node_cnt = len(nodes)

    probs = np.empty(node_cnt)
    for n in range(node_cnt):
        temp = wgs84.inv(lon, lat, nodes[n][1], nodes[n][0], radians=False)
        az, rng = temp[0], temp[2] / 1000.0
        def node_pdf(tloss):
            return model.eval(rng, tloss, az)
        probs[n] = tloss_ccdf(node_pdf, nodes[n][2])

    result = 1.0 - np.prod(1.0 - probs)
    for n in range(1, det_cnt_min):
        for indices in itertools.combinations(list(range(node_cnt)), r=n):
            mask = np.ones(node_cnt, dtype=bool)
            mask[list(indices)] = False
            result -= np.prod(probs[np.invert(mask)]) * np.prod(1.0 - probs[mask])

    return result


def detection_probability_wrapper(args):
    return detection_probability(*args)


