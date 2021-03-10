# infrapy.performance.network.py
#
# Methods to estimate network performance using
# the various propagation models included in the
# infrapy toolkit.
#
# Author            Philip Blom (pblom@lanl.gov)


import itertools
import imp 

import numpy as np

from scipy.interpolate import interp1d
from scipy.integrate import simps

from ..characterization import spye

from pyproj import Geod


sph_proj = Geod(ellps='sphere')
ref_dB = 10.0 * np.log10(20.e-6)

# ########################## #
#      IMS Noise Models      #
# ########################## #
ims_low_ns = np.loadtxt(imp.find_module('infrapy')[1] + '/resources/noise_models/IMSnoisemodel_low.txt')
ims_high_ns = np.loadtxt(imp.find_module('infrapy')[1] + '/resources/noise_models/IMSnoisemodel_high.txt')

low_ns_interp = interp1d(ims_low_ns[:, 0], 10.0 * ims_low_ns[:,1] - 2.0 * ref_dB)
high_ns_interp = interp1d(ims_high_ns[:, 0], 10.0 * ims_high_ns[:,1] - 2.0 * ref_dB)

def low_noise(freq):
    """
        Interpolation of the IMS low noise PSD curve
                
        Parameters
        ----------
        freq: float
            Frequency [Hz]
            
        Returns
        -------
        Noise PSD: float
            Power spectral density of the noise [Pa^2/Hz]
    """

    return low_ns_interp(freq)


def high_noise(freq):
    """
        Interpolation of the IMS high noise PSD curve
                
        Parameters
        ----------
        freq: float
            Frequency [Hz]
            
        Returns
        -------
        Noise PSD: float
            Power spectral density of the noise [Pa^2/Hz]
    """

    return high_ns_interp(freq)


def med_noise(freq):
    """
        Mean of the IMS high and low noise PSD curves
                
        Parameters
        ----------
        freq: float
            Frequency [Hz]
            
        Returns
        -------
        Noise PSD: float
            Power spectral density of the noise [Pa^2/Hz]
    """
    return (low_ns_interp(freq) + high_ns_interp(freq)) / 2.0

# ########################## #
#    Detection Probability   #
#        Calculations        #
# ########################## #
def signal_duration(rng, cel_lims=[0.28, 0.34]):
    """
        Compute the duration of an analysis window for tropospheri
        and stratospheric returns
                
        Parameters
        ----------
        rng: float
            Propagation Range [km]
        cel_lims: iterable
            Min and max celerity to define the window (defaults to 280 and 340 m/s)
            
        Returns
        -------
        duration: float
            Duration of the analysis window
    """
    return max(5.0, rng * (cel_lims[1] - cel_lims[0]) / (cel_lims[0] * cel_lims[1]))


def loc_prob(rng, az, tloss_model, f0, W=100.0e3, noise_lvl="medium", resol=500):
    dur = signal_duration(rng)

    r0 = 0.035 * W**(1.0 / 3.0)
    src_amp = 10.0 * np.log10(spye.blastwave_spectrum(f0, W, r0) * r0) - ref_dB

    tloss_vals =  np.linspace(10.0 * np.log10(1.0 / rng**2.0), 10.0 * np.log10(1.0 / np.sqrt(rng)), resol)
    tloss_pdf = tloss_model.eval(np.array([rng] * len(tloss_vals)), tloss_vals, np.array([az] * len(tloss_vals)))
    
    if noise_lvl=="high":
        mask = (src_amp + tloss_vals) > high_noise(f0) / 2.0 + 5.0 * np.log10(dur) 
    elif noise_lvl == "medium":
        mask = (src_amp + tloss_vals) > med_noise(f0) / 2.0 + 5.0 * np.log10(dur) 
    elif noise_lvl == "low":
        mask = (src_amp + tloss_vals) > low_noise(f0) / 2.0 + 5.0 * np.log10(dur) 
    else:
        print("warning...")

    if len(tloss_pdf[mask]) > 2:
        return simps(tloss_pdf[mask], src_amp + tloss_vals[mask])
    else:
        return 0.0


def _loc_prob_wrapper(args):
    return loc_prob(*args)


def detection_prob(src_locs, array_locs, tloss_model, f0, W=100.0e3, noise_lvl="medium", resol=500, det_cnt_min=2, pool=None):
    src_locs = np.atleast_2d(src_locs)
    LA1, LA2 = np.meshgrid(src_locs[:, 0], array_locs[:, 0])
    LO1, LO2 = np.meshgrid(src_locs[:, 1], array_locs[:, 1])

    SRC_LATS, SRC_LONS = LA1.flatten(), LO1.flatten()
    RCVR_LATS, RCVR_LONS = LA2.flatten(), LO2.flatten()

    temp = np.array(sph_proj.inv(SRC_LONS, SRC_LATS, RCVR_LONS, RCVR_LATS, radians=False)).T
    azs, rngs = temp[:, 0], temp[:, 2] / 1000.0
    rngs[rngs < 10.0] = 10.0

    if pool:
        args = [[rngs[n], azs[n], tloss_model, f0, W, noise_lvl, resol] for n in range(len(rngs))]
        probs = np.array(pool.map(_loc_prob_wrapper, args))
    else:
        probs = np.array([loc_prob(rngs[n], azs[n], tloss_model, f0, W=W, noise_lvl=noise_lvl, resol=resol) for n in range(len(rngs))])
    
    probs = probs.reshape(len(array_locs), len(src_locs)).T

    array_cnt = len(array_locs)
    result = 1.0 - np.prod(1.0 - probs, axis=1)
    for n in range(1, det_cnt_min):
        for indices in itertools.combinations(list(range(array_cnt)), r=n):
            mask = np.ones_like(probs, dtype=bool)            
            mask[:, list(indices)] = False

            P1, P2 = 1.0 * probs, 1.0 - np.copy(probs)
            P1[mask] = 1.0
            P2[np.invert(mask)] = 1.0

            result = result - np.prod(P1, axis=1) * np.prod(P2, axis=1)

    return result