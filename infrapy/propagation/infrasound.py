# infrapy.propagation.infrasound.py
#
# Infrasound propagation models for association and localization.
# General models are valid at all locations.  Propagation-based
# stochastic models are location, month, and time of day specific
#
# Author            Philip Blom (pblom@lanl.gov)

import warnings
import json

import numpy as np

from pathlib import Path

from scipy.stats import norm

from ..utils import prog_bar
from ..utils import skew_norm

np.seterr(over='ignore', divide='ignore')
warnings.simplefilter(action='ignore', category=FutureWarning)

# ############################ #
#      General (Canonical)     #
#      Propagation Models      #
# ############################ #
canon_rcel_wts = np.array([0.0539, 0.0899, 0.8562])
canon_rcel_mns = np.array([1.0 / 0.327, 1.0 / 0.293, 1.0 / 0.26])
canon_rcel_vrs = np.array([0.066, 0.08, 0.33])

def canonical_rcel(rcel):
    if len(np.atleast_1d(rcel)) == 1:
        vals = canon_rcel_wts / canon_rcel_vrs * norm.pdf((rcel - canon_rcel_mns) / canon_rcel_vrs)
        return np.sum(vals)
    else:
        vals = np.asarray([canon_rcel_wts] * len(rcel)) / np.asarray([canon_rcel_vrs] * len(rcel)) * norm.pdf((np.asarray([rcel] * 3).T - np.asarray([canon_rcel_mns] * len(rcel))) / np.asarray([canon_rcel_vrs] * len(rcel)))
        return np.sum(vals, axis=1)


# ########################### #
#     Accessible Celerity     #
#      Statistics Models      #
# ########################### #
def _load_celerity_model(option):

    global canon_rcel_wts
    global canon_rcel_mns
    global canon_rcel_vrs

    if option in ["regional_hf", "regional_lf", "infGEM"]:
        rcel_gmm_file = str(Path(__file__).parent.parent)
        rcel_gmm_file = rcel_gmm_file + "/resources/travelTimeTables/" +  option + ".rcg.json"
    else: 
        rcel_gmm_file = option 

    with open(rcel_gmm_file, 'r') as infile:
        rcel_gmm = json.load(infile)
        canon_rcel_wts = np.array(rcel_gmm["weights"])
        canon_rcel_mns = np.array(rcel_gmm["means"])
        canon_rcel_vrs = np.array(rcel_gmm["stdevs"])


canon_tloss_rates = np.array([-0.87, -0.835, -0.81])
canon_tloss_shifts = np.array([1.5, 1.5, 1.5])
canon_tloss_widths = np.array([7.5, 7.5, 7.5])
canon_tloss_skews = np.array([-6.6, -6.6, -6.6])
canon_tloss_weights = np.array([0.0, 0.5, 0.0])

def canonical_tloss(rng, tloss):
    rng = max(rng, 0.01)
    widths = canon_tloss_widths + (0.33 - canon_tloss_widths) * np.exp(-rng / 25.0)
    vals = canon_tloss_weights * skew_norm.pdf(tloss, 20.0 * np.log10((rng)**canon_tloss_rates) + canon_tloss_shifts, widths, canon_tloss_skews)
    return np.sum(vals)

