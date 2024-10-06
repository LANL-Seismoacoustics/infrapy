# infrapy.propagation.infrasound.py
#
# Infrasound propagation models for association and localization.
# General models are valid at all locations.  Propagation-based
# stochastic models are location, month, and time of day specific
#
# Author            Philip Blom (pblom@lanl.gov)

import sys
import pickle
import time
import itertools
import warnings

import numpy as np

from sklearn import mixture
from pathlib import Path

from scipy.interpolate import interp1d, interp2d, RectBivariateSpline
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from scipy.stats import norm
from scipy.signal import savgol_filter

import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

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


# ############################ #
#          Stochastic          #
#      Propagation Models      #
# ############################ #
def find_azimuth_bin(az, bin_cnt=8):
    reduced = np.degrees(np.arctan2(np.sin(np.radians(az)), np.cos(np.radians(az))))
    bins = np.arange(-180.0, 180.0, 360.0 / (bin_cnt * 2.0))

    result = np.asarray(np.digitize(reduced, bins) / 2)
    result[result >= bin_cnt] = 0
    return result.astype(int)


# Path geometry models (range-celerity and azimuth deviation)
class PathGeometryModel(object):
    az_bin_cnt = 8
    az_bin_wdth = 60.0
    bnc_max = 10
    rng_max = 1000.0

    default_az_dev_var = 4.0
    min_az_dev_var = 0.9

    tropo_strat_bnd = 1.0 / 0.31
    strat_therm_bnd = 1.0 / 0.26
    bnd_overlap = 0.075

    rcel_vrs_min = 0.05

    def __init__(self):
        self.rngs = np.array([])

        self.rcel_wts = []
        self.rcel_mns = []
        self.rcel_vrs = []

        self.az_dev_mns = []
        self.az_dev_vrs = []


    def eval_rcel_gmm(self, rng, rcel, az):
        rng_eval = np.array(rng)
        rng_eval[rng_eval > self.rng_max] = self.rng_max

        if len(np.atleast_1d(rng)) == 1:
            n_az = find_azimuth_bin(az, self.az_bin_cnt)
            fit_rcel_wts = np.array([func(rng_eval) for func in self.rcel_wts[n_az]])
            fit_rcel_mns = np.array([func(rng_eval) for func in self.rcel_mns[n_az]])
            fit_rcel_vrs = np.array([func(rng_eval) for func in self.rcel_vrs[n_az]])
            result = np.sum(fit_rcel_wts / fit_rcel_vrs * norm.pdf((rcel - fit_rcel_mns) / fit_rcel_vrs))
        else:
            mn = np.empty((len(rng_eval), 3))
            vr = np.empty((len(rng_eval), 3))
            wt = np.empty((len(rng_eval), 3))

            az_indices = find_azimuth_bin(az, self.az_bin_cnt)
            for n_az in range(self.az_bin_cnt):
                mask = [az_indices == n_az]
                if np.any(mask):
                    mn[mask] = np.array([func(rng_eval[mask]) for func in self.rcel_mns[n_az]]).T
                    vr[mask] = np.array([func(rng_eval[mask]) for func in self.rcel_vrs[n_az]]).T
                    wt[mask] = np.array([func(rng_eval[mask]) for func in self.rcel_wts[n_az]]).T

            result = np.sum(wt / vr * norm.pdf((np.array([rcel] * 3).T - mn) / vr), axis=1)
        return result


    def eval_az_dev_mn(self, rng, az):
        if len(np.atleast_1d(rng)) == 1:
            return self.az_dev_mns[find_azimuth_bin(az, self.az_bin_cnt)](min(rng, self.rng_max))
        else:
            rng_eval = np.array(rng)
            rng_eval[rng_eval > self.rng_max] = self.rng_max
            az_indices = find_azimuth_bin(az, self.az_bin_cnt)

            mn = np.empty_like(rng_eval)
            for n_az in range(self.az_bin_cnt):
                mask = az_indices==n_az
                if np.any(mask):
                    mn[mask] = self.az_dev_mns[n_az](rng_eval[mask])
            return mn

    def eval_az_dev_vr(self, rng, az):
        if len(np.atleast_1d(rng)) == 1:
            return self.az_dev_vrs[find_azimuth_bin(az, self.az_bin_cnt)](min(rng, self.rng_max))
        else:
            rng_eval = np.array(rng)
            rng_eval[rng_eval > self.rng_max] = self.rng_max
            az_indices = find_azimuth_bin(az, self.az_bin_cnt)

            vr = np.empty_like(rng_eval)
            for n_az in range(self.az_bin_cnt):
                mask = az_indices==n_az
                if np.any(mask):
                    vr[mask] = self.az_dev_vrs[n_az](rng_eval[mask])
            return vr

    def build(self, results_file, model_file, show_fits=False, file_id=None, verbose_output=False, rng_width=40.0, rng_spacing=10.0, data_format="new", abs_lim = -100.0, trn_ht_min = 2.0):
        print('-' * 75)
        print('Builing celerity and azimuth priors from file:', results_file)

        az_dirs = ['S', 'SW', 'W', 'NW', 'N', 'NE', 'E', 'SE']

        # define range bins and parameter arrays
        rng_bins = np.arange(0.0, self.rng_max, rng_spacing)
        rng_cnt = len(rng_bins)

        az_dev_mns = np.empty((self.az_bin_cnt, rng_cnt))
        az_dev_vrs = np.empty((self.az_bin_cnt, rng_cnt))

        rcel_wts = np.empty((self.az_bin_cnt, rng_cnt, 3))
        rcel_mns = np.empty((self.az_bin_cnt, rng_cnt, 3))
        rcel_vrs = np.empty((self.az_bin_cnt, rng_cnt, 3))

        arrival_cnt = np.empty(self.az_bin_cnt, dtype=int)

        # load Cartesian formatted infraGA/GeoAc predictions
        print('\t', "Loading propagation modeling results...")
        if data_format == "new":
            theta, phi, n, x, y, t, cel, z_max, incl, back_az, amp_geo, amp_atmo = np.loadtxt(results_file, unpack=True)
        else:
            theta, phi, n, x, y, t, z_max, incl, back_az, amp_geo, amp_atmo = np.loadtxt(results_file, unpack=True)

        rng = np.sqrt(x**2 + y**2)
        rcel = t / rng

        az = 90.0 - np.degrees(np.arctan2(y, x))
        az_dev = (90.0 - np.degrees(np.arctan2(-y, -x))) - back_az

        # wrap angles to +/- 180 degrees
        phi[phi > 180.0] -= 360.0
        phi[phi < -180.0] += 360.0

        az[az > 180.0] -= 360.0
        az[az < -180.0] += 360.0

        az_dev[az_dev > 180.0] -= 360.0
        az_dev[az_dev < -180.0] += 360.0

        # Cycle through azimuth bins creating fit
        for n_az in range(self.az_bin_cnt):
            print('\t', "Fitting distribution (" + az_dirs[n_az] + ")", '\t', end=' ')
            prog_bar.prep(int(np.ceil(rng_cnt / 3.0)))

            if n_az == 0:
                az_mask = np.logical_or(az > 180.0 - self.az_bin_wdth / 2.0, az < -180.0 + self.az_bin_wdth / 2.0)
            else:
                center = -180 + 360.0 / self.az_bin_cnt * n_az
                az_mask = np.logical_and(center - self.az_bin_wdth / 2.0 <= az, az <= center + self.az_bin_wdth / 2.0)

            combo_mask = np.logical_and((z_max >= trn_ht_min), az_mask)

            # create multi-bounce arrival arrays for arrivals within range and without excessive absorption
            rngs_multi_bounce = np.array([])
            rcel_multi_bounce = np.array([])
            az_dev_multi_bounce = np.array([])

            for n_bnc in range(1, self.bnc_max + 1):
                combo_mask = np.logical_and(combo_mask, n_bnc * amp_atmo > abs_lim)
                combo_mask = np.logical_and(combo_mask, float(n_bnc) * rng < self.rng_max * 1.1)

                rngs_multi_bounce = np.concatenate((rngs_multi_bounce, float(n_bnc) * rng[combo_mask]))
                rcel_multi_bounce = np.concatenate((rcel_multi_bounce, rcel[combo_mask]))
                az_dev_multi_bounce = np.concatenate((az_dev_multi_bounce, az_dev[combo_mask]))

            arrival_cnt[n_az] = len(rngs_multi_bounce)

            # display work if show_fits is enabled
            # if show_fits:
            #     f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 11))

            #     ax1.set_xlim([0.0, self.rng_max])
            #     ax1.set_ylim([0.2, 0.4])

            #     ax3.set_xlim([0.0, self.rng_max])
            #     ax3.set_ylim([-12.0, 12.0])

            #     ax1.set_xlabel('Range [km]')
            #     ax1.set_ylabel('Celerity [km/s]')

            #     ax2.set_xlabel('Celerity [km/s]')
            #     ax2.set_ylabel('Probability')

            #     ax3.set_xlabel('Range [km]')
            #     ax3.set_ylabel('Azimuth Deviation [degrees]')

            #     ax4.set_xlabel('Azimuth Deviation [degrees]')
            #     ax4.set_ylabel('Probability')

            #     plt.suptitle('Propagation-Based Priors (' + az_dirs[n_az] + ')', fontsize=18)
            #     plt.show(block=False)

            #     ax1.plot(rngs_multi_bounce, 1.0 / rcel_multi_bounce, 'k.', markersize=2.0)
            #     ax1.set_title('Celerity-range scatter')
            #     plt.pause(0.01)

            #     ax3.plot(rngs_multi_bounce, az_dev_multi_bounce, 'k.', markersize=2.0)
            #     ax3.set_title('Back-azimuth deviation')
            #     plt.pause(0.01)

            # Compute fits (1D interpolations of parameters)
            for nr in range(rng_cnt):
                rng_mask = np.logical_and(rng_bins[nr] <= rngs_multi_bounce, rngs_multi_bounce <= rng_bins[nr] + rng_width)

                tropo_mask = np.logical_and(rcel_multi_bounce < self.tropo_strat_bnd + self.bnd_overlap, rng_mask)
                strat_mask = np.logical_and(np.logical_and(self.tropo_strat_bnd - self.bnd_overlap <= rcel_multi_bounce, rcel_multi_bounce <= self.strat_therm_bnd + self.bnd_overlap), rng_mask)
                therm_mask = np.logical_and(self.strat_therm_bnd - self.bnd_overlap < rcel_multi_bounce, rng_mask)

                min_pts = 100

                if show_fits:
                    window1 = ax1.axvspan(rng_bins[nr], rng_bins[nr] + rng_width, facecolor='b', alpha=0.5)
                    window2 = ax3.axvspan(rng_bins[nr], rng_bins[nr] + rng_width, facecolor='b', alpha=0.5)
                    plt.pause(0.1)

                # set azimuth deviation fit
                if az_dev_multi_bounce[rng_mask].shape[0] > min_pts:
                    az_dev_mns[n_az][nr] = np.mean(az_dev_multi_bounce[rng_mask], dtype=np.float64)
                    az_dev_vrs[n_az][nr] = max(np.sqrt(np.var(az_dev_multi_bounce[rng_mask], dtype=np.float64)), self.min_az_dev_var)
                else:
                    if nr == 0:
                        az_dev_mns[n_az][nr] = 0.0
                        az_dev_vrs[n_az][nr] = self.default_az_dev_var
                    else:
                        az_dev_mns[n_az][nr] = az_dev_mns[n_az][nr - 1] * 0.75
                        az_dev_vrs[n_az][nr] = self.default_az_dev_var + (az_dev_vrs[n_az][nr - 1] - self.default_az_dev_var) * 0.5

                # set tropospheric contribution to reciprocal celerity fit
                if rcel_multi_bounce[tropo_mask].shape[0] > min_pts:
                    rcel_mns[n_az][nr][0] = np.mean(rcel_multi_bounce[tropo_mask], dtype=np.float64)
                    rcel_vrs[n_az][nr][0] = max(np.sqrt(np.var(rcel_multi_bounce[tropo_mask], dtype=np.float64)) * 1.5, self.rcel_vrs_min)
                else:
                    if nr == 0:
                        rcel_mns[n_az][nr][0] = canon_rcel_mns[0]
                        rcel_vrs[n_az][nr][0] = canon_rcel_vrs[0]
                    else:
                        rcel_mns[n_az][nr][0] = canon_rcel_mns[0] + (rcel_mns[n_az][nr - 1][0] - canon_rcel_mns[0]) * 0.5
                        rcel_vrs[n_az][nr][0] = canon_rcel_vrs[0] + (rcel_vrs[n_az][nr - 1][0] - canon_rcel_vrs[0]) * 0.5

                # set stratospheric contribution to reciprocal celerity fit
                if rcel_multi_bounce[strat_mask].shape[0] > min_pts:
                    rcel_mns[n_az][nr][1] = np.mean(rcel_multi_bounce[strat_mask], dtype=np.float64)
                    rcel_vrs[n_az][nr][1] = max(np.sqrt(np.var(rcel_multi_bounce[strat_mask], dtype=np.float64)) * 1.5, self.rcel_vrs_min)
                else:
                    if nr == 0:
                        rcel_mns[n_az][nr][1] = canon_rcel_mns[1]
                        rcel_vrs[n_az][nr][1] = canon_rcel_vrs[1]
                    else:
                        rcel_mns[n_az][nr][1] = canon_rcel_mns[1] + (rcel_mns[n_az][nr - 1][1] - canon_rcel_mns[1]) * 0.5
                        rcel_vrs[n_az][nr][1] = canon_rcel_vrs[1] + (rcel_vrs[n_az][nr - 1][1] - canon_rcel_vrs[1]) * 0.5

                # set thermospheric contribution to reciprocal celerity fit
                if rcel_multi_bounce[therm_mask].shape[0] > min_pts:
                    rcel_mns[n_az][nr][2] = np.mean(rcel_multi_bounce[therm_mask], dtype=np.float64)
                    rcel_vrs[n_az][nr][2] = max(np.sqrt(np.var(rcel_multi_bounce[therm_mask], dtype=np.float64)) * 1.5, self.rcel_vrs_min)
                else:
                    if nr == 0:
                        rcel_mns[n_az][nr][2] = canon_rcel_mns[2]
                        rcel_vrs[n_az][nr][2] = canon_rcel_vrs[2]
                    else:
                        rcel_mns[n_az][nr][2] = canon_rcel_mns[2] + (rcel_mns[n_az][nr - 1][2] - canon_rcel_mns[2]) * 0.5
                        rcel_vrs[n_az][nr][2] = canon_rcel_vrs[2] + (rcel_vrs[n_az][nr - 1][2] - canon_rcel_vrs[2]) * 0.5

                # set weights of reciprocal celerity fit and scale normalization by point count in range bin
                # divide by number of points at all azimuths later using: (n_rng / n_az) * (n_az / n_all) = n_rng / n_all
                rcel_wt_min = 1.0e-2
                if len(az_dev_multi_bounce[rng_mask]) > min_pts:
                    wt1 = max(float(len(rngs_multi_bounce[tropo_mask])) / len(az_dev_multi_bounce[rng_mask]), rcel_wt_min)
                    wt2 = max(float(len(rngs_multi_bounce[strat_mask])) / len(az_dev_multi_bounce[rng_mask]), rcel_wt_min)
                    wt3 = max(float(len(rngs_multi_bounce[therm_mask])) / len(az_dev_multi_bounce[rng_mask]), rcel_wt_min)

                    rcel_wts[n_az][nr][0] = wt1 / (wt1 + wt2 + wt3) * len(az_dev_multi_bounce[rng_mask])
                    rcel_wts[n_az][nr][1] = wt2 / (wt1 + wt2 + wt3) * len(az_dev_multi_bounce[rng_mask])
                    rcel_wts[n_az][nr][2] = wt3 / (wt1 + wt2 + wt3) * len(az_dev_multi_bounce[rng_mask])
                else:
                    if nr == 0:
                        rcel_wts[n_az][nr][0] = canon_rcel_wts[0] * min_pts
                        rcel_wts[n_az][nr][1] = canon_rcel_wts[1] * min_pts
                        rcel_wts[n_az][nr][2] = canon_rcel_wts[2] * min_pts
                    else:
                        rcel_wts[n_az][nr][0] = canon_rcel_wts[0] * min_pts + (rcel_wts[n_az][nr - 1][0] - canon_rcel_wts[1] * min_pts) * 0.5
                        rcel_wts[n_az][nr][1] = canon_rcel_wts[1] * min_pts + (rcel_wts[n_az][nr - 1][1] - canon_rcel_wts[1] * min_pts) * 0.5
                        rcel_wts[n_az][nr][2] = canon_rcel_wts[2] * min_pts + (rcel_wts[n_az][nr - 1][2] - canon_rcel_wts[1] * min_pts) * 0.5

                if show_fits:
                    az_dev_vals = np.linspace(-12.0, 12.0, 1000)
                    cel_vals = np.linspace(0.2, 0.4, 1000)

                    cel_dist1 = (rcel_wts[n_az][nr][0] / len(rngs_multi_bounce)) / rcel_vrs[n_az][nr][0] * norm.pdf((1.0 / cel_vals - rcel_mns[n_az][nr][0]) / rcel_vrs[n_az][nr][0])
                    cel_dist2 = (rcel_wts[n_az][nr][1] / len(rngs_multi_bounce)) / rcel_vrs[n_az][nr][1] * norm.pdf((1.0 / cel_vals - rcel_mns[n_az][nr][1]) / rcel_vrs[n_az][nr][1])
                    cel_dist3 = (rcel_wts[n_az][nr][2] / len(rngs_multi_bounce)) / rcel_vrs[n_az][nr][2] * norm.pdf((1.0 / cel_vals - rcel_mns[n_az][nr][2]) / rcel_vrs[n_az][nr][2])

                    ax2.cla()
                    ax4.cla()

                    ax2.set_xlim([0.2, 0.4])
                    ax2.plot(cel_vals, cel_dist1, linewidth=2.0, color='Green')
                    ax2.plot(cel_vals, cel_dist2, linewidth=2.0, color='Green')
                    ax2.plot(cel_vals, cel_dist3, linewidth=2.0, color='Green')
                    ax2.plot(cel_vals, cel_dist1 + cel_dist2 + cel_dist3, linewidth=4.0, color='Blue')

                    ax4.set_xlim([-12.0, 12.0])
                    ax4.plot(az_dev_vals, 1.0 / az_dev_vrs[n_az][nr] * norm.pdf((az_dev_vals - az_dev_mns[n_az][nr]) / az_dev_vrs[n_az][nr]), linewidth=4.0, color='Blue')

                    plt.pause(0.01)
                    window1.remove()
                    window2.remove()

                if nr % 3 == 0:
                    prog_bar.increment(1)
            prog_bar.close()
            plt.close('all')

        # Normalize weights by total arrivals at all azimuths
        rcel_wts /= np.sum(arrival_cnt)

        print('\n' + 'Azimuth weights are:')
        for n_az in range(self.az_bin_cnt):
            print('\t', float(arrival_cnt[n_az]) / np.sum(arrival_cnt))
            rcel_wts[n_az] /= float(np.sum(arrival_cnt))

        priors = [0] * 6
        priors[0] = rng_bins
        priors[1] = az_dev_mns
        priors[2] = az_dev_vrs
        priors[3] = rcel_mns
        priors[4] = rcel_vrs
        priors[5] = rcel_wts

        pickle.dump(priors, open(model_file, "wb"))
        print('')

    def load(self, model_file, smooth=None):
        fit_params = pickle.load(open(model_file, "rb"), encoding='latin1')
        self.az_bin_cnt = len(fit_params[1])

        self.rng_max = max(fit_params[0])

        self.az_dev_mns = [0] * self.az_bin_cnt
        self.az_dev_vrs = [0] * self.az_bin_cnt

        self.rcel_wts = [0] * self.az_bin_cnt
        self.rcel_mns = [0] * self.az_bin_cnt
        self.rcel_vrs = [0] * self.az_bin_cnt

        for n_az in range(self.az_bin_cnt):
            self.rcel_mns[n_az] = [0] * 3
            self.rcel_vrs[n_az] = [0] * 3
            self.rcel_wts[n_az] = [0] * 3

        if smooth:
            print("Loading propagation model parameters from " + model_file + " with smoothing.")
            for n_az in range(self.az_bin_cnt):
                self.az_dev_mns[n_az] = interp1d(fit_params[0], savgol_filter(fit_params[1][n_az], 9, 3), kind='cubic', bounds_error=False, fill_value=(fit_params[1][n_az][0], fit_params[1][n_az][-1]))
                self.az_dev_vrs[n_az] = interp1d(fit_params[0], savgol_filter(fit_params[2][n_az], 9, 3), kind='cubic', bounds_error=False, fill_value=(fit_params[2][n_az][0], fit_params[2][n_az][-1]))

                self.rcel_mns[n_az] = [0] * 3
                self.rcel_vrs[n_az] = [0] * 3
                self.rcel_wts[n_az] = [0] * 3

                for j in range(3):
                    self.rcel_mns[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[3][n_az][:, j], 9, 3), kind='cubic', bounds_error=False, fill_value=(fit_params[3][n_az][:, j][0], fit_params[3][n_az][:, j][-1]))
                    self.rcel_vrs[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[4][n_az][:, j], 9, 3), kind='cubic', bounds_error=False, fill_value=(fit_params[4][n_az][:, j][0], fit_params[4][n_az][:, j][-1]))
                    self.rcel_wts[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[5][n_az][:, j], 9, 3), kind='cubic', bounds_error=False, fill_value=(fit_params[5][n_az][:, j][0], fit_params[5][n_az][:, j][-1]))
        else:
            print("Loading propagation model parameters from " + model_file + " without smoothing.")
            for n_az in range(self.az_bin_cnt):
                self.az_dev_mns[n_az] = interp1d(fit_params[0], fit_params[1][n_az], kind='cubic', bounds_error=False, fill_value=(fit_params[1][n_az][0], fit_params[1][n_az][-1]))
                self.az_dev_vrs[n_az] = interp1d(fit_params[0], fit_params[2][n_az], kind='cubic', bounds_error=False, fill_value=(fit_params[2][n_az][0], fit_params[2][n_az][-1]))

                self.rcel_mns[n_az] = [0] * 3
                self.rcel_vrs[n_az] = [0] * 3
                self.rcel_wts[n_az] = [0] * 3

                for j in range(3):
                    self.rcel_mns[n_az][j] = interp1d(fit_params[0], fit_params[3][n_az][:, j], kind='cubic', bounds_error=False, fill_value=(fit_params[3][n_az][:, j][0], fit_params[3][n_az][:, j][-1]))
                    self.rcel_vrs[n_az][j] = interp1d(fit_params[0], fit_params[4][n_az][:, j], kind='cubic', bounds_error=False, fill_value=(fit_params[4][n_az][:, j][0], fit_params[4][n_az][:, j][-1]))
                    self.rcel_wts[n_az][j] = interp1d(fit_params[0], fit_params[5][n_az][:, j], kind='cubic', bounds_error=False, fill_value=(fit_params[5][n_az][:, j][0], fit_params[5][n_az][:, j][-1]))

    def display(self, file_id=None, hold_fig=None):
        resol = 100
        rngs = np.linspace(0.0, 1000.0, resol)
        bias = np.empty([resol])
        width = np.empty([resol])

        bias_color = 'Blue'
        var_color = 'LightBlue'

        compass_file = str(Path(__file__).parent / "propagation" / "compass.png")

        f1, ax = plt.subplots(3, 3, figsize=(12, 9))

        for n1, n2 in itertools.product(list(range(3)), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([-10, 10])
                ax[n1, n2].set_xlim([0, 1000])
                ax[n1, n2].set_xticks([0, 250, 500, 750, 1000])
            if n2 != 0:
                ax[n1, n2].set_yticklabels([])
            if n1 != 2:
                ax[n1, n2].set_xticklabels([])

        img = mpimg.imread(compass_file)
        ax[1, 1].axis('off')
        ax[1, 1].imshow(img)

        ax[2, 1].set_xlabel('Range [km]')
        ax[1, 0].set_ylabel('Azimuth Deviation [deg]')

        plt.suptitle('Azimuth Deviation Priors', fontsize=22)
        plt.show(block=False)

        bias = self.eval_az_dev_mn(rngs, [-45.0] * len(rngs))
        width = self.eval_az_dev_vr(rngs, [-45.0] * len(rngs))
        ax[0, 0].fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
        ax[0, 0].plot(rngs, bias, linewidth=2.0, color=bias_color)
        ax[0, 0].plot(rngs, [0] * len(rngs), 'k:')
        plt.pause(0.1)

        bias = self.eval_az_dev_mn(rngs, [0.0] * len(rngs))
        width = self.eval_az_dev_vr(rngs, [0.0] * len(rngs))
        ax[0, 1].fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
        ax[0, 1].plot(rngs, bias, linewidth=2.0, color=bias_color)
        ax[0, 1].plot(rngs, [0] * len(rngs), 'k:')
        plt.pause(0.1)

        bias = self.eval_az_dev_mn(rngs, [45.0] * len(rngs))
        width = self.eval_az_dev_vr(rngs, [45.0] * len(rngs))
        ax[0, 2].fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
        ax[0, 2].plot(rngs, bias, linewidth=2.0, color=bias_color)
        ax[0, 2].plot(rngs, [0] * len(rngs), 'k:')
        plt.pause(0.1)

        bias = self.eval_az_dev_mn(rngs, [-90.0] * len(rngs))
        width = self.eval_az_dev_vr(rngs, [-90.0] * len(rngs))
        ax[1, 0].fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
        ax[1, 0].plot(rngs, bias, linewidth=2.0, color=bias_color)
        ax[1, 0].plot(rngs, [0] * len(rngs), 'k:')
        plt.pause(0.1)

        bias = self.eval_az_dev_mn(rngs, [90.0] * len(rngs))
        width = self.eval_az_dev_vr(rngs, [90.0] * len(rngs))
        ax[1, 2].fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
        ax[1, 2].plot(rngs, bias, linewidth=2.0, color=bias_color)
        ax[1, 2].plot(rngs, [0] * len(rngs), 'k:')
        plt.pause(0.1)

        bias = self.eval_az_dev_mn(rngs, [-135.0] * len(rngs))
        width = self.eval_az_dev_vr(rngs, [-135.0] * len(rngs))
        ax[2, 0].fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
        ax[2, 0].plot(rngs, bias, linewidth=2.0, color=bias_color)
        ax[2, 0].plot(rngs, [0] * len(rngs), 'k:')
        plt.pause(0.1)

        bias = self.eval_az_dev_mn(rngs, [-180.0] * len(rngs))
        width = self.eval_az_dev_vr(rngs, [-180.0] * len(rngs))
        ax[2, 1].fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
        ax[2, 1].plot(rngs, bias, linewidth=2.0, color=bias_color)
        ax[2, 1].plot(rngs, [0] * len(rngs), 'k:')
        plt.pause(0.1)

        bias = self.eval_az_dev_mn(rngs, [135.0] * len(rngs))
        width = self.eval_az_dev_vr(rngs, [135.0] * len(rngs))
        ax[2, 2].fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
        ax[2, 2].plot(rngs, bias, linewidth=2.0, color=bias_color)
        ax[2, 2].plot(rngs, [0] * len(rngs), 'k:')
        plt.pause(0.1)

        if file_id:
            plt.savefig(file_id + "_az-dev.png", bbox_inches='tight', dpi=200)

        # Plot celerity-range statistics
        cels = np.linspace(0.2, 0.4, resol)
        rngs = np.linspace(255, 265, resol)
        R, V = np.meshgrid(rngs, cels)
        R = R.flatten()
        V = V.flatten()

        pdf = np.empty([resol, resol])

        palette = cm.nipy_spectral_r


        az_dirs = [-180.0, -135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0]
        pdf_max = 0.0
        for az in az_dirs:
            pdf_max = max(pdf_max, max(self.eval_rcel_gmm(R, 1.0 / V, [az] * len(R))))

        f2, ax = plt.subplots(3, 3, figsize=(12, 9))

        for n1, n2 in itertools.product(list(range(3)), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([0.2, 0.4])
                ax[n1, n2].set_xlim([0, 1000])
                ax[n1, n2].set_xticks([0, 250, 500, 750, 1000])
            if n2 != 0:
                ax[n1, n2].set_yticklabels([])
            if n1 != 2:
                ax[n1, n2].set_xticklabels([])

        img = mpimg.imread(compass_file)
        ax[1, 1].axis('off')
        ax[1, 1].imshow(img)

        ax[2, 1].set_xlabel('Range [km]')
        ax[1, 0].set_ylabel('Celerity [km/s]')

        plt.suptitle('Celerity-Range Priors', fontsize=22)
        plt.show(block=False)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [-45.0] * len(R))

        ax[0, 0].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        plt.pause(0.1)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [0.0] * len(R))
        ax[0, 1].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        plt.pause(0.1)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [45.0] * len(R))
        ax[0, 2].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        plt.pause(0.1)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [-90.0] * len(R))
        ax[1, 0].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        plt.pause(0.1)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [90.0] * len(R))
        ax[1, 2].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        plt.pause(0.1)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [-135.0] * len(R))
        ax[2, 0].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        plt.pause(0.1)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [-180.0] * len(R))
        ax[2, 1].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        plt.pause(0.1)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [135.0] * len(R))
        ax[2, 2].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        plt.pause(0.1)

        if file_id:
            plt.savefig(file_id + "_cel-rng.png", bbox_inches='tight', dpi=200)

        if hold_fig:
            plt.show(block=True)
        else:
            time.sleep(10)
            plt.close('all')


# Propagation-based, stochastic transmission loss model
class TLossModel(object):
    az_bin_cnt = 8
    az_bin_wdth = 60.0
    
    def __init__(self):
        self.rng_vals = [0]
        self.tloss_vals = [0]
        self.pdf_vals = [0] * self.az_bin_cnt
        self.pdf_fits = [0] * self.az_bin_cnt
    
    def build(self, results_file, model_file, show_fits=False, file_id=None, pool=None):
        print('-' * 75)
        print('Builing transmission loss priors from file:', results_file)
        
        az_dirs = ['S', 'SW', 'W', 'NW', 'N', 'NE', 'E', 'SE']
        
        # read in data, convert tloss to dB relative to 1 km, and wrap azimuths to [-180.0:180.0]
        print('\t' + "Reading in data...")
        rngs, az, tloss = np.loadtxt(results_file, unpack=True)
        
        output_rngs = np.sort(np.unique(rngs)[::5])
        
        az[az > 180.0] -= 360.0
        az[az < -180.0] += 360.0
        
        tloss = 10.0 * np.log10(tloss)
        tloss[np.isneginf(tloss)] = min(tloss[np.isfinite(tloss)])
        tloss[np.isposinf(tloss)] = max(tloss[np.isfinite(tloss)])
        
        tloss_vals = np.linspace(min(tloss) - 2.5, max(tloss) + 2.5, len(output_rngs) * 2)
        pdf_vals = np.empty((self.az_bin_cnt, len(output_rngs), len(tloss_vals)))
        
        for az_index in range(self.az_bin_cnt):
            center = -180 + 360.0 / self.az_bin_cnt * az_index
            if az_index == 0:
                az_mask = np.logical_or(az >= 180.0 - self.az_bin_wdth / 2.0, az <= -180.0 + self.az_bin_wdth / 2.0)
            else:
                az_mask = np.logical_and(center - self.az_bin_wdth / 2.0 <= az, az <= center + self.az_bin_wdth / 2.0)
        
            if show_fits:
                f, ((ax1, ax2)) = plt.subplots(2, 1, figsize=(7.5, 10))
                
                ax1.set_xlabel('Range [km]')
                ax1.set_ylabel('Transmission Loss [dB]')
                ax1.set_xlim([0.0, 1000.0])
                ax1.set_ylim([min(tloss) - 5.0, max(tloss) + 5.0])
                
                ax2.set_xlabel('Range [km]')
                ax2.set_ylabel('Transmission Loss [dB]')
                ax2.set_xlim([0.0, 1000.0])
                ax2.set_ylim([min(tloss) - 5.0, max(tloss) + 5.0])
                
                plt.suptitle("Stochastic Transmission Loss Model \n Azimuth: " + az_dirs[az_index], fontsize=18)
                plt.show(block=False)
                
                ax1.plot(rngs[az_mask][::11], tloss[az_mask][::11], 'ko', markersize=1)
                plt.pause(0.001)
            
            print('\t' + "Propagation direction (" + az_dirs[az_index] + ")..." + '\t', end=' ')
            prog_bar.prep(50)
            
            # Define tloss pdf at each range point from KDE
            for nr, rng_val in enumerate(output_rngs):
                masked_tloss = tloss[np.logical_and(az_mask, rngs == rng_val)]
                
                if np.std(masked_tloss) < 0.01:
                    pdf_vals[az_index][nr] = norm.pdf(tloss_vals, loc=np.mean(masked_tloss), scale=0.01)
                else:
                    kernel = gaussian_kde(masked_tloss)
                    pdf_vals[az_index][nr] = kernel.evaluate(tloss_vals)
            
                prog_bar.increment( int(np.floor((50.0 * (nr + 1)) / len(output_rngs)) - np.floor((50.0 * nr) / len(output_rngs))))
                
                if show_fits:
                    ax2.scatter([rng_val] * len(tloss_vals), tloss_vals, c=pdf_vals[az_index][nr], cmap=cm.nipy_spectral_r, marker='o', s=[12.5] * len(tloss_vals), alpha=0.5, edgecolor='none')
                    plt.pause(0.001)
    
            prog_bar.close()
            if show_fits:
                plt.close()
        
        priors = [0] * 3
        priors[0] = output_rngs
        priors[1] = tloss_vals
        priors[2] = pdf_vals
        
        pickle.dump(priors, open(model_file, "wb"))
        print(' ')
    
    
    def load(self, model_file):
        fit_params = pickle.load(open(model_file, "rb"), encoding='latin1')
        
        self.rng_vals = fit_params[0]
        self.tloss_vals = fit_params[1]
        
        self.az_bin_cnt = len(fit_params[2])
        self.pdf_vals = [0] * self.az_bin_cnt
        self.pdf_fits = [0] * self.az_bin_cnt
        for az_index in range(self.az_bin_cnt):
            self.pdf_vals[az_index] = fit_params[2][az_index]
            self.pdf_fits[az_index] = RectBivariateSpline(self.rng_vals, self.tloss_vals, self.pdf_vals[az_index])


    def eval(self, rng, tloss, az):
        az_index = find_azimuth_bin(az, self.az_bin_cnt)
    
        if len(np.atleast_1d(rng)) == 1:
            result = self.pdf_fits[az_index].ev(rng, tloss)
        else:
            result = np.empty(len(rng))
            for n_az in range(self.az_bin_cnt):
                mask = az_index==n_az
                if np.any(mask):
                    result[mask] = self.pdf_fits[n_az].ev(np.array(rng)[mask], np.array(tloss)[mask])
    
        return result


    def display(self, title="Transmission Loss Statistics", file_id=None, hold_fig=None):
        scale_max = 0.1
        
        compass_file = str(Path(__file__).parent / "propagation" / "compass.png")

        resol = 100
        rngs = np.linspace(0.0, 1000.0, resol)
        tloss_min, tloss_max = -60.0, 0.0
        tloss = np.linspace(tloss_min, tloss_max, resol)
        
        R, TL = np.meshgrid(rngs, tloss)
        R = R.flatten()
        TL = TL.flatten()
        
        palette = cm.nipy_spectral_r
        f1, ax = plt.subplots(3, 3, figsize=(12, 9))
        
        for n1, n2 in itertools.product(range(3), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([tloss_min, tloss_max])
                ax[n1, n2].set_xlim([0, 1000])
                ax[n1, n2].set_xticks([0, 250, 500, 750, 1000])
            if n2 != 0:
                ax[n1, n2].set_yticklabels([])
            if n1 != 2:
                ax[n1, n2].set_xticklabels([])
    
        img = mpimg.imread(compass_file)
        ax[1, 1].axis('off')
        ax[1, 1].imshow(img)
        
        ax[2, 1].set_xlabel('Range [km]')
        ax[1, 0].set_ylabel('Transmission Loss [dB]')
        
        if title:
            plt.suptitle(title, fontsize=22)
        plt.show(block=False)
        
        pdf = self.eval(R, TL, np.array([-45.0] * len(R)))
        tloss_plot = ax[0, 0].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        f1.colorbar(tloss_plot, ax=[ax[0,2], ax[1,2], ax[2,2]], label="Probability")
        plt.pause(0.1)

        pdf = self.eval(R, TL, np.array([0.0] * len(R)))
        ax[0, 1].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        plt.pause(0.1)

        pdf = self.eval(R, TL, np.array([45.0] * len(R)))
        ax[0, 2].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        plt.pause(0.1)

        pdf = self.eval(R, TL, np.array([-90.0] * len(R)))
        ax[1, 0].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        plt.pause(0.1)

        pdf = self.eval(R, TL, np.array([90.0] * len(R)))
        ax[1, 2].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        plt.pause(0.1)

        pdf = self.eval(R, TL, np.array([-135.0] * len(R)))
        ax[2, 0].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        plt.pause(0.1)

        pdf = self.eval(R, TL, np.array([-180.0] * len(R)))
        ax[2, 1].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        plt.pause(0.1)

        pdf = self.eval(R, TL, np.array([135.0] * len(R)))
        ax[2, 2].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        plt.pause(0.1)

        if file_id:
            plt.savefig(file_id + "_tloss.png", bbox_inches='tight')
        
        if hold_fig:
            plt.show(block=True)
        else:
            plt.pause(5)
            plt.close('all')
