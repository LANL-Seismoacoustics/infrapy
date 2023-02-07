#!/usr/bin/env python

# test.py
#
# Tutorial for ObsPy stream and trace manipulations
#
# Philip Blom (pblom@lanl.gov)
# Created 10/27/2022
# Last Modified 10/27/2022

import os 
import numpy as np
import json 

import pathos.multiprocessing as mp

import matplotlib.pyplot as plt
from matplotlib import cm

from obspy.core import read, UTCDateTime

from scipy.integrate import simps
from scipy.signal import spectrogram, savgol_filter
from scipy.stats import gaussian_kde, norm, skewnorm
from scipy.optimize import curve_fit, minimize_scalar

from sklearn.cluster import DBSCAN

from infrapy.detection import spectral as spec_det
from infrapy.utils import data_io as infrapy_data_io
from infrapy.detection import visualization as det_vis



if __name__ == '__main__':
    # ######################### #
    #     Define Parameters     #
    # ######################### #
    data_file = "data/YJ.BRP1..EDF.SAC"

    freq_min, freq_max = 0.2, 30.0
    spec_overlap = 0.75

    p_val = 0.01
    threshold_window = 900.0
    threshold_overlap = 0.5
    smoothing_factor = 4

    run_clustering = False

    clustering_freq_dist = 35.0
    clustering_eps = 10.0
    clustering_min_samples = 40

    freq_plot = None
    file_out = "test.json"
    fig_out = "test.png"

    # pl = None
    pl = mp.Pool(14)


    # ######################### #
    #       Read data and       #
    #    compute spectrogram    #
    # ######################### #

    tr = read(data_file)[0]
    dt = tr.stats.delta
    nperseg = int((4.0 / freq_min) / dt) 

    f, t, Sxx = spectrogram(tr.data, 1.0 / dt, nperseg=nperseg, noverlap = int(nperseg * spec_overlap))
    freq_band_mask = np.logical_and(freq_min < f, f < freq_max)
    Sxx_log = 10.0 * np.log10(Sxx)

    # ########################## #
    #     Compute thresholds     #
    # ########################## #
    thresh_history = []
    peaks_history = []
    times_history = []

    spec_dets = []
    for window_start in np.arange(t[0], t[-1], threshold_window * min(1.0, max(0.25, 1.0 - threshold_overlap))):
        print("Running analysis at " + str(window_start + threshold_window / 2.0) + "...")
        window_mask = np.logical_and(window_start <= t, t <= window_start + threshold_window)
        Sxx_window = Sxx_log[:, window_mask]
        t_window = t[window_mask]

        if pl is not None:
            args = [[Sxx_window[fn], p_val] if freq_min < f[fn] and f[fn] < freq_max else [None, False] for fn in range(len(f))]
            temp = pl.map(spec_det.calc_thresh_wrapper, args)
        else:
            temp = np.array([spec_det.calc_thresh(Sxx_window[fn], p_val) if (freq_min < f[fn] and f[fn] < freq_max) else 0.0 for fn in range(len(f))])

        threshold = np.array(temp)[:, 0]
        peaks = np.array(temp)[:, 1]

        if smoothing_factor > 2:
            threshold[freq_band_mask] = savgol_filter(threshold[freq_band_mask], smoothing_factor * 2, smoothing_factor)
            peaks[freq_band_mask] = savgol_filter(peaks[freq_band_mask], smoothing_factor * 2, smoothing_factor)

        thresh_history = thresh_history + [threshold]
        peaks_history = peaks_history + [peaks]
        times_history = times_history + [UTCDateTime(tr.stats.starttime) + (window_start + threshold_window / 2.0)]

        for fn in range(len(f)):
            if freq_min < f[fn] and f[fn] < freq_max:
                spec_dets = spec_dets + [[t_window[tk], f[fn], Sxx_window[fn][tk]] for tk in range(len(t_window)) if Sxx_window[fn][tk] >= threshold[fn]]                

        if window_start + threshold_window > t[-1]:
            break

    spec_dets = np.unique(np.array(spec_dets), axis=0)
    thresh_history = np.array(thresh_history)
    peaks_history = np.array(peaks_history)

    fig, a = plt.subplots(3, sharex=True)
    a[0].plot(tr.times(), tr.data, '-k')

    f_grid, t_grid = np.meshgrid(f, t)
    a[1].scatter(t_grid.flatten(), f_grid.flatten(), c=Sxx_log.T.flatten(), marker="s", s=2.5, cmap=cm.jet)
    a[1].set_yscale('log')
    a[2].set_yscale('log')

    a[2].set_xlabel("Time [s]")
    a[2].set_ylabel("Frequency [Hz]")
    a[1].set_ylabel("Frequency [Hz]")
    a[0].set_ylabel("Amplitude")

    a[1].axhline(freq_min, color='0.5')
    a[1].axhline(freq_max, color='0.5')

    # ######################### #
    #  Cluster into detections  #
    # ######################### #
    spec_dets_logf = np.stack((spec_dets[:, 0], clustering_freq_dist * np.log10(spec_dets[:, 1]))).T
    clustering = DBSCAN(eps=clustering_eps, min_samples=clustering_min_samples).fit(spec_dets_logf)
    print('\n' + "Identified " + str(max(clustering.labels_) + 1) + " detections.")

    det_list = []
    a[2].plot(spec_dets[:, 0][clustering.labels_ == -1], spec_dets[:, 1][clustering.labels_ == -1], '.', color='0.5', markersize=1.5)
    for k in range(max(clustering.labels_) + 1):
        det_list = det_list + [spec_det.det2dict(f, t, Sxx_log, spec_dets[clustering.labels_ == k], tr, peaks_history, thresh_history, times_history)]
        a[2].plot(spec_dets[:, 0][ clustering.labels_ == k], spec_dets[:, 1][ clustering.labels_ == k], '.', markersize=1.5)            

    a[2].set_yscale('log')
    a[2].set_ylim(freq_min, freq_max)

    # Save figure and detection info
    plt.savefig(fig_out, dpi=300)
    with open(file_out, 'w') as of:
        json.dump(det_list, of, indent=4, cls=infrapy_data_io.Infrapy_Encoder)

    # Load and plot one of the detections...
    det_list = json.load(open(file_out))
    det_vis.plot_sd_single(tr, det_list[3], show_fig=True)

    if pl is not None:
        pl.close()
    