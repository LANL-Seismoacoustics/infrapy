"""
infrapy.detection.spectral.py

Methods for analyzing a single sensor stream and
identifying detections via a spectrogram

Author            Philip Blom (pblom@lanl.gov)

"""

import numpy as np

from obspy.core import UTCDateTime

from scipy.integrate import simps
from scipy.signal import spectrogram, savgol_filter
from scipy.stats import gaussian_kde, norm, skewnorm
from scipy.optimize import curve_fit, minimize_scalar

from sklearn.cluster import DBSCAN

from ..utils import prog_bar


def calc_thresh(Sxx_vals, p_val):
    if Sxx_vals is not None:
        kernel = gaussian_kde(Sxx_vals)

        spec_spread = np.max(Sxx_vals) - np.min(Sxx_vals)
        spec_vals = np.linspace(np.min(Sxx_vals) - 0.25 * spec_spread, np.max(Sxx_vals) + 0.25 * spec_spread, 100)

        mean0 = simps(spec_vals * kernel(spec_vals), spec_vals)
        stdev0 = np.sqrt(simps((spec_vals - mean0)**2 * kernel(spec_vals), spec_vals))
        thresh0 = norm.ppf(1.0 - p_val, loc=mean0, scale=stdev0)

        mask = np.logical_and(mean0 - 2.0 * stdev0 < spec_vals, spec_vals < mean0 + 2.0 * stdev0)

        def temp(x, sk, A0, x0, sig0):
            return A0 * skewnorm.pdf(x, sk, loc=x0, scale=sig0)

        popt, _ = curve_fit(temp, spec_vals[mask], kernel(spec_vals[mask]), p0=(0.0, 1.0, mean0, stdev0))
        thresh_fit = skewnorm.ppf(1.0 - p_val, popt[0], loc=popt[2], scale=popt[3])
        thresh = min(thresh0, thresh_fit)

        def temp2(x):
            return -skewnorm.pdf(x, popt[0], loc=popt[2], scale=popt[3])
        peak = minimize_scalar(temp2, bracket=(popt[2] - 2.0 * popt[3], popt[2] + 2.0 * popt[3])).x

        return thresh, peak
    else:
        return 0.0, 0.0


def calc_thresh_wrapper(args):
    return calc_thresh(*args)


def det2dict(f, t, Sxx_log, det_pnts, trace, peaks_history, thresh_history, times_history):

        t0 = trace.stats.starttime

        t_mean = UTCDateTime(t0) + np.mean(det_pnts[:, 0])
        t1 = UTCDateTime(t0) + min(det_pnts[:, 0])
        t2 = UTCDateTime(t0) + max(det_pnts[:, 0])

        t_mid = UTCDateTime(t0) + np.mean(det_pnts[:, 0])
        tm_index = np.argmin([abs(tn - t_mid) for tn in times_history])
        bg_freqs = f[peaks_history[tm_index] != 0]
        bg_peaks = peaks_history[tm_index][peaks_history[tm_index] != 0]
        bg_thresh = thresh_history[tm_index][peaks_history[tm_index] != 0]

        det_info = dict()
        det_info['Time (UTC)'] = str(t_mean)
        det_info['Start'] = t1 - t_mean 
        det_info['End'] = t2 - t_mean
        det_info['Freq Range'] = [np.round(min(det_pnts[:, 1]), 2),
                                  np.round(max(det_pnts[:, 1]), 2)]

        try:
            det_info['Latitude'] = float(trace.stats.sac['stla'])
            det_info['Longitude'] = float(trace.stats.sac['stlo'])
        except:
            print("Lat/Lon info not in trace header, omitting from detection file.")

        det_info['Network'] = trace.stats.network
        det_info['Station'] = trace.stats.station
        det_info['Channel'] = trace.stats.channel

        det_info['Sxx_points'] = det_pnts

        SXX_det_mask = np.logical_and(min(det_pnts[:, 0]) < t, t < max(det_pnts[:, 0]))
        det_info['Sxx_det_mean'] = [f, np.mean(Sxx_log[:, SXX_det_mask], axis=1)]
        det_info['Sxx_det_max'] = [f, np.max(Sxx_log[:, SXX_det_mask], axis=1)]

        det_info['Background Peaks'] = [bg_freqs, bg_peaks]
        det_info['Background Threshold'] = [bg_freqs, bg_thresh]

        return det_info


def run_sd(trace, freq_band, spec_overlap, p_val, adaptive_window_length, adaptive_window_step, smoothing_factor, 
            clustering_freq_scaling, clustering_eps, clustering_min_samples, pl):
    """Run the spectral detection (sd) methods

        trace: obspy.core.Trace
            Obspy trace containing single channel data
        freq_band: 1darray
            Iterable with minimum and maximum frequencies for analysis
        spec_overlap: float
            Overlap factor for computing spectrogram (noverlap = nperseg * spec_overlap)
        p_val: float
            P-value for spectrogram background analysis
        adaptive_window_length: float
            Adaptive window length in seconds
        adaptive_window_step: float
            Adaptive window step in seconds (np.unique used to remove duplicated above-threshold points)
        smoothing_factor: float
            Smoothing factor (not currently used)
        clustering_freq_scaling: float
            Mapping from frequency to psuedo-time (\tau = S*log10(f))
        clustering_eps: float
            Linkage distance for DBSCAN (eps)
        clustering_min_sample: int
            Count of required members in a cluster in DBSCAN
        pl: multiprocessing.Pool
            Multiprocessing pool for simulatenous analysis of windows


        Returns:
        ----------
        dets: iterable of dicts
            List of dictionaries containing detection info
        """


    print('\n' + "Running spectral detection (sd) analysis...")

    # Compute spectrogram from the trace
    dt = trace.stats.delta
    nperseg = int((4.0 / freq_band[0]) / dt)

    f, t, Sxx = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap = int(nperseg * spec_overlap))
    freq_band_mask = np.logical_and(freq_band[0] < f, f < freq_band[1])
    Sxx_log = 10.0 * np.log10(Sxx)

    if freq_band[1] > f[-1]:
        print("Warning!  Maximum frequency is above Nyquist (" + str(f[-1]) + ")")

    # Scan through adaptive windows to identify above-background spectrogram points
    thresh_history, peaks_history, times_history = [], [], []
    spec_dets = []

    prog_bar_len, win_cnt = 50, np.ceil((t[-1] - t[0]) / adaptive_window_step)
    print('\t' + "Progress: ", end = '')
    prog_bar.prep(prog_bar_len)

    for win_n, window_start in enumerate(np.arange(t[0], t[-1], adaptive_window_step)):
        window_mask = np.logical_and(window_start <= t, t <= window_start + adaptive_window_length)
        Sxx_window = Sxx_log[:, window_mask]
        t_window = t[window_mask]

        if pl is not None:
            args = [[Sxx_window[fn], p_val] if freq_band[0] < f[fn] and f[fn] < freq_band[1] else [None, False] for fn in range(len(f))]
            temp = pl.map(calc_thresh_wrapper, args)
        else:
            temp = np.array([calc_thresh(Sxx_window[fn], p_val) if (freq_band[0] < f[fn] and f[fn] < freq_band[1]) else (0.0, 0.0) for fn in range(len(f))])

        threshold = np.array(temp)[:, 0]
        peaks = np.array(temp)[:, 1]

        if smoothing_factor is not None:
            if smoothing_factor > 2:
                threshold[freq_band_mask] = savgol_filter(threshold[freq_band_mask], smoothing_factor * 2, smoothing_factor)
                peaks[freq_band_mask] = savgol_filter(peaks[freq_band_mask], smoothing_factor * 2, smoothing_factor)

        thresh_history = thresh_history + [threshold]
        peaks_history = peaks_history + [peaks]
        times_history = times_history + [UTCDateTime(trace.stats.starttime) + (window_start + adaptive_window_length / 2.0)]

        for fn in range(len(f)):
            if freq_band[0] < f[fn] and f[fn] < freq_band[1]:
                spec_dets = spec_dets + [[t_window[tk], f[fn], Sxx_window[fn][tk]] for tk in range(len(t_window)) if Sxx_window[fn][tk] >= threshold[fn]]                

        prog_bar.increment(prog_bar.set_step(win_n, win_cnt, prog_bar_len))

    prog_bar.close()

    # Remove duplicate above-threshold points and convert histories to numpy arrays
    spec_dets = np.unique(np.array(spec_dets), axis=0)
    thresh_history = np.array(thresh_history)
    peaks_history = np.array(peaks_history)

    # Cluster into detections
    spec_dets_logf = np.stack((spec_dets[:, 0], clustering_freq_scaling * np.log10(spec_dets[:, 1]))).T
    clustering = DBSCAN(eps=clustering_eps, min_samples=clustering_min_samples).fit(spec_dets_logf)
    print("Identified " + str(max(clustering.labels_) + 1) + " detections." + '\n')

    det_list = [det2dict(f, t, Sxx_log, spec_dets[clustering.labels_ == k], trace, peaks_history, thresh_history, times_history) for k in range(max(clustering.labels_) + 1)]

    return det_list 
