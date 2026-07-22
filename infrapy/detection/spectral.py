"""
infrapy.detection.spectral.py

Methods for analyzing a single sensor stream and
identifying detections via a spectrogram

Author            Philip Blom (pblom@lanl.gov)

"""

import pywt
import numpy as np

from obspy.core import UTCDateTime

from scipy.integrate import simpson
from scipy.signal import spectrogram, stft
from scipy.stats import gaussian_kde, norm, skewnorm
from scipy.optimize import curve_fit, minimize_scalar

from sklearn.cluster import DBSCAN

from ..utils import prog_bar


def calc_thresh(Sxx_vals, p_val):
    if Sxx_vals is not None:
        kernel = gaussian_kde(Sxx_vals)
 
        spec_spread = np.max(Sxx_vals) - np.min(Sxx_vals)
        spec_vals = np.linspace(np.min(Sxx_vals) - 0.25 * spec_spread, np.max(Sxx_vals) + 0.25 * spec_spread, 100)
 
        mean0 = simpson(spec_vals * kernel(spec_vals), spec_vals)
        stdev0 = np.sqrt(simpson((spec_vals - mean0)**2 * kernel(spec_vals), spec_vals))
        thresh0 = norm.ppf(1.0 - p_val, loc=mean0, scale=stdev0)
 
        mask = np.logical_and(mean0 - 2.0 * stdev0 < spec_vals, spec_vals < mean0 + 2.0 * stdev0)
 
        try:
            # Try fitting with a skew normal
            def skew_fit(x, sk, A0, x0, sig0):
                return A0 * skewnorm.pdf(x, sk, loc=x0, scale=sig0)
 
            popt, _ = curve_fit(skew_fit, spec_vals[mask], kernel(spec_vals[mask]), p0=(0.0, 1.0, mean0, stdev0))
            thresh_fit = skewnorm.ppf(1.0 - p_val, popt[0], loc=popt[2], scale=popt[3])
            thresh = min(thresh0, thresh_fit)
 
            def temp2(x):
                return -skewnorm.pdf(x, popt[0], loc=popt[2], scale=popt[3])
            peak = minimize_scalar(temp2, bracket=(popt[2] - 2.0 * popt[3], popt[2] + 2.0 * popt[3])).x
        except:
            # Fit using a standard normal distribution if the skew fit fails to converge
            def norm_fit(x, A0, x0, sig0):
                return A0 * norm.pdf(x, loc=x0, scale=sig0)
           
            popt, _ = curve_fit(norm_fit, spec_vals[mask], kernel(spec_vals[mask]), p0=(1.0, mean0, stdev0))
            thresh_fit = norm.ppf(1.0 - p_val, loc=popt[1], scale=popt[2])
 
            thresh = min(thresh0, thresh_fit)
            peak = popt[1]
 
        return thresh, peak
    else:
        return 0.0, 0.0

def calc_thresh_wrapper(args):
    return calc_thresh(*args)


def det2dict(f, t, Sxx_log, det_pnts, trace, peaks_history, thresh_history, times_history):
        t0 = UTCDateTime(trace.stats.starttime)

        dt_det = np.mean(det_pnts[:, 0])
        tm_det = t0 + dt_det

        dt_start = min(det_pnts[:, 0])
        dt_end = max(det_pnts[:, 0])
        
        if dt_end == dt_start:
            dt_start = dt_start - 1
            dt_end = dt_end + 1

        det_buffer = (dt_end - dt_start) * 0.15
        det_buffer = min(det_buffer, 60.0)
        det_buffer = max(det_buffer, 15.0)

        dt1 = dt_start - det_buffer
        dt2 = dt_end + det_buffer

        # write detection info into dictionary
        det_info = dict()
        det_info['peak f-stat time'] = str(tm_det)

        # band passed wvform
        tr_bandpass = trace.copy()
        tr_bandpass.detrend()
        tr_bandpass.filter('bandpass', freqmin=min(det_pnts[:, 1]), freqmax=max(det_pnts[:, 1]))
        tr_bandpass.trim(t0 + dt1, t0 + dt2)
        
        det_info["waveform"] = [tr_bandpass.times() - (dt_det - dt1), tr_bandpass.data]

        # spectrogram
        det_pnts[:, 0] = det_pnts[:, 0] - dt_det
        det_info['spec pnts'] = det_pnts

        SXX_det_mask = np.logical_and(dt1 < t, t < dt2)
        det_info['spectrogram'] = [f, t[SXX_det_mask] - dt_det, Sxx_log[:, SXX_det_mask]]

        # spectral curves
        SXX_det_mask = np.logical_and(dt_start < t, t < dt_end)
        spec_mean = np.mean(Sxx_log[:, SXX_det_mask], axis=1)
        spec_max = np.max(Sxx_log[:, SXX_det_mask], axis=1)

        det_info['spec'] = [f, spec_mean, spec_max]
        tm_index = np.argmin([abs(tn - dt_det) for tn in times_history])

        bg_freqs = f[peaks_history[tm_index] != 0]
        bg_peaks = peaks_history[tm_index][peaks_history[tm_index] != 0]
        bg_thresh = thresh_history[tm_index][peaks_history[tm_index] != 0]
        det_info['bg spec'] = [bg_freqs, bg_peaks, bg_thresh]

        return det_info


def run_sd(f, t, Sxx_log, freq_band, p_val, adaptive_window_length, adaptive_window_step,
            clustering_freq_scaling, clustering_eps, clustering_min_samples, clustering_window_len,
            pl, t_skip, verbose=False):

    """Run the spectral detection (sd) methods
 
        NEED TO UPDATE THIS NOW THAT WE'VE SEPARATED FUNCTIONS
 
        trace: obspy.core.Trace
            Obspy trace containing single channel data
        spec_option: str
            Spectrogram method ('spectrogra', 'stft', or 'cwt')
        morelet_omega0: float
            Frequency scalar for Morlet waveleth used in 'cwt' option
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
        clustering_freq_scaling: float
            Mapping from frequency to psuedo-time (\tau = S*log10(f))
        clustering_eps: float
            Linkage distance for DBSCAN (eps)
        clustering_min_sample: int
            Count of required members in a cluster in DBSCAN
        clustering_window_len : float
            Length of window used in DBSCAN (avoids memory issues)
        pl: multiprocessing.Pool
            Multiprocessing pool for simultaneous analysis of windows
 
        Returns:
        ----------
        dets: iterable of dicts
            List of dictionaries containing detection info
    """
 
    if verbose:
        print('\n' + "Running spectral detection analysis...")
 
    if freq_band[1] > f[-1]:
        print("Warning!  Maximum frequency is above Nyquist (" + str(f[-1]) + ")")
    freq_band_mask = np.logical_and(freq_band[0] < f, f < freq_band[1])
 
    # Scan through adaptive windows to identify above-background spectrogram points
    thresh_history, peaks_history, times_history = [], [], []
    spec_dets = []
 
    prog_bar_len, win_cnt = 50, np.ceil((t[-1] - t[0]) / adaptive_window_step)
    if verbose:
        print("  Analyzing spectrogram... ", end = '\t')
        prog_bar.prep(prog_bar_len)
 
    for win_n, window_start in enumerate(np.arange(t[0], t[-1], adaptive_window_step)):
        window_mask = np.logical_and(window_start <= t, t <= window_start + adaptive_window_length)
        Sxx_window = Sxx_log[:, window_mask]
        t_window = t[window_mask]
 
        if pl is not None:
            args = [[Sxx_window[:, ::t_skip][fn], p_val] if freq_band[0] < f[fn] and f[fn] < freq_band[1] else [None, False] for fn in range(len(f))]
            temp = pl.map(calc_thresh_wrapper, args)
        else:
            temp = np.array([calc_thresh(Sxx_window[:, ::t_skip][fn], p_val) if freq_band_mask[fn] else (0.0, 0.0) for fn in range(len(f))])
 
        threshold = np.array(temp)[:, 0]
        peaks = np.array(temp)[:, 1]
 
        thresh_history = thresh_history + [threshold]
        peaks_history = peaks_history + [peaks]
        times_history = times_history + [(window_start + adaptive_window_length / 2.0)]
 
        _, thresh_grid = np.meshgrid(t_window, threshold)
        spec_dets = spec_dets + [[t_window[k], f[fn], Sxx_window[fn][k]] for fn, k in np.argwhere(Sxx_window > thresh_grid) if freq_band_mask[fn]]      
        if verbose:
            prog_bar.increment(prog_bar.set_step(win_n, win_cnt, prog_bar_len))
 
    if verbose:
        prog_bar.close()
 
    # Remove duplicate above-threshold points and convert histories to numpy arrays
    spec_dets = np.unique(np.array(spec_dets), axis=0)
    thresh_history = np.array(thresh_history)
    peaks_history = np.array(peaks_history)
 
    history_info = [peaks_history, thresh_history, times_history]
 
    # Cluster into detections
    if verbose:
        print("  Clustering into detections...", end = '\t')
        prog_bar.prep(prog_bar_len)

    win_cnt = np.ceil(np.max(spec_dets[:, 0]) / clustering_window_len)

    cluster_results = []
    for dt in np.arange(t[0], t[-1], clustering_window_len):
        t1 = dt
        t2 = dt + clustering_window_len * 1.2
        
        tm_mask = np.logical_and(t1 <= spec_dets[:, 0], spec_dets[:, 0] <= t2)
        if spec_dets[tm_mask].shape[0] > clustering_min_samples:
            # spec_dets_logf = np.stack((spec_dets[tm_mask, 0], clustering_freq_scaling * np.log10(spec_dets[tm_mask, 1]))).T
            spec_dets_logf = np.stack((spec_dets[tm_mask, 0], clustering_freq_scaling * np.log10(spec_dets[tm_mask, 1]) + 0.0 * spec_dets[tm_mask, 1])).T

            clustering = DBSCAN(eps=clustering_eps, min_samples=clustering_min_samples).fit(spec_dets_logf)
            cluster_results += [spec_dets[tm_mask][clustering.labels_ == k] for k in range(0, max(clustering.labels_) + 1)]
        
            if verbose:
                prog_bar.increment(prog_bar.set_step(win_n, win_cnt, prog_bar_len))

    if verbose:
        prog_bar.close()

    cluster_cnt = len(cluster_results)
    for n1 in range(cluster_cnt):
        for n2 in range(n1 + 1, cluster_cnt):
            if len(list(cluster_results[n1])) > 0 and len(list(cluster_results[n2])) > 0:
                set1 = set([tuple(x) for x in cluster_results[n1]])
                set2 = set([tuple(x) for x in cluster_results[n2]])
                overlap = np.array([x for x in set1 & set2]).shape[0]

                if overlap > 5:
                    cluster_results[n1] = np.array(list(set1.union(set2)))
                    cluster_results[n2] = np.array([])

    cluster_results = [cl for cl in cluster_results if len(list(cl)) > 0]

    if verbose:
        print('\nIdentified ' + str(len(cluster_results)) + " detections." + '\n')
 
    return spec_dets, cluster_results, history_info
 

def cli_sd(trace, spec_option, morlet_omega0, freq_band, spec_overlap, p_val, adaptive_window_length, adaptive_window_step, clustering_freq_scaling, clustering_eps, clustering_min_samples, cluster_window_len, pl):

    # Compute spectrogram from the trace
    dt = trace.stats.delta
    nperseg = int((8.0 / freq_band[0]) / dt)
    t_skip = 1
    cwt_t_skip = 4
 
    if spec_option == "spectrogram":
        f, t, Sxx = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * spec_overlap))
        Sxx_log = 10.0 * np.log10(Sxx)
    elif spec_option == "stft":
        f, t, Sxx = stft(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * spec_overlap))
        Sxx_log = 10.0 * np.log10(abs(Sxx))
    elif spec_option == "cwt":
        f, _, _ = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap=int(nperseg * spec_overlap))

        wavelet_name = 'cmor1.0-' + str(morlet_omega0 / (2 * np.pi))
        scales = pywt.frequency2scale(wavelet_name, f[1:] * dt)      
        coefficients, f = pywt.cwt(trace.data, scales, wavelet_name, sampling_period=dt)

        t = trace.times()[::cwt_t_skip]
        t_skip = max(1, int(nperseg * (1.0 - spec_overlap) / 4.0))

        Sxx_log = 10.0 * np.log10(abs(coefficients[:, ::cwt_t_skip]))
    else:
        print("Error: unrecognized spectrogram option: " + spec_option + ".")
        return []
    
    if (t[1] - t[0]) > clustering_eps:
        print("Warning!!  clustering_eps less than spectrogram temporal resolution, adjusting to " + str(clustering_eps))
        clustering_eps = (t[1] - t[0]) * 1.1
   
    _, cluster_results, history = run_sd(f, t, Sxx_log, freq_band, p_val, adaptive_window_length, adaptive_window_step, clustering_freq_scaling, clustering_eps, clustering_min_samples, cluster_window_len, pl, t_skip, verbose=True)
 
    times_history = [UTCDateTime(trace.stats.starttime) + tn for tn in history[2]]
    det_list = [det2dict(f, t, Sxx_log, cluster_results[k], trace, history[0], history[1], times_history) for k in range(len(cluster_results))]
 
    return det_list, [f, t, Sxx_log], history


##########################
## Dictionary Extaction ##
##########################
def spec_det_dict(trace, sd_params, pl):
        return cli_sd(trace,
                      sd_params["spectral_option"],
                      sd_params["morlet_omega0"],
                      [sd_params["freq_min"], sd_params["freq_max"]],
                      0.9,  # hard coded 90% overlap of spectrogram
                      sd_params["p_value"],
                      sd_params["window_len"],
                      sd_params["window_step"],
                      sd_params["freq_tm_factor"],
                      sd_params["cluster_eps"],
                      sd_params["cluster_min_samples"],
                      sd_params["cluster_window_len"], pl)
        