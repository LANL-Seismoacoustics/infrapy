# visualization.py
#
# Visualization methods for beamforming (fk) and detection (fd) results
#
# Philip Blom (pblom@lanl.gov)


import numpy as np

from obspy import UTCDateTime

from scipy.signal import spectrogram

import matplotlib.pyplot as plt 
from matplotlib import cm

from . import beamforming_new


def plot_fk1(stream, latlon, times, peaks, detections=None, title=None, output_path=None, show_fig=True, det_thresh=None):
    '''
    Visualize beamforming (fk) results with waveform data included

    '''   
    x, t, t0, _ = beamforming_new.stream_to_array_data(stream, latlon)

    f, a = plt.subplots(4, figsize=(10, 6), sharex=True)
    a[3].set_xlabel("Time")
    a[0].set_ylabel("F-stat")    
    a[1].set_ylabel("Tr. Vel. [m/s]")
    a[2].set_ylabel("Back Az. [deg.]")
    a[3].set_ylabel("Pr. [Pa]")

    a[3].plot(np.array([t0 + np.timedelta64(int(tn * 1000.0), 'ms') for tn in t]), x[0,:], '-k')
    a[2].plot(times, peaks[:, 0], '.k', markersize=4)
    a[1].plot(times, peaks[:, 1], '.k', markersize=4)
    a[0].plot(times, peaks[:, 2], '.k', markersize=4)

    if detections:
        for det in detections:
            t1 = det.peakF_UTCtime + np.timedelta64(int(det.start * 1000.0), 'ms')
            t2 = det.peakF_UTCtime + np.timedelta64(int(det.end * 1000.0), 'ms')
            for n in range(4):
                a[n].axvspan(t1, t2, color="steelblue")

    if det_thresh is not None:
        a[0].plot(det_thresh[0], det_thresh[1], '--k', linewidth=0.5)

    if title:
        a[0].set_title(title)

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()

def plot_fk2(times, peaks, detections=None, title=None, output_path=None, show_fig=True):
    '''
    Visualize beamforming (fk) results without waveform data
    '''

    f, a = plt.subplots(3, figsize=(10, 6), sharex=True)
    a[2].set_xlabel("Time")
    a[0].set_ylabel("F-stat")    
    a[1].set_ylabel("Tr. Vel. [m/s]")
    a[2].set_ylabel("Back Az. [deg.]")

    a[2].plot(times, peaks[:, 0], '.k', markersize=4)
    a[1].plot(times, peaks[:, 1], '.k', markersize=4)
    a[0].plot(times, peaks[:, 2], '.k', markersize=4)

    if detections:
        for det in detections:
            t1 = det.peakF_UTCtime + np.timedelta64(int(det.start * 1000.0), 'ms')
            t2 = det.peakF_UTCtime + np.timedelta64(int(det.end * 1000.0), 'ms')
            for n in range(3):
                a[n].axvspan(t1, t2, color="steelblue")

    if title:
        a[0].set_title(title)
       
    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_sd(trace, det_list, freq_band, title=None, output_path=None, show_fig=False):
    '''
    Visualize multiple spectral detection (sd) results

    '''  
    trace.filter('bandpass', freqmin=min([det['Freq Range'][0] for det in det_list]),
                             freqmax=max([det['Freq Range'][1] for det in det_list]))

    dt = trace.stats.delta
    nperseg = int((4.0 / freq_band[0]) / dt) 
    f, t, Sxx = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap = int(nperseg * 0.75))
    Sxx_log = 10.0 * np.log10(Sxx)

    fig, a = plt.subplots(3, sharex=True, figsize=(9, 5))
    a[0].plot(trace.times(), trace.data, '-k')

    f_grid, t_grid = np.meshgrid(f, t)
    a[1].scatter(t_grid.flatten(), f_grid.flatten(), c=Sxx_log.T.flatten(), marker="s", s=2.5, cmap=cm.jet)

    a[2].set_xlabel("Time [s]")
    a[2].set_ylabel("Frequency [Hz]")
    a[1].set_ylabel("Frequency [Hz]")
    a[0].set_ylabel("Amplitude")

    a[1].axhline(freq_band[0], color='0.5')
    a[1].axhline(freq_band[1], color='0.5')

    a[2].set_ylim(freq_band[0], freq_band[1])
    for det in det_list:
        Sxx_pnts = np.array(det['Sxx_points'])
        a[2].plot(Sxx_pnts[:, 0], Sxx_pnts[:, 1], '.', markersize=1.5)

    a[1].set_yscale('log')
    a[2].set_yscale('log')

    if title:
        a[0].set_title(title)
       
    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_sd_single(trace, det_info, freq_band, title=None, output_path=None, show_fig=False):
    '''
    Visualize a single spectral detection (sd) result

    '''  
    t_shift = UTCDateTime(det_info['Time (UTC)']) - UTCDateTime(trace.stats.starttime)

    dt = trace.stats.delta
    nperseg = int((4.0 / freq_band[0]) / dt) 
    f, t, Sxx = spectrogram(trace.data, 1.0 / dt, nperseg=nperseg, noverlap = int(nperseg * 0.95))
    f_grid, t_grid = np.meshgrid(f, t)
    Sxx_log = 10.0 * np.log10(Sxx)
    Sxx_pnts = np.array(det_info['Sxx_points'])

    fig = plt.figure(figsize=(10, 5), layout="constrained")
    spec = fig.add_gridspec(3, 5)

    ax1 = fig.add_subplot(spec[1, :3])
    ax1.scatter(t_grid.flatten() - t_shift, f_grid.flatten(), c=Sxx_log.T.flatten(), marker="s", s=2.5, cmap=cm.jet)
    ax1.set_xlim([det_info['Start'] - 150.0, det_info['End'] + 150.0])
    ax1.axhline(freq_band[0], color='0.5')
    ax1.axhline(freq_band[1], color='0.5')
    ax1.set_yscale('log')
    ax1.set_ylabel("Frequency [Hz]")

    ax0 = fig.add_subplot(spec[2, :3], sharex=ax1, sharey=ax1)
    ax0.plot(Sxx_pnts[:, 0] - t_shift, Sxx_pnts[:, 1], '.', markersize=1.5, color='black')
    ax0.set_yscale('log')
    ax0.set_xlabel("Time (rel. " + det_info['Time (UTC)'] + ") [s]")
    ax0.set_ylabel("Frequency [Hz]")
    
    trace.filter('bandpass', freqmin=det_info['Freq Range'][0], freqmax=det_info['Freq Range'][1])

    ax2 = fig.add_subplot(spec[0, :3], sharex=ax1)
    ax2.plot(trace.times() - t_shift, trace.data, '-k')
    ax2.axvspan(det_info['Start'], det_info['End'], color='green', alpha=0.5)
    ax2.set_xlabel("")
    ax2.set_ylabel("Press. [Pa]")

    ax3 = fig.add_subplot(spec[1:, 3:])
    ax3.semilogx(det_info['Background Peaks'][0], det_info['Background Peaks'][1], '-k', linewidth=1.0, label="Background (Mean)")
    ax3.semilogx(det_info['Background Threshold'][0], det_info['Background Threshold'][1], '--k', linewidth=1.0, label="Detection Threshold")
    ax3.semilogx(det_info['Sxx_det_mean'][0], det_info['Sxx_det_mean'][1], '-b', linewidth=1.5, label="Detection (mean)")
    ax3.semilogx(det_info['Sxx_det_max'][0], det_info['Sxx_det_max'][1], '-r', linewidth=1.0, label="Detection (max)")
    ax3.axvspan(det_info['Freq Range'][0], det_info['Freq Range'][1], color='green', alpha=0.5, edgecolor=None)
    ax3.set_xlabel("Frequency [Hz]")
    ax3.set_ylabel("Power Spectral Density [Pa^2/Hz]")
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.set_ticks_position("right")

    fig.legend(loc="upper right")

    if title:
        ax2.set_title(title)
       
    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()