# visualization.py
#
# Visualization methods for beamforming (fk) and detection (fd) results
#
# Philip Blom (pblom@lanl.gov)


import numpy as np

import matplotlib.pyplot as plt 
from matplotlib import cm

from . import beamforming_new


def plot_fk1(stream, latlon, times, peaks, detections=None, title=None, output_path=None, show_fig=True, det_thresh=None):
    '''
    Visualize beamforming (fk) results with waveform data included

    '''   
    x, t, t0, _ = beamforming_new.stream_to_array_data(stream, latlon)

    f, a = plt.subplots(4, sharex=True)
    a[3].set_xlabel("Time")
    a[0].set_ylabel("log10(F-value)")    
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

    f, a = plt.subplots(3, sharex=True)
    a[2].set_xlabel("Time")
    a[0].set_ylabel("log10(F-value)")    
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


