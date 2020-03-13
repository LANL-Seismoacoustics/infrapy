#!/usr/bin/env python  -W ignore::DeprecationWarning

# test.py
#
# Tutorial for ObsPy stream and trace manipulations
#
# Philip Blom (pblom@lanl.gov)
# Created 01/06/2016
# Last Modified 01/06/2016

import numpy as np

import pathos.multiprocessing as mp
from multiprocessing import cpu_count
from multiprocess import Pool
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from obspy.core import read

from infrapy.detection import beamforming_new

if __name__ == '__main__':
    # ######################### #
    #     Define Parameters     #
    # ######################### #
    sac_glob = "data/*.SAC"

    freq_min, freq_max = 0.5, 2.5
    window_length, window_step = 5.0, 2.5

    ns_start, ns_end = 100.0, 400.0
    sig_start, sig_end = 600, 800

    back_az_vals = np.arange(-180.0, 180.0, 1.5)
    trc_vel_vals = np.arange(250.0, 600.0, 2.5)

    method="capon"

    #p = mp.ProcessingPool(cpu_count() - 1)
    p = Pool(cpu_count() - 1)

    # ######################### #
    #         Read and          #
    #        Filter Data        #
    # ######################### #
    x, t, t0, geom = beamforming_new.stream_to_array_data(read(sac_glob))
    M, N = x.shape

    # ######################### #
    #         View Data         #
    # ######################### #
    plt.figure(1)
    for m in range(M):
        plt.subplot(M, 1, m + 1)
        plt.xlim([sig_start, sig_end])
        plt.plot(t, x[m], '-k')
        if m < (M - 1) : plt.setp(plt.subplot(M, 1, m + 1).get_xticklabels(), visible=False)
    plt.show(block=False)

    plt.figure(2)
    plt.show(block=False)

    plt.figure(3)
    plt.show(block=False)

    window_shade = list(range(M))
    palette = cm.jet

    # ######################### #
    #        Run Methods        #
    # ######################### #
    # define slowness
    slowness = beamforming_new.build_slowness(back_az_vals, trc_vel_vals)
    delays = beamforming_new.compute_delays(geom, slowness)

    # define the noise covariance if using generalized least squares method
    if method == "gls":
        _, _, S = beamforming_new.fft_array_data(x, t, window=[ns_start, ns_end], sub_window_len=window_length)
        for n in range(len(S)):
            S[n] += 1.0e-3 * np.mean(np.diag(S[n])) * np.eye(S[n].shape[0])
        ns_covar_inv = np.linalg.inv(S)
    else:
        ns_covar_inv = None

    for window_start in np.arange(sig_start, sig_end - window_length, window_step):
        print("Running analysis on time window " + str(window_start) + " - " + str(window_start + window_length) + " seconds in frequency band " + str(freq_min) + " - " + str(freq_max) + " Hz...")
        plt.figure(1)
        for m in range(M):
            window_shade[m] = plt.subplot(M, 1, m + 1).axvspan(xmin = window_start, xmax = window_start + window_length, alpha = 0.33, color = 'blue')
        plt.pause(0.1)

        X, S, f = beamforming_new.fft_array_data(x, t, window=[window_start, window_start + window_length])
        beam_power = beamforming_new.run(X, S, f, geom, delays, [freq_min, freq_max], method=method, signal_cnt=1, pool=p, ns_covar_inv=ns_covar_inv, normalize_beam=True)

        avg_beam_power = np.average(beam_power, axis=0)
        #avg_beam_power = beamforming_new.multi_freq_beam(beam_power)

        plt.figure(2)
        plt.clf()
        plt.xlim([min(slowness[:, 0]), max(slowness[:, 0])])
        plt.ylim([min(slowness[:, 1]), max(slowness[:, 1])])
        if method == "bartlett_covar" or method == "bartlett" or method == "gls":
            plt.scatter(slowness[:, 0], slowness[:, 1], c=avg_beam_power, cmap=palette, marker="o", s=[12.5] * len(slowness), edgecolor='none', vmin=0.0, vmax=1.0)
        else:
            plt.scatter(slowness[:, 0], slowness[:, 1], c=avg_beam_power, cmap=palette, marker="o", s=[12.5] * len(slowness), edgecolor='none', vmin=0.0, vmax=np.max(avg_beam_power))
        plt.pause(1.0)

        # Compute back azimuth projection of distribution
        az_proj, tv_proj = beamforming_new.project_beam(beam_power, back_az_vals, trc_vel_vals, method="mean")

        plt.figure(3)
        plt.clf()
        plt.xlim([min(back_az_vals), max(back_az_vals)])
        if method == "bartlett_covar" or method == "bartlett" or method == "gls":
            plt.ylim([0.0, 1.0])
        else:
            plt.ylim([0.0, np.max(avg_beam_power)])
        plt.plot(back_az_vals, az_proj, '-k', linewidth=2.5)
        plt.pause(0.2)

        for m in range(M):
            window_shade[m].remove()

    plt.close()
    p.close()
