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
import matplotlib.ticker as mtick

from obspy.core import read

from scipy import signal

from infrapy.detection import beamforming_new

if __name__ == '__main__':
    # ######################### #
    #     Define Parameters     #
    # ######################### #
    sac_glob = "data/*.SAC"

    freq_min, freq_max = 0.5, 2.5
    window_length, window_step = 10.0, 2.5

    ns_start, ns_end = 100.0, 400.0
    sig_start, sig_end = 600, 800

    back_az_vals = np.arange(-180.0, 180.0, 1.5)
    trc_vel_vals = np.arange(300.0, 600.0, 2.5)

    method="bartlett"

    #p = mp.ProcessingPool(cpu_count())
    p = Pool(cpu_count() - 1)
    # ######################### #
    #  Read, Shift Start Time,  #
    #      and Filter Data      #
    # ######################### #
    x, t, t0, geom = beamforming_new.stream_to_array_data(read(sac_glob))
    M, N = x.shape

    # ######################### #
    #         View Data         #
    # ######################### #
    plt.figure(1)
    for m in range(M):
        plt.subplot(M, 1, m + 1)
        plt.xlim([0, t[-1]])
        plt.plot(t, x[m], 'k-')
        plt.axvspan(xmin = sig_start , xmax = sig_end, alpha = 0.25, color = 'blue')
        if method == "gls":
            plt.axvspan(xmin = ns_start , xmax = ns_end, alpha = 0.25, color = 'red')
        if m < (M - 1) : plt.setp(plt.subplot(M, 1, m + 1).get_xticklabels(), visible=False)

    if method == "gls":
        plt.suptitle("Data windows for signal (blue) and noise (red) \n Filtered in frequency range: " + str(freq_min) + " - " + str(freq_max) + "  Hz \n ")
    else:
        plt.suptitle("Data window for analysis \n Filtered in frequency range: " + str(freq_min) + " - " + str(freq_max) + "  Hz \n ")

    plt.show(block=False)
    plt.pause(0.1)

    # ######################### #
    #        Run Methods        #
    # ######################### #

    # define slowness and delays
    slowness = beamforming_new.build_slowness(back_az_vals, trc_vel_vals)
    delays = beamforming_new.compute_delays(geom, slowness)

    # define the noise covariance if using generalized least squares method
    if method == "gls":
        _, S, _ = beamforming_new.fft_array_data(x, t, window=[ns_start, ns_end], sub_window_len=window_length)

        ns_covar_inv = np.empty_like(S)
        for n in range(S.shape[2]):
            S[:, :, n] += 1.0e-3 * np.mean(np.diag(S[:, :, n])) * np.eye(S.shape[0])
            ns_covar_inv[:, :, n] = np.linalg.inv(S[:, :, n])
    else:
        ns_covar_inv = None

    # Prep figure
    f, a = plt.subplots(4, sharex=True)
    plt.xlim([sig_start, sig_end])
    a[3].set_xlabel("Time [s]")
    a[3].set_ylabel("Pr. [Pa]")
    a[2].set_ylabel("Back Az. [deg.]")
    a[1].set_ylabel("Tr. Vel. [m/s]")
    if method == "music":
        a[0].set_ylabel("Beam Power")
    else:
        a[0].set_ylabel("log10(F-value)")

    a[3].plot(t, x[1,:], '-k')
    plt.suptitle("Frequency range: " + str(freq_min) + " - " + str(freq_max) + " Hz \n window size " + str(window_length) + " seconds, window step " + str(window_step) +  " seconds")
    plt.show(block=False)

    # Run beamforming in windowed data and write to file
    times, beam_results = [],[]
    for window_start in np.arange(sig_start, sig_end, window_step):
        if window_start + window_length > sig_end:
            break

        print("Running analysis on time window " + str(window_start) + " - " + str(window_start + window_length), end=' ')
        print(" seconds in frequency band " + str(freq_min) + " - " + str(freq_max) + " Hz...")
        times = times + [[t0 + np.timedelta64(int(window_start), 's')]]
        X, S, f = beamforming_new.fft_array_data(x, t, window=[window_start, window_start + window_length])
        beam_power = beamforming_new.run(X, S, f, geom, delays, [freq_min, freq_max], method="bartlett", pool=p, normalize_beam=True, ns_covar_inv=ns_covar_inv)
        peaks = beamforming_new.find_peaks(beam_power, back_az_vals, trc_vel_vals, signal_cnt=1)
        beam_results = beam_results + [[peaks[0][0], peaks[0][1], peaks[0][2] / (1.0 - peaks[0][2]) * (x.shape[0] - 1)]]

        if method == "music":
            a[2].plot([window_start + 1.0 / 2.0 * window_length], [peaks[0][0]], 'ok', markersize=3.3)
            a[1].plot([window_start + 1.0 / 2.0 * window_length], [peaks[0][1]], 'ok', markersize=3.3)
            a[0].plot([window_start + 1.0 / 2.0 * window_length], [peaks[0][2]], 'ok', markersize=3.3)
            plt.pause(0.1)
        else:
            fisher_val = peaks[0][2] / (1.0 - peaks[0][2]) * (M - 1)
            a[2].plot([window_start + 1.0 / 2.0 * window_length], [peaks[0][0]], 'ok', markersize=3.3)
            a[1].plot([window_start + 1.0 / 2.0 * window_length], [peaks[0][1]], 'ok', markersize=3.3)
            a[0].plot([window_start + 1.0 / 2.0 * window_length], [np.log10(fisher_val)], 'ok', markersize=3.3)
            plt.pause(0.1)
        plt.pause(0.1)
    times = np.array(times)[:, 0]
    beam_results = np.array(beam_results)

    np.save("data/times", times)
    np.save("data/beam_results", beam_results)
    # Define best beam time series and residuals
    back_az = beam_results[np.argmax(beam_results[:, 2]), 0]
    tr_vel = beam_results[np.argmax(beam_results[:, 2]), 1]

    X, S, f = beamforming_new.fft_array_data(x, t, window=[sig_start, sig_end], fft_window="boxcar")
    sig_est, residual = beamforming_new.extract_signal(X, f, np.array([back_az, tr_vel]), geom)

    plt.figure(3)
    plt.loglog(f, abs(sig_est), '-b', linewidth=1.0)
    plt.loglog(f, np.mean(abs(residual), axis=0), '-k', linewidth=0.5)

    signal_wvfrm = np.fft.irfft(sig_est) / (t[1] - t[0])
    resid_wvfrms = np.fft.irfft(residual, axis=1) / (t[1] - t[0])
    t_mask = np.logical_and(sig_start < t, t < sig_end)

    plt.figure(4)
    for m in range(M):
        plt.subplot(M + 1, 1, m + 1)
        plt.xlim([t[t_mask][0], t[t_mask][-1]])
        plt.plot(t[t_mask], x[m, t_mask], '0.5')
        plt.plot(t[t_mask], resid_wvfrms[m, :len(t[t_mask])], 'k-')
        plt.setp(plt.subplot(M + 1, 1, m + 1).get_xticklabels(), visible=False)
    plt.subplot(M + 1, 1, M + 1)
    plt.xlim([t[t_mask][0], t[t_mask][-1]])
    plt.plot(t[t_mask], signal_wvfrm[:len(t[t_mask])], 'b-')
    plt.pause(30.0)

    plt.close()
    p.close()
