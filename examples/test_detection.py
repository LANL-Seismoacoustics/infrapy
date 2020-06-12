#!/usr/bin/env python  -W ignore::DeprecationWarning

# test_detection.py
#
# Tutorial for detection utilizing the AFD on a series of beamforming results, saved in a file as
#
# @fkdd (fransiska at lanl dot gov)
# @pblom (pblom at lanl dot gov)
# Last Modified 12/19/2019

import numpy as np

import pathos.multiprocessing as mp
from multiprocessing import cpu_count

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from obspy.core import read

from scipy import signal

from infrapy.detection import beamforming_new

# ######################### #
#     Define Parameters     #
# ######################### #

# Detection params
# times_file, beam_results_file = None, None
times_file, beam_results_file = "data/times.npy", "data/beam_results.npy"

det_win_len = 60 * 5
det_thresh = 0.99
min_seq = 5
TB_prod = 40 * 10
back_az_lim = 10
channel_cnt = 4

if __name__ == '__main__':
    ######################################
    ##  Load data and prepare analysis  ##
    ######################################

    if times_file and beam_results_file:
        times = np.load(times_file)
        beam_results = np.load(beam_results_file)
    else:
        print('No beamforming input provided')

    ######################################
    ##      Run detection analysis      ##
    ######################################
    dets = beamforming_new.detect_signals(times, beam_results, det_win_len, TB_prod, channel_cnt, det_thresh=det_thresh, min_seq=min_seq, back_az_lim=back_az_lim)

    print('\n' + "Detection Summary:")
    for det in dets:
        print("Detection time:", det[0], '\t', "Rel. detection onset:", det[1], '\t',"Rel. detection end:", det[2], '\t',end=' ')
        print("Back azimuth:", np.round(det[3], 2), '\t', "Trace velocity:", np.round(det[4], 2), '\t', "F-stat:", np.round(det[5], 2), '\t', "Array dim:", channel_cnt)
