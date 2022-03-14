# test_yield.py
#
# Test computation of the yield using Spectral Yield Estimation (SpYE)
#
# Author:   Philip Blom (pblom@lanl.gov)
#

from obspy.core import read

import numpy as np

import matplotlib.pyplot as plt 

from infrapy.propagation import likelihoods as lklhds
from infrapy.propagation import infrasound

from infrapy.characterization import spye

if __name__ == '__main__':
    # ######################### #
    #     Define Parameters     #
    # ######################### #
    det_file = "data/HRR-5.dets.json"
    data_path = "../infrapy-data/hrr-5/"
    data_ids = ["W220/HR5.W220*.sac", "W240/HR5.W240*.sac", "W340/HR5.W340*.sac", "W420/HR5.W420*.sac", "W460/HR5.W460*.sac"]

    ns_opt = "post"
    win_buffer = 0.2
    
    src_loc = np.array([33.5377, -106.333961])
    freq_band = np.array([0.25, 2.0])
    yld_rng = np.array([1.0e3, 1000.0e3])
    ref_rng = 1.0

    grnd_truth=None
    resol = 200

    # ############################# #
    #     Define the detections     #
    #          and spectra          #
    # ############################# #
    det_list = lklhds.json_to_detection_list(det_file)
    st_list = [0] * len(det_list)
    for j in range(len(st_list)):
        st_list[j] = read(data_path + data_ids[j] )
    smn_specs = spye.extract_spectra(det_list, st_list, win_buffer=win_buffer, ns_opt=ns_opt)
    
    # ######################### #
    #     Load TLoss Models     #
    # ######################### #
    tloss_f_min, tloss_f_max, tloss_f_cnt = 0.025, 2.5, 25

    models = [0] * 2
    models[0] = list(np.logspace(np.log10(tloss_f_min), np.log10(tloss_f_max), tloss_f_cnt))
    models[1] = [0] * tloss_f_cnt
    for n in range(tloss_f_cnt):
        models[1][n] = infrasound.TLossModel()
        models[1][n].load("../infrapy/propagation/priors/tloss/2007_08-" + "%.3f" % models[0][n] + "Hz.pri")

    # ######################## #
    #         Run Yield        #
    #    Estimation Methods    #
    # ######################## #
    yld_vals, yld_pdf, conf_bnds = spye.run(det_list, smn_specs, src_loc, freq_band, models, yld_rng=yld_rng, ref_src_rng=ref_rng, resol=resol)

    print('\nResults:')
    print('\t' + "Maximum a Posteriori Yield:", yld_vals[np.argmax(yld_pdf)])
    print('\t' + "68% Confidence Bounds:", conf_bnds[0])
    print('\t' + "95% Confidence Bounds:", conf_bnds[1])

    plt.semilogx(yld_vals, yld_pdf)
    plt.fill_between(yld_vals, yld_pdf, where=np.logical_and(conf_bnds[0][0] <= yld_vals, yld_vals <= conf_bnds[0][1]), color='g', alpha=0.25)
    plt.fill_between(yld_vals, yld_pdf, where=np.logical_and(conf_bnds[1][0] <= yld_vals, yld_vals <= conf_bnds[1][1]), color='g', alpha=0.25)

    plt.show(block=False)
    plt.pause(5.0)
    plt.close()

