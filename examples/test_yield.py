# test_yield.py
#
# Test computation of the yield using Spectral Yield Estimation (SpYE)
#
# Author:   Philip Blom (pblom@lanl.gov)
#

from obspy import read, Stream

import numpy as np

import matplotlib.pyplot as plt 

from infrapy.utils import data_io
from infrapy.propagation import infrasound

from infrapy.characterization import spye

if __name__ == '__main__':
    # ######################### #
    #     Define Parameters     #
    # ######################### #
    det_file = "data/HRR-5.dets.json"
    wvfrm_path = "../infrapy-data/hrr-5/*/*.sac"
    tloss_path = "../infrapy/propagation/priors/tloss/2007_08-"

    ns_opt = "post"
    win_buffer = 0.2
    
    src_loc = np.array([33.5377, -106.333961])
    freq_band = np.array([0.25, 2.0])
    yld_rng = np.array([1.0e3, 1000.0e3])
    ref_rng = 1.0

    grnd_truth = None
    resol = 200

    # ############################# #
    #     Define the detections     #
    #          and spectra          #
    # ############################# #
    det_list = data_io.json_to_detection_list(det_file)
    st_list = [Stream([tr for tr in read(wvfrm_path) if det.station in tr.stats.station]) for det in det_list]
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
        models[1][n].load(tloss_path + "%.3f" % models[0][n] + "Hz.pri")

    # ######################## #
    #         Run Yield        #
    #    Estimation Methods    #
    # ######################## #
    yld_results = spye.run(det_list, smn_specs, src_loc, freq_band, models, yld_rng=yld_rng, ref_src_rng=ref_rng, resol=resol)

    print('\nResults:')
    print('\t' + "Maximum a Posteriori Yield:", yld_results['yld_vals'][np.argmax(yld_results['yld_pdf'])])
    print('\t' + "68% Confidence Bounds:", yld_results['conf_bnds'][0])
    print('\t' + "95% Confidence Bounds:", yld_results['conf_bnds'][1])

    plt.semilogx(yld_results['yld_vals'], yld_results['yld_pdf'])
    plt.fill_between(yld_results['yld_vals'], yld_results['yld_pdf'], where=np.logical_and(yld_results['conf_bnds'][0][0] <= yld_results['yld_vals'], yld_results['yld_vals'] <= yld_results['conf_bnds'][0][1]), color='g', alpha=0.25)
    plt.fill_between(yld_results['yld_vals'], yld_results['yld_pdf'], where=np.logical_and(yld_results['conf_bnds'][1][0] <= yld_results['yld_vals'], yld_results['yld_vals'] <= yld_results['conf_bnds'][1][1]), color='g', alpha=0.25)

    plt.show()

