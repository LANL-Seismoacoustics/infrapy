#!/usr/bin/env python

# test.py
#
# Tutorial for ObsPy stream and trace manipulations
#
# Philip Blom (pblom@lanl.gov)
# Created 10/27/2022
# Last Modified 10/27/2022

import json
import gzip

import numpy as np
import matplotlib.pyplot as plt

from infrapy.detection import spectral, visualization
from infrapy.utils import config, data_io


if __name__ == '__main__':
    # ######################### #
    #     Define Parameters     #
    # ######################### #
    data_file = "data/YJ.BRP1..EDF.SAC"
    output_label = "spectral_test"
    cpu_cnt = 1

    sd_params = dict(config.defaults['SD'])
    sd_params["freq_min"] = 1.0
    sd_params["freq_max"] = 20.0
    sd_params["window_len"] = 900.0
    sd_params["window_step"] = 450.0
    sd_params["p_value"] = 0.01
    sd_params["freq_tm_factor"] = 35.0
    sd_params["cluster_eps"] = 10.0
    sd_params["cluster_min_samples"] = 40
    sd_params["cluster_window_len"] = 600.0

    # ##################### #
    #     Read data and     #
    #     run detector      #
    # ##################### #
    stream, latlon = data_io.set_stream(data_file, None, None)
    wvfrm_info = data_io.wvfrm_info(stream, latlon)

    det_list, spectrogram, history = spectral.spec_det_dict(stream[0], sd_params, cpu_cnt)
    det_output = {'wvfrm_info' : wvfrm_info, 'sd_params' : sd_params, 'spectrogram': spectrogram, 'history': history, 'det_info' : det_list}
    with gzip.open(output_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

    # ################### #
    #    Summarize and    #
    #     Plot Results    #
    # ################### #
    det_data = json.load(gzip.open(output_label + ".dets.json.gz", 'rt'))

    print("Detection summary:")
    for nd, det in enumerate(det_data['det_info']):
        spec_pnts = np.array(det['spec pnts'])
        t1, t2 = min(spec_pnts[:, 0]), max(spec_pnts[:, 0])
        f1, f2 = min(spec_pnts[:, 1]), max(spec_pnts[:, 1])

        print("   time: " + det['peak f-stat time'], end='\t')
        print("   duration [s]: " + str(np.round(t2 - t1, 2)), end='\t')
        print("   frequency range [Hz]: " + str(f1) + ", " + str(f2))

visualization.plot_sd_json(det_data, log_scale_freq=False, show_fig=False)
for k in range(len(det_list)):
    visualization.plot_sd_single_json(det_data, k, log_scale_freq=False, show_fig=False)
plt.show()

