#!/usr/bin/env python  -W ignore::DeprecationWarning

# test.py
#
# Tutorial for ObsPy stream and trace manipulations
#
# Philip Blom (pblom@lanl.gov)
# Created 01/06/2016
# Last Modified 01/06/2016

import json
import gzip

import numpy as np
import matplotlib.pyplot as plt

from infrapy.detection import beam, visualization
from infrapy.utils import config, data_io

if __name__ == '__main__':

    # ######################### #
    #     Define Parameters     #
    # ######################### #
    sac_glob = "data/*.SAC"
    output_label = "beam_test"
    cpu_cnt = 4

    fk_params = dict(config.defaults['FK'])
    fk_params['freq_min'] = 0.5
    fk_params['freq_max'] = 2.5

    fk_params['window_len'] = 10.0
    fk_params['window_step'] = 2.5

    fk_params['back_az_min'] = -180.0
    fk_params['back_az_max'] = 180.0
    fk_params['back_az_step'] = 2.5

    fk_params['trace_vel_min'] = 300.0
    fk_params['trace_vel_max'] = 600.0
    fk_params['trace_vel_step'] = 2.5

    det_params = dict(config.defaults['FD'])
    det_params['window_len'] = 3600.0
    det_params['p_value'] = 0.02
    det_params['back_az_width'] = 15.0
    det_params['min_duration'] = 20.0
    det_params['merge_dets'] = True

    # convert string 'None' to Python None for logic flow
    for key in det_params.keys():
        if 'none' in str(det_params[key]).lower():
            det_params[key] = None

    # ######################### #
    #  Read Data and Summarize  #
    # ######################### #
    stream, latlon = data_io.set_stream(sac_glob, None, None)
    wvfrm_info = data_io.wvfrm_info(stream, latlon)

    print('\n' + "Data summary:")
    for tr in stream:
        print(tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    # ####################### #
    #     Run Beamforming     #
    # ####################### #
    beam_times, beam_peaks = beam.run_fk_dict(stream, latlon, fk_params, None, cpu_cnt)
    dets, thresh_vals = beam.run_afd_dict(beam_times, beam_peaks, fk_params, det_params, len(stream))

    # extract fk results for the full duration
    fk = {}
    fk['time'] = np.array([(tn - np.datetime64(stream[0].stats.starttime)).astype('m8[ms]').astype(float) * 1.0e-3 for tn in beam_times])
    fk['back az'] = beam_peaks[:, 0]
    fk['tr vel'] = beam_peaks[:, 1]
    fk['f-stat'] = beam_peaks[:, 2]
    fk['thresh'] = thresh_vals

    dets = [beam.det2dict(stream, latlon, beam_times, beam_peaks, fk_params, det_info) for det_info in dets]
    dets_dict = {'wvfrm_info' : [wvfrm_info],
                 'fk_params' : [fk_params],
                 'det_params' : [det_params],
                 'fk' : fk,
                 'det_info' : dets}

    with gzip.open(output_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(dets_dict, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

    # ################### #
    #    Summarize and    #
    #     Plot Results    #
    # ################### #
    det_data = json.load(gzip.open(output_label + ".dets.json.gz", 'rt'))

    print('\n' + "Detection summary:")
    for k, det in enumerate(det_data['det_info']):
        print("  time: " + str(det['peak f-stat time']).split(".")[0], end='\t')
        print("  back az [deg]: " + f"{det["back az"]:>8.2f}", end='\t')
        print("  tr vel [m/s]: " + f"{det["tr vel"]:.1f}", end='\t')
        print("  f-stat: " + f"{det["f-stat"]:>5.1f}")
    print('')

    # Plot beam results
    visualization.plot_fk_json(det_data, show_fig=False)
    for k in range(len(dets)):
        visualization.plot_det_json(det_data, k, show_fig=False)
    plt.show()