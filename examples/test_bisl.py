#!/usr/bin/env python

# test_bisl.py
#
# Author    Philip Blom (pblom@lanl.gov)
import json
import gzip

from infrapy.event import bisl, visualization
from infrapy.propagation import infrasound
from infrapy.utils import config, data_io



if __name__ == '__main__':
    # ######################### #
    #       Define Inputs       #
    # ######################### #

    # Event file from Blom et al. (2024) is used to test BISL
    ev_file = "data/HRR-5.ev.json.gz"

    # Define localization parameters
    loc_params = dict(config.defaults['LOC'])
    loc_params['celerity_model'] = 'regional_lf'
    loc_params['back_az_width'] = 10.0
    loc_params['range_max'] = 1000.0
    loc_params['grid_resol'] = 180

    for key in loc_params.keys():
        if 'none' in str(loc_params[key]).lower():
            loc_params[key] = None

    # Load detection list from json file
    ev_data = data_io._load_dets_json(ev_file)[0]
    det_list = [data_io._det_dict_to_likelihood(dict) for dict in ev_data['det_info']]

    infrasound._load_celerity_model(loc_params["celerity_model"])
    result = bisl.run_dict(det_list, loc_params, None, None)

    ev_data['location'] = [{'params' : loc_params, 'result' : result}]
    with gzip.open(ev_file, 'wt', encoding='UTF-8') as zipfile:
        json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

    ev_data = data_io._load_dets_json(ev_file)[0]
    loc_info = ev_data['location'][0]

    print('\nLocalization Summary:\n' + bisl.summarize(loc_info['result']))
    visualization.plot_localization(det_list, loc_info, {}, range_max=loc_params['range_max'], show_fig=True)

