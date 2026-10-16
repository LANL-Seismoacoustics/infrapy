#!/usr/bin/env python

# test_bisl.py
#
# Author    Philip Blom (pblom@lanl.gov)
import json
import gzip

from infrapy.location import bisl, visualization
from infrapy.propagation import infrasound
from infrapy.utils import config, data_io



if __name__ == '__main__':
    # ######################### #
    #       Define Inputs       #
    # ######################### #

    # Define ground_truth if known (41.131, -112.896 for UTTR; Test includes show in June 2004)
    ev_file = 'data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz'

    # Define localization parameters
    loc_params = dict(config.defaults['LOC'])
    loc_params['celerity_model'] = 'regional_lf'
    loc_params['back_az_width'] = 10.0
    loc_params['range_max'] = 1500.0
    loc_params['grid_resol'] = 180

    for key in loc_params.keys():
        if 'none' in str(loc_params[key]).lower():
            loc_params[key] = None

    # Load detection list from json file
    ev_data = data_io._load_dets_json(ev_file)[0]
    det_list = [data_io._det_dict_to_likelihood(dict) for dict in ev_data['det_info']]
    result = bisl.run_dict(det_list, loc_params, None, None)

    print('\n' + "Localization Summary:")
    print(bisl.summarize(result))
