# test_yield.py
#
# Test computation of the yield using Spectral Yield Estimation (SpYE)
#
# Author:   Philip Blom (pblom@lanl.gov)
#

import os
import sys
import json
import gzip
import fnmatch

import numpy as np

from importlib.util import find_spec

from infrapy.characterization import spye
from infrapy.location import visualization
from infrapy.utils import config, data_io

if find_spec('stochprop'):
    import stochprop.propagation as sp_prop

if __name__ == '__main__':

    # ######################### #
    #     Define Parameters     #
    # ######################### #
    ev_file = "data/HRR-5.ev.json.gz"

    char_params = dict(config.defaults['YIELD'])
    char_params['det_mask'] = '0, 0, 1, 1, 0, 0, 1, 1'
    char_params['tlm_label'] = 'data/TLMs/2007_08-'

    char_params['freq_min'] = 0.2
    char_params['freq_max'] = 1.0
    char_params['yld_min'] = 1.0
    char_params['yld_max'] = 1.0e3
    char_params['resolution'] = 200

    char_params['ref_rng'] = 1.0
    char_params['amb_press'] = 101.325
    char_params['amb_temp'] = 288.15
    char_params['grnd_burst'] = True

    char_params['det_mask'] = [int(val) for val in char_params['det_mask'].strip(' ()[]').split(',')]

    for key in char_params:
        print(key, '\t', char_params[key])

    # #################### #
    #     Run Analysis     #
    # #################### #
    ev_data = data_io._load_dets_json(ev_file)[0]
    if len(ev_data['location']) > 0:
        src_loc = [ev_data['location'][0]['result']['lat_mean'],
                   ev_data['location'][0]['result']['lon_mean']]
    else:
        print("Location needed for transmission loss estimation.  Run test_bisl.py first.")
        sys.exit()

    det_info = [dict for k, dict in enumerate(ev_data['det_info']) if bool(char_params['det_mask'][k])]
    det_list = [data_io._det_dict_to_likelihood(dict) for dict in det_info]

    det_specs = [spye.extract_json_spectra(det['spec']) for det in det_info]
    det_specs = [np.array([spec[0], 10.0 * np.log10(spec[1])]) for spec in det_specs]

    # Load TLoss Models
    if not find_spec('stochprop'):
        print('\n' + "Can't load transmission loss model (TLM) file without stochprop installed." + '\n')
        sys.exit()

    print('\nLoading transmission loss statistics...')
    tlm_dir = os.path.dirname(char_params['tlm_label'])
    tlm_pattern = char_params['tlm_label'].split("/")[-1]
    tlm_files = [file_name for file_name in np.sort(os.listdir(tlm_dir)) if fnmatch.fnmatch(file_name, tlm_pattern + "*")]

    models = [0] * 2
    models[0] = [float(file_name.split("Hz")[0][len(tlm_pattern):]) for file_name in tlm_files]
    models[1] = [0] * len(tlm_files)
    for n in range(len(tlm_files)):
        print("  " + tlm_files[n])
        models[1][n] = sp_prop.TLossModel()
        models[1][n].load(tlm_dir + "/" + tlm_files[n], verbose=False)
    print('')

    spye_result = spye.run(det_list, det_specs, src_loc, np.array([char_params['freq_min'], char_params['freq_max']]), models,
                        yld_rng=np.array([char_params['yld_min'] * 1.0e3, char_params['yld_max'] * 1.0e3]),
                        ref_src_rng=char_params['ref_rng'], resol=char_params['resolution'], grnd_brst= char_params['grnd_burst'],
                        p_amb= char_params['amb_press'], T_amb= char_params['amb_temp'], exp_type= char_params['exp_type'])

    ev_data['characterization'] = ev_data['characterization'] + [{'params' : char_params, 'result' : spye_result}]
    with gzip.open(ev_file, 'wt', encoding='UTF-8') as zipfile:
        json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

    print('\n' + 'Results Summary (tons eq. TNT):')
    print('\t' + "Maximum a Posteriori Yield: " + str(spye_result['yld_vals'][np.argmax(spye_result['yld_pdf'])]))
    print('\t' + "68% Confidence Bounds: " + str(spye_result['conf_bnds'][0]))
    print('\t' + "95% Confidence Bounds: " + str(spye_result['conf_bnds'][1]))
    print('')

    visualization.plot_characterization(ev_data["det_info"], ev_data['location'][0], {'params' : char_params, 'result' : spye_result}, {},
                                        range_max=1000.0, show_fig=True)






