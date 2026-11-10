#!/usr/bin/env python -W ignore::DeprecationWarning

# test_assoc.py
#
# Author    Philip Blom (pblom@lanl.gov)

import json
import gzip

from infrapy.event import hjl
from infrapy.propagation import infrasound
from infrapy.utils import config, data_io

if __name__ == '__main__':
    #########################
    ### Define parameters ###
    #########################

    # Read in detections from file
    dets_files = 'data/Blom_etal2020_GJI/SY*dets.json.gz'
    output_label = 'ev_bld_test'
    cpu_cnt = 4

    # define joint-likelihood calculation parameters
    ev_bld_params = dict(config.defaults['ASSOC'])

    ev_bld_params['celerity_model'] = 'regional_lf'
    ev_bld_params['back_az_width'] = 10.0
    ev_bld_params['range_max'] = 1500.0
    ev_bld_params['resolution'] = 180
    ev_bld_params['distance_matrix_max'] = 8.0
    ev_bld_params['cluster_threshold'] = 5.0
    ev_bld_params['trimming_threshold'] = 3.6
    ev_bld_params['ev_population_min'] = 3
    ev_bld_params['ev_station_min'] = 2

    for key in ev_bld_params.keys():
        if 'none' in str(ev_bld_params[key]).lower():
            ev_bld_params[key] = None


    ######################
    #### Run analysis ####
    ######################
    det_data = data_io._load_dets_json(dets_files)
    det_dicts = []
    for entry in det_data:
        for det in entry["det_info"]:
            det_dicts = det_dicts + [det]
            if 'fk_params' in entry.keys():
                det_dicts[-1]["wvfrm_info"] = entry["wvfrm_info"]
                det_dicts[-1]["fk_params"] = entry["fk_params"]
                det_dicts[-1]["det_params"] = entry["det_params"]


    det_list = [data_io._det_dict_to_likelihood(dict) for dict in det_dicts]
    infrasound._load_celerity_model(ev_bld_params['celerity_model'])
    events, ev_qls = hjl.id_events_dict(det_list, ev_bld_params, cpu_cnt)

    # Summarize clusters
    print("Identified " + str(len(events)) + " event(s).")
    for j, ev in enumerate(events):
        dist_mat = hjl.build_distance_matrix([det_list[k] for k in ev], bm_width=ev_bld_params['back_az_width'], rng_max=ev_bld_params['range_max'],
                                            rad_min=100.0, rad_max=(ev_bld_params['range_max'] / 4.0), resol=ev_bld_params['resolution'],  cpu_cnt=cpu_cnt, progress=False)

        ev_output = {'ground truth' : {}, 'det_info' : [det_dicts[k] for k in ev], 'assoc_params' : ev_bld_params, 'dist_matrix' : dist_mat, 'location' : [], 'characterization' : []}
        with gzip.open(output_label + "-" + str(j) + ".ev.json.gz", 'wt', encoding='UTF-8') as zipfile:
            json.dump(ev_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

        print('\nEvent (' + str(j) + ') deteciton summary:')
        ev_data = data_io._load_dets_json(output_label + "-" + str(j) + ".ev.json.gz")[0]
        for det in ev_data['det_info']:
            print(det['wvfrm_info'][0][0]['trace id'], end='\t')
            print("  loc: " + f"{det['wvfrm_info'][0][0]['latitude']:>7.4f}" + ", " + f"{det['wvfrm_info'][0][0]['longitude']:>8.4f}" , end='\t')
            print("  time: " + det['peak f-stat time'].split(".")[0], end='\t')
            print("  back az [deg]: " + f"{det["back az"]:>8.2f}", end='\t')
            print("  tr vel [m/s]: " + f"{det["tr vel"]:.1f}" , end='\t')
            print("  f-stat: " + f"{det["f-stat"]:>5.1f}")





