#!/usr/bin/env python
import os 
import click
import json
import gzip
import fnmatch
import tempfile
import warnings

import configparser as cnfg
import numpy as np

from multiprocessing import Pool
from importlib.util import find_spec

import matplotlib.pyplot as plt 

from ..utils import config
from ..utils import data_io

from ..association import hjl

from ..location import bisl, tribl
from ..propagation import infrasound

from ..characterization import spye


@click.command('build', short_help="Associate detections into events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--detect-files", help="Detection path and pattern", default=None)
@click.option("--event-label", help="Path for event info output", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--celerity-model", help="Celerity model option or file (default: '" + config.defaults['ASSOC']['celerity_model'], default=None)
@click.option("--back-az-width", help="Width of beam projection (default: " + config.defaults['ASSOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Maximum source-receiver range (default: " + config.defaults['ASSOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['ASSOC']['resolution'] + ")", default=None, type=int)
@click.option("--distance-matrix-max", help="Distance matrix maximum (default: " + config.defaults['ASSOC']['distance_matrix_max'] + ")", default=None, type=float)
@click.option("--cluster-linkage", help="Linkage method for clustering (default: " + config.defaults['ASSOC']['cluster_linkage'] + ")", default=None)
@click.option("--cluster-threshold", help="Cluster linkage threshold (default: " + config.defaults['ASSOC']['cluster_threshold'] + ")", default=None, type=float)
@click.option("--trimming-threshold", help="Mishapen cluster threshold (default: " + config.defaults['ASSOC']['trimming_threshold'] + ")", default=None, type=float)
@click.option("--event-population-min", help="Minimum detection count in event (default: " + config.defaults['ASSOC']['event_population_min'] + ")", default=None, type=int)
@click.option("--event-station-min", help="Minimum station count in event (default: " + config.defaults['ASSOC']['event_station_min'] + ")", default=None, type=int)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def build(config_file, detect_files, event_label, starttime, endtime, celerity_model, back_az_width, range_max, resolution, distance_matrix_max,  
                cluster_linkage, cluster_threshold, trimming_threshold, event_population_min, event_station_min, cpu_cnt):
    '''
    Run association analysis to identify events in a detection set

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy event build --detect-files 'data/Blom_etal2020_GJI/SY*dets.json.gz' --event-label Blom_etal2020_GJI --range-max 1500.0 --cpu-cnt 4
    '''

    click.echo("")
    click.echo("####################################")
    click.echo("##                                ##")
    click.echo("##             InfraPy            ##")
    click.echo("##         Event Building         ##")
    click.echo("##                                ##")
    click.echo("####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None


    # Data IO parameters
    detect_files = config.set_param(user_config, 'DETECTION IO', 'detect_files', detect_files, 'string')
    event_label = config.set_param(user_config, 'DETECTION IO', 'event_label', event_label, 'string')

    # Data IO parameters
    click.echo('\n' + "Data summary:")
    click.echo("  detect_files: " + str(detect_files))
    click.echo("  event_file: " + str(event_label))

    if detect_files is None or event_label is None:
        msg = "Association analysis requires detection input (--detect-files) and output path (--event-file)"
        warnings.warn(msg)
        return 0

    # Algorithm parameters
    assoc_params = {}
    assoc_params['starttime'] = config.set_param(user_config, 'ASSOC', 'starttime', starttime, 'string')
    assoc_params['endtime'] = config.set_param(user_config, 'ASSOC', 'endtime', endtime, 'string')

    assoc_params['celerity_model'] = config.set_param(user_config, 'ASSOC', 'celerity_model', celerity_model, 'string')

    assoc_params['back_az_width'] = config.set_param(user_config, 'ASSOC', 'back_az_width', back_az_width, 'float')
    assoc_params['range_max'] = config.set_param(user_config, 'ASSOC', 'range_max', range_max, 'float')
    assoc_params['resolution'] = config.set_param(user_config, 'ASSOC', 'resolution', resolution, 'int')

    assoc_params['distance_matrix_max'] = config.set_param(user_config, 'ASSOC', 'distance_matrix_max', distance_matrix_max, 'float')
    assoc_params['cluster_linkage'] = config.set_param(user_config, 'ASSOC', 'cluster_linkage', cluster_linkage, 'string')
    assoc_params['cluster_threshold'] = config.set_param(user_config, 'ASSOC', 'cluster_threshold', cluster_threshold, 'float')

    assoc_params['trimming_threshold'] = config.set_param(user_config, 'ASSOC', 'trimming_threshold', trimming_threshold, 'float')
    assoc_params['event_population_min'] = config.set_param(user_config, 'ASSOC', 'event_population_min', event_population_min, 'int')
    assoc_params['event_station_min'] = config.set_param(user_config, 'ASSOC', 'event_station_min', event_station_min, 'int')
    assoc_params['cpu_cnt'] = config.set_param(user_config, 'ASSOC', 'cpu_cnt', cpu_cnt, 'int')

    click.echo('\n' + "association parameters:")
    for key in assoc_params.keys():
        click.echo("  " + key + ": " + str(assoc_params[key]))

    if assoc_params['cpu_cnt'] is not None:
        pl = Pool(assoc_params['cpu_cnt'])
    else:
        pl = None

    click.echo("")

    det_data = data_io._load_dets_json(detect_files)
    det_dicts = []
    for entry in det_data:
        for det in entry["det_info"]:
            det_dicts = det_dicts + [det]
            if 'fk_params' in entry.keys():
                det_dicts[-1]["wvfrm_info"] = entry["wvfrm_info"]
                det_dicts[-1]["fk_params"] = entry["fk_params"]
                det_dicts[-1]["det_params"] = entry["det_params"]

    # Check if an event file exists with matching parameter configuration and detections from this list
    result_check = False
    if os.path.isfile(event_label + "-0.ev.json.gz"):
        ev0_data = data_io._load_dets_json(event_label + "-0.ev.json.gz")[0]

        param_check = False
        try:
            np.testing.assert_equal(assoc_params, ev0_data['assoc_params'])
            param_check = True 
        except:
            pass

        dets_check = True
        for ev_det in ev0_data['det_info']:
            check = np.any([np.all([ev_det[key] == list_det[key] for key in ['peak f-stat time', 'f-stat', 'back az', 'tr vel']]) for list_det in det_dicts])
            dets_check = dets_check and check 

        result_check = param_check and dets_check

        if not result_check:
            warning_message = "Event result exists, but parameters or detections don't match.  I hope you're meaning to overwrite existing results from another event building run!"
            warnings.warn((warning_message))
        
    if result_check:
        click.echo("Event results found matching this parameter configuration and detection set.  Skipping event building to avoid overwriting existing results.")
    else:
        det_list = [data_io._det_dict_to_likelihood(dict) for dict in det_dicts]

        infrasound._load_celerity_model(assoc_params['celerity_model'])
        events, event_qls = hjl.id_events(det_list, assoc_params['cluster_threshold'], starttime=assoc_params['starttime'], endtime=assoc_params['endtime'], dist_max=assoc_params['distance_matrix_max'], 
                                        bm_width=assoc_params['back_az_width'], rng_max=assoc_params['range_max'], rad_min=100.0, rad_max=(assoc_params['range_max'] / 4.0), 
                                        resol=assoc_params['resolution'], linkage_method=assoc_params['cluster_linkage'], trimming_thresh=assoc_params['trimming_threshold'], 
                                        cluster_det_population=assoc_params['event_population_min'], cluster_array_population=assoc_params['event_station_min'], pool=pl)

        click.echo("Identified " + str(len(events)) + " event(s)." + '\n')
        for j, ev in enumerate(events):
            dist_mat = hjl.build_distance_matrix([det_list[k] for k in ev], bm_width=assoc_params['back_az_width'], rng_max=assoc_params['range_max'],
                                                rad_min=100.0, rad_max=(assoc_params['range_max'] / 4.0), resol=assoc_params['resolution'],  pool=pl, progress=False)

            ev_output = {'ground truth' : {}, 'det_info' : [det_dicts[k] for k in ev], 'assoc_params' : assoc_params, 'dist_matrix' : dist_mat, 'location' : [], 'characterization' : []}
            with gzip.open(event_label + "-" + str(j) + ".ev.json.gz", 'wt', encoding='UTF-8') as zipfile:
                json.dump(ev_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

    if pl is not None:
        pl.terminate()
        pl.close()


@click.command('localize', short_help="Estimate source location and origin time")
@click.option("--event-file", help="Event JSON file to be analyzed", default=None)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--back-az-width", help="Width of beam projection (default: " + config.defaults['LOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)

@click.option("--grid-resol", help="Grid resolution (number of points) (default: " + config.defaults['LOC']['grid_resol'] + ")", default=None, type=int)

@click.option("--ll-corner", help="Lower left corner of region (lat, lon)", default=None)
@click.option("--ur-corner", help="Upper right corner of region (lat, lon)", default=None)
@click.option("--latlon-resol", help="Resolution of latitude/longitude grid (degrees)", default=None, type=float)

@click.option("--tm-min", help="Minimum origin time", default=None)
@click.option("--tm-max", help="Maximum origin time", default=None)
@click.option("--tm-resol", help="Resolution of origin time grid (seconds)", default=None, type=float)

@click.option("--celerity-model", help="Celerity model option or file (default: '" + config.defaults['LOC']['celerity_model'], default=None)
@click.option("--pgm-file", help="Path geometry model (PGM) file (optional)", default=None)

@click.option("--atmo-data", help="Atmosphere data if using TRIBL", default=None)
@click.option("--alt-lims", help="Altitude limits if using TRIBL", default=None)
@click.option("--alt-resol", help="Altitude resolution if using TRIBL", default=None)
@click.option("--grnd-snd-spd", help="Sound speed at the ground", default=None)
@click.option("--c0-stdev", help="Sound speed uncertainty if using TRIBL", default=None)
@click.option("--det-tm-stdev", help="Detection time uncertainty", default=None)
@click.option("--az-limit", help="Azimuth resolution limit", default=None)
@click.option("--local-temp-dir", help="Local temporary directory if using TRIBL", default=None)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)

def localize(event_file, config_file, back_az_width, range_max, grid_resol, ll_corner, ur_corner, latlon_resol, tm_min, tm_max, tm_resol, celerity_model, pgm_file, atmo_data, alt_lims, alt_resol, grnd_snd_spd, c0_stdev, det_tm_stdev, az_limit, local_temp_dir, cpu_cnt):
    '''
    Run Bayesian Infrasonic Source Localization (BISL) methods to estimate the source location and origin time for an event

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy event localize --event-file data/Blom_etal2020_GJI/Blom_etal2020_GJI-0.ev.json.gz
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##      Localization Analysis      ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Data IO parameters
    event_file = config.set_param(user_config, 'DATA IO', 'event_file', event_file, 'string')
            
    click.echo('\n' + "Data summary:")
    click.echo("  event_file: " + str(event_file))

    ev_data = data_io._load_dets_json(event_file)[0]
    if range_max is None:
        range_max = ev_data['assoc_params']['range_max']

    # Algorithm parameters
    loc_params = {}
    loc_params['back_az_width'] = config.set_param(user_config, 'LOC', 'back_az_width', back_az_width, 'float')
    loc_params['range_max'] = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    loc_params['grid_resol'] = config.set_param(user_config, 'LOC', 'grid_resol', grid_resol, 'int')

    loc_params['ll_corner'] = config.set_param(user_config, 'LOC', 'll_corner', ll_corner, 'str')
    loc_params['ur_corner'] = config.set_param(user_config, 'LOC', 'ur_corner', ur_corner, 'str')
    loc_params['latlon_resol'] = config.set_param(user_config, 'LOC', 'latlon_resol', latlon_resol, 'float')

    loc_params['tm_min'] = config.set_param(user_config, 'LOC', 'tm_min', tm_min, 'str')
    loc_params['tm_max'] = config.set_param(user_config, 'LOC', 'tm_max', tm_max, 'str')
    loc_params['tm_resol'] = config.set_param(user_config, 'LOC', 'tm_resol', tm_resol, 'float')

    loc_params['celerity_model'] = config.set_param(user_config, 'LOC', 'celerity_model', celerity_model, 'str')
    loc_params['pgm_file'] = config.set_param(user_config, 'LOC', 'pgm_file', pgm_file, 'str')

    loc_params['atmo_data'] = config.set_param(user_config, 'LOC', 'atmo_data', atmo_data, 'str')
    loc_params['alt_lims'] = config.set_param(user_config, 'LOC', 'alt_lims', alt_lims, 'str')
    loc_params['alt_resol'] = config.set_param(user_config, 'LOC', 'alt_resol', alt_resol, 'float')
    loc_params['grnd_snd_spd'] = config.set_param(user_config, 'LOC', 'grnd_snd_spd', grnd_snd_spd, 'float')
    loc_params['c0_stdev'] = config.set_param(user_config, 'LOC', 'c0_stdev', c0_stdev, 'float')
    loc_params['det_tm_stdev'] = config.set_param(user_config, 'LOC', 'det_tm_stdev', det_tm_stdev, 'float')

    loc_params['az_limit'] = config.set_param(user_config, 'LOC', 'az_limit', az_limit, 'float')

    loc_params['local_temp_dir'] = config.set_param(user_config, 'LOC', 'local_temp_dir', local_temp_dir, 'str')
    loc_params['cpu_cnt'] = config.set_param(user_config, 'LOC', 'cpu_cnt', cpu_cnt, 'int')

    if loc_params['cpu_cnt'] is not None:
        pl = Pool(loc_params['cpu_cnt'])
    else:
        pl = None

    if loc_params['ll_corner'] is not None:
        # set grid from corners
        loc_params['ll_corner'] = np.array([float(val) for val in loc_params['ll_corner'].replace(" ","").split(",")])
        loc_params['ur_corner'] = np.array([float(val) for val in loc_params['ur_corner'].replace(" ","").split(",")])
        tm_lims = (np.datetime64(loc_params['tm_min']), np.datetime64(loc_params['tm_max']))

        if loc_params['alt_lims'] is not None:
            loc_params['alt_lims'] = np.array([float(val) for val in loc_params['alt_lims'].replace(" ","").split(",")])

        loc_params['back_az_width'] = None
        loc_params['range_max'] = None
        loc_params['grid_resol'] = None 
    else:
        # set automatically
        loc_params['alt_resol'] = None
        loc_params['c0_stdev'] = None
        loc_params['az_limit'] = None
        
        tm_lims = None

    if loc_params['atmo_data'] is None:               
        if loc_params['pgm_file'] is not None:
            click.echo("  pgm_file: " + str(loc_params['pgm_file']))
            pgm = infrasound.PathGeometryModel()
            pgm.load(loc_params['pgm_file'])
        else:
            infrasound._load_celerity_model(loc_params["celerity_model"])
            pgm = None
    else:
        loc_params["celerity_model"] = None

    click.echo('\n' + "localization parameters:")
    for key in loc_params.keys():
        # click.echo("  " + key + ": " + str(loc_params[key]))
        if loc_params[key] is not None:
            click.echo("  " + key + ": " + str(loc_params[key]))
    click.echo("")


    # Check if results already exist for this parameter set
    new_param_set = True
    param_index = 0
    for k, result_set in enumerate(ev_data['location']):        
        try:
            new_param_set = False 
            param_index = k
            np.testing.assert_equal(loc_params, result_set['params'])
            break
        except:
            new_param_set = True
            pass

    if new_param_set:
        det_list = [data_io._det_dict_to_likelihood(dict) for dict in ev_data['det_info']]

        click.echo("")
        if loc_params['atmo_data'] is None:           
            result = bisl.run(det_list, bm_width=loc_params['back_az_width'], rng_max=loc_params['range_max'], grid_resol=loc_params['grid_resol'],
                              ll_corner=loc_params['ll_corner'], ur_corner=loc_params['ur_corner'], latlon_resol=loc_params['latlon_resol'],
                              tm_lims=tm_lims, tm_resol=loc_params['tm_resol'], path_geo_model=pgm)
        else:
            if find_spec('infraga'):
                if not os.path.isfile(find_spec('infraga').submodule_search_locations[0] + "/bin/infraga-sph"):
                    click.echo("InfraGA methods not compiled.  Run 'infraga compile' and try again.")
                else:                               
                    with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
                        if loc_params['local_temp_dir'] is not None:
                            if not os.path.isdir(loc_params['local_temp_dir']):
                                os.mkdir(loc_params['local_temp_dir'])
                            tmpdirname = loc_params['local_temp_dir']

                        temp_path = tmpdirname + "/temp"

                        if "*" in loc_params['atmo_data']:
                            if len(os.path.dirname(loc_params['atmo_data'])) > 0:
                                file_path = os.path.dirname(loc_params['atmo_data']) + "/"
                            else:
                                file_path = ""

                            if "/" in loc_params['atmo_data']:
                                dir_files = os.listdir(os.path.dirname(loc_params['atmo_data']))
                            else:
                                dir_files = os.listdir(".")

                            file_list = np.sort([file for file in dir_files if fnmatch.fnmatch(file, os.path.basename(loc_params['atmo_data']))])

                            click.echo('\n' + "Computing localization using atmosphere ensemble:")
                            norms = []
                            for k, file_name in enumerate(file_list):                            
                                print('\t' + str(k + 1) + '/' + str(len(file_list)) + '\t' + file_path + file_name + '\t', end='')
                                temp = tribl.run(det_list, file_path + file_name, temp_path + "-" + str(k), bm_width=loc_params['back_az_width'], rng_max=loc_params['range_max'], grid_resol=loc_params['grid_resol'],
                                                ll_corner=loc_params['ll_corner'], ur_corner=loc_params['ur_corner'], latlon_resol=loc_params['latlon_resol'], tm_lims=tm_lims, tm_resol=loc_params['tm_resol'],
                                                alt_lims=loc_params['alt_lims'], alt_resol=loc_params['alt_resol'], grnd_snd_spd=loc_params['grnd_snd_spd'], c0_stdev=loc_params['c0_stdev'],
                                                det_time_stdev=loc_params['det_tm_stdev'], az_limit=loc_params['az_limit'], verbose=False, show_prog=True, pool=pl) 
                                norms = norms + [temp['norm']]

                            norms = norms / np.sum(norms)

                            print('\n' + "Merging PDFs across the ensemble...")
                            tmp_0 = np.load(temp_path + "-0.pdf.npz")

                            lat_grid, lon_grid, alt_grid, tm_grid = np.meshgrid(tmp_0['lat_vals'], tmp_0['lon_vals'], tmp_0['alt_vals'], tmp_0['tm_vals'], indexing='ij')
                            lat_grid = np.squeeze(lat_grid)
                            lon_grid = np.squeeze(lon_grid)
                            alt_grid = np.squeeze(alt_grid)
                            tm_grid = np.squeeze(tm_grid)

                            pdf = tmp_0['pdf']

                            print('\t1/' + str(len(file_list)) + '\t' + file_path + file_list[0] + '\t' + str(norms[0]))

                            for k, file_name in enumerate(file_list[1:]):                            
                                tmp_k = np.load(temp_path + "-" + str(k + 1) + ".pdf.npz")
                                pdf = pdf + tmp_k['pdf']
                                print('\t' + str(k + 2) + '/' + str(len(file_list)) + '\t' + file_path + file_name + '\t' + str(norms[k + 1]))
                        
                            click.echo('\n' + "Analyzing combined localization PDF...")
                            result = bisl.analyze_pdf(pdf, lat_grid, lon_grid, tm_grid, verbose=True)

                            # Add ensemble information tp location dictionary
                            result['atmo ensemble'] = file_list
                            result['atmo norms'] = norms 

                        else:
                            result = tribl.run(det_list, loc_params['atmo_data'], temp_path, bm_width=loc_params['back_az_width'], rng_max=loc_params['range_max'], grid_resol=loc_params['grid_resol'], 
                                                ll_corner=loc_params['ll_corner'], ur_corner=loc_params['ur_corner'], latlon_resol=loc_params['latlon_resol'], tm_lims=tm_lims, tm_resol=loc_params['tm_resol'], 
                                                alt_lims=loc_params['alt_lims'], alt_resol=loc_params['alt_resol'], grnd_snd_spd=loc_params['grnd_snd_spd'], c0_stdev=loc_params['c0_stdev'],
                                                det_time_stdev=loc_params['det_tm_stdev'], az_limit=loc_params['az_limit'], verbose=True, pool=pl)
            else:
                click.echo('\n' + "Can't run TRIBL methods without infraGA installed for ray tracing")
                return

        # Determine output format for BISL results
        if result['norm'] > 0.0:
            click.echo('\n' + "Localization Summary:")
            click.echo(bisl.summarize(result))

            ev_data['location'] = ev_data['location'] + [{'params' : loc_params, 'result' : result}]
            with gzip.open(event_file, 'wt', encoding='UTF-8') as zipfile:
                json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)
    else:
        click.echo("Localization result already exists in this event file for this parameter set.")
        click.echo('\n' + "Localization Summary:")
        click.echo(bisl.summarize(ev_data['location'][param_index]['result']))


@click.command('characterize', short_help="Characterize the source spectra for an event")
@click.option("--event-file", help="Event JSON file to be analyzed", default=None)
@click.option("--config-file", help="Configuration file", default=None)

@click.option("--det-mask", help="Mask to select detections for analysis", default=None)
@click.option("--loc-index", help="Localization index to use in analysis", default=None)
@click.option("--tlm-label", help="Transmission loss model (TLM) path", default=None)


@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['YIELD']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['YIELD']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--yld-min", help="Minimum yield (default: " + config.defaults['YIELD']['yld_min'] + " [tons eq. TNT])", default=1.0, type=float)
@click.option("--yld-max", help="Maximum yield (default: " + config.defaults['YIELD']['yld_max'] + " [tons eq. TNT])", default=1000.0, type=float)
@click.option("--ref-rng", help="Reference range for blastwave model (default " + config.defaults['YIELD']['ref_rng'] + " km)", default=1.0, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['YIELD']['resolution'] + ")", default=None, type=int)
@click.option("--amb-press", help="Ambient pressure (default: " + config.defaults['YIELD']['amb_press'] + " [Pa])", default=None, type=float)
@click.option("--amb-temp", help="Ambient temperature (default: " + config.defaults['YIELD']['amb_temp'] + " [K])", default=None, type=float)
@click.option("--grnd-burst", help="Ground burst assumption (default: " + config.defaults['YIELD']['grnd_burst'] + " [Hz])", default=None, type=bool)
@click.option("--exp-type", help="Explosion type ('chemical' or 'nuclear')", default=None)
def characterize(event_file, config_file, det_mask, loc_index, tlm_label, freq_min, freq_max, yld_min, yld_max, ref_rng, resolution, amb_press, amb_temp, grnd_burst, exp_type):
    '''
    Run Bayesian Infrasonic Source Localization (BISL) methods to estimate the source location and origin time for an event

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy characterize --event-file GJI_example-ev0
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##    Characterization Analysis    ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None    

    event_file = config.set_param(user_config, 'DATA IO', 'event_file', event_file, 'string')
    det_mask = config.set_param(user_config, 'YIELD', 'det_mask', det_mask, 'string')
    loc_index = config.set_param(user_config, 'YIELD', 'loc_index', loc_index, 'str')

    ev_data = data_io._load_dets_json(event_file)[0]

    if det_mask is None:
        det_mask = np.ones(len(ev_data['det_info']))
    else:
        det_mask = [int(val) for val in det_mask.strip(' ()[]').split(',')]
            
    if loc_index is not None:
        src_loc = [loc_lats[loc_index], loc_lons[loc_index]]
    else:
        loc_std = [loc['result']['NS_stdev'] * loc['result']['NS_stdev'] for loc in ev_data['location']]
        loc_lats = [loc['result']['lat_mean'] for loc in ev_data['location']]
        loc_lons = [loc['result']['lon_mean'] for loc in ev_data['location']]

        loc_index = np.argmin(loc_std)
        src_loc = [loc_lats[loc_index], loc_lons[loc_index]]

    click.echo('\n' + "Data summary:")
    click.echo("  event_file: " + str(event_file))

    # Set analysis parameter dictionary
    char_params = {}
    char_params['det_mask'] = det_mask
    char_params['loc_index'] = loc_index

    char_params['tlm_label'] = config.set_param(user_config, 'YIELD', 'tlm_label', tlm_label, 'str')

    char_params['freq_min'] = config.set_param(user_config, 'YIELD', 'freq_min', freq_min, 'float')
    char_params['freq_max'] = config.set_param(user_config, 'YIELD', 'freq_max', freq_max, 'float')
    char_params['yld_min'] = config.set_param(user_config, 'YIELD', 'yld_min', yld_min, 'float')
    char_params['yld_max'] = config.set_param(user_config, 'YIELD', 'yld_max', yld_max, 'float')
    char_params['resolution'] = config.set_param(user_config, 'YIELD', 'resolution', resolution, 'int')

    char_params['ref_rng'] = config.set_param(user_config, 'YIELD', 'ref_rng', ref_rng, 'float')
    char_params['amb_press'] = config.set_param(user_config, 'YIELD', 'amb_press', amb_press, 'float')
    char_params['amb_temp'] = config.set_param(user_config, 'YIELD', 'amb_temp', amb_temp, 'float')
    char_params['grnd_burst'] = config.set_param(user_config, 'YIELD', 'grnd_burst', grnd_burst, 'bool')
    char_params['exp_type'] = config.set_param(user_config, 'YIELD', 'exp_type', exp_type, 'str')

    char_params['det_mask'] = str(det_mask)

    click.echo('\n' + "characterization parameters:")
    for key in char_params.keys():
        if char_params[key] is not None:
            click.echo("  " + key + ": " + str(char_params[key]))
    click.echo("")

    # ######################### #
    #      Load Detections      #
    # ######################### #
    det_info = [dict for k, dict in enumerate(ev_data['det_info']) if bool(det_mask[k])]

    det_list = [data_io._det_dict_to_likelihood(dict) for dict in det_info]
    det_specs = [spye.extract_json_spectra(det['spec']) for det in det_info]

    click.echo("=" * 17 + '\n' + "Detection Summary" + '\n' + "=" * 17 + '\n')
    for k, det in enumerate(det_info):
        freq = det_specs[k][0]
        mask = (det_specs[k][1] / det_specs[k][2]) > 2.0

        click.echo(det['wvfrm_info'][0][0]['trace id'])
        click.echo("  location: " + str(det['wvfrm_info'][0][0]['latitude']) + ", " + str(det['wvfrm_info'][0][0]['longitude']))
        click.echo("  detection time: " + det['peak f-stat time'])
        click.echo("  back azimuth [deg]: " + str(np.round(det["back az"], 2)))
        click.echo("  tface velocity [m/s]: " + str(np.round(det["tr vel"], 2)))
        click.echo("  f-stat: " + str(np.round(det["f-stat"], 2)))
        click.echo("  high snr band: " + str(np.round(freq[mask][0], 2)) + " - " + str(np.round(freq[mask][-1], 2)) + ' Hz\n')

    # drop residual spectra and scale to dB
    det_specs = [np.array([spec[0], 10.0 * np.log10(spec[1])]) for spec in det_specs]

    # Check if a result for this parameter set already exists.
    new_param_set = True
    param_index = 0
    for k, result_set in enumerate(ev_data['characterization']):        
        try:
            new_param_set = False 
            param_index = k
            np.testing.assert_equal(char_params, result_set['params'])
            break
        except:
            new_param_set = True
            pass

    if new_param_set:
        # ######################### #
        #     Load TLoss Models     #
        # ######################### #
        click.echo("Loading transmission loss statistics...")
        tlm_dir = os.path.dirname(char_params['tlm_label'])
        tlm_pattern = char_params['tlm_label'].split("/")[-1]
        tlm_files = [file_name for file_name in np.sort(os.listdir(tlm_dir)) if fnmatch.fnmatch(file_name, tlm_pattern + "*")]

        models = [0] * 2
        models[0] = [float(file_name.split("Hz")[0][len(tlm_pattern):]) for file_name in tlm_files]
        models[1] = [0] * len(tlm_files)
        for n in range(len(tlm_files)):
            models[1][n] = infrasound.TLossModel()
            models[1][n].load(tlm_dir + "/" + tlm_files[n])
        

        # ######################## #
        #         Run Yield        #
        #    Estimation Methods    #
        # ######################## #
        spye_result = spye.run(det_list, det_specs, src_loc, np.array([char_params['freq_min'], char_params['freq_max']]), models, 
                            yld_rng=np.array([char_params['yld_min'] * 1.0e3, char_params['yld_max'] * 1.0e3]),
                            ref_src_rng=char_params['ref_rng'], resol=char_params['resolution'], grnd_brst= char_params['grnd_burst'],
                            p_amb= char_params['amb_press'], T_amb= char_params['amb_temp'], exp_type= char_params['exp_type'])


        ev_data['characterization'] = ev_data['characterization'] + [{'params' : char_params, 'result' : spye_result}]
        with gzip.open(event_file, 'wt', encoding='UTF-8') as zipfile:
            json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

        click.echo('\n' + 'Results Summary (tons eq. TNT):')
        click.echo('\t' + "Maximum a Posteriori Yield: " + str(spye_result['yld_vals'][np.argmax(spye_result['yld_pdf'])]))
        click.echo('\t' + "68% Confidence Bounds: " + str(spye_result['conf_bnds'][0]))
        click.echo('\t' + "95% Confidence Bounds: " + str(spye_result['conf_bnds'][1]))
        click.echo('')

    else:
        spye_result = ev_data['characterization'][param_index]['result']

        click.echo("Characterization result already exists in this event file for this parameter set:")
        click.echo('\t' + "Maximum a Posteriori Yield: " + str(spye_result['yld_vals'][np.argmax(spye_result['yld_pdf'])]))
        click.echo('\t' + "68% Confidence Bounds: " + str(spye_result['conf_bnds'][0]))
        click.echo('\t' + "95% Confidence Bounds: " + str(spye_result['conf_bnds'][1]))
        click.echo('')
