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

from ..utils import config
from ..utils import data_io

from ..association import hjl

from ..location import bisl, tribl
from ..propagation import infrasound


@click.command('build', short_help="Associate detections into events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--detect-files", help="Detection path and pattern", default=None)
@click.option("--event-label", help="Path for event info output", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
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
def build(config_file, detect_files, event_label, starttime, endtime, back_az_width, range_max, resolution, distance_matrix_max, cluster_linkage, 
                cluster_threshold, trimming_threshold, event_population_min, event_station_min, cpu_cnt):
    '''
    Run association analysis to identify events in a detection set

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy event build --detect-files 'data/Blom_etal2020_GJI/*' --event-label GJI_example --cpu-cnt 4
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

    if cpu_cnt is not None:
        pl = Pool(cpu_cnt)
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

    det_list = [data_io._det_dict_to_likelihood(dict) for dict in det_dicts]

    events, event_qls = hjl.id_events(det_list, assoc_params['cluster_threshold'], starttime=assoc_params['starttime'], endtime=assoc_params['endtime'], dist_max=assoc_params['distance_matrix_max'], 
                                    bm_width=assoc_params['back_az_width'], rng_max=assoc_params['range_max'], rad_min=100.0, rad_max=(assoc_params['range_max'] / 4.0), 
                                    resol=assoc_params['resolution'], linkage_method=assoc_params['cluster_linkage'], trimming_thresh=assoc_params['trimming_threshold'], 
                                    cluster_det_population=assoc_params['event_population_min'], cluster_array_population=assoc_params['event_station_min'], pool=pl)

    click.echo("Identified " + str(len(events)) + " events." + '\n')
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

@click.option("--celerity-model", help="Use included celerity model (default: '" + config.defaults['LOC']['celerity_model'], default=None, hidden=True)
@click.option("--rcel-wts", help="Custom reciprocal celerity model weights", default=None, hidden=True)
@click.option("--rcel-mns", help="Custom reciprocal celerity model means", default=None, hidden=True)
@click.option("--rcel-sds", help="Custom reciprocal celerity model standard deviations", default=None, hidden=True)
@click.option("--pgm-file", help="Path geometry model (PGM) file (optional)", default=None)

@click.option("--atmo-data", help="Atmosphere data if using TRIBL", default=None)
@click.option("--alt-lims", help="Altitude limits if using TRIBL", default=None)
@click.option("--alt-resol", help="Altitude resolution if using TRIBL", default=None)
@click.option("--grnd-snd-spd", help="Sound speed at the ground", default=None)
@click.option("--c0-stdev", help="Sound speed uncertainty if using TRIBL", default=None)
@click.option("--det-tm-stdev", help="Detection time uncertainty", default=None)
@click.option("--local-temp-dir", help="Local temporary directory if using TRIBL", default=None)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)

def localize(event_file, config_file, back_az_width, range_max, grid_resol, ll_corner, ur_corner, latlon_resol, tm_min, tm_max, tm_resol, celerity_model, rcel_wts, rcel_mns, rcel_sds, pgm_file, atmo_data, alt_lims, alt_resol, grnd_snd_spd, c0_stdev, det_tm_stdev, local_temp_dir, cpu_cnt):
    '''
    Run Bayesian Infrasonic Source Localization (BISL) methods to estimate the source location and origin time for an event

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy localize --event-file GJI_example-ev0
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
    loc_params['rcel_wts'] = config.set_param(user_config, 'LOC', 'rcel_wts', rcel_wts, 'str')
    loc_params['rcel_mns'] = config.set_param(user_config, 'LOC', 'rcel_mns', rcel_mns, 'str')
    loc_params['rcel_sds'] = config.set_param(user_config, 'LOC', 'rcel_sds', rcel_sds, 'str')

    loc_params['pgm_file'] = config.set_param(user_config, 'LOC', 'pgm_file', pgm_file, 'str')

    loc_params['atmo_data'] = config.set_param(user_config, 'LOC', 'atmo_data', atmo_data, 'str')
    loc_params['alt_lims'] = config.set_param(user_config, 'LOC', 'alt_lims', alt_lims, 'str')
    loc_params['alt_resol'] = config.set_param(user_config, 'LOC', 'alt_resol', alt_resol, 'float')
    loc_params['grnd_snd_spd'] = config.set_param(user_config, 'LOC', 'grnd_snd_spd', grnd_snd_spd, 'float')
    loc_params['c0_stdev'] = config.set_param(user_config, 'LOC', 'c0_stdev', c0_stdev, 'float')
    loc_params['det_tm_stdev'] = config.set_param(user_config, 'LOC', 'det_tm_stdev', det_tm_stdev, 'float')
    loc_params['local_temp_dir'] = config.set_param(user_config, 'LOC', 'local_temp_dir', local_temp_dir, 'str')
    loc_params['cpu_cnt'] = config.set_param(user_config, 'ASSOC', 'cpu_cnt', cpu_cnt, 'int')


    if loc_params['ll_corner'] is not None:
        # set grid from corners
        loc_params['ll_corner'] = np.array([float(val) for val in ll_corner.replace(" ","").split(",")])
        loc_params['ur_corner'] = np.array([float(val) for val in ur_corner.replace(" ","").split(",")])
        tm_lims = (np.datetime64(loc_params['tm_min']), np.datetime64(loc_params['tm_max']))

        loc_params['back_az_width'] = None
        loc_params['range_max'] = None
    else:
        # set automatically
        loc_params['alt_resol'] = None
        loc_params['c0_stdev'] = None
        tm_lims = None


    click.echo('\n' + "localization parameters:")
    for key in loc_params.keys():
        click.echo("  " + key + ": " + str(loc_params[key]))
        #if loc_params[key] is not None:
        #    click.echo("  " + key + ": " + str(loc_params[key]))
    click.echo("")

    infrasound.set_celerity_model(loc_params["celerity_model"], rcel_wts=loc_params['rcel_wts'], rcel_mns=loc_params['rcel_mns'], rcel_sds=loc_params['rcel_sds'])

    if loc_params['pgm_file'] is not None:
        click.echo("  pgm_file: " + str(loc_params['pgm_file']))
        pgm = infrasound.PathGeometryModel()
        pgm.load(loc_params['pgm_file'])
    else:
        pgm = None

    ev_data = data_io._load_dets_json(event_file)[0]

    # Check if results already exist for this parameter set
    new_param_set = True
    param_index = 0
    for k, result_set in enumerate(ev_data['location']):
        if loc_params == result_set['params']:
            new_param_set = False
            param_index = k

    if new_param_set:
        det_list = [data_io._det_dict_to_likelihood(dict) for dict in ev_data['det_info']]

        click.echo("")
        if atmo_data is None:           
            result = bisl.run(det_list, bm_width=loc_params['back_az_width'], rng_max=loc_params['range_max'], grid_resol=loc_params['grid_resol'],
                              ll_corner=loc_params['ll_corner'], ur_corner=loc_params['ur_corner'], latlon_resol=loc_params['latlon_resol'],
                              tm_lims=tm_lims, tm_resol=loc_params['tm_resol'], path_geo_model=pgm)
        else:
            if find_spec('infraga'):
                if not os.path.isfile(find_spec('infraga').submodule_search_locations[0] + "/bin/infraga-sph"):
                    click.echo("InfraGA methods not compiled.  Run 'infraga compile' and try again.")
                else:                               
                    with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
                        if local_temp_dir is not None:
                            if not os.path.isdir(local_temp_dir):
                                os.mkdir(local_temp_dir)
                            tmpdirname = local_temp_dir

                        temp_path = tmpdirname + "/temp"

                        if "*" in atmo_data:
                            if len(os.path.dirname(atmo_data)) > 0:
                                file_path = os.path.dirname(atmo_data) + "/"
                            else:
                                file_path = ""

                            if "/" in atmo_data:
                                dir_files = os.listdir(os.path.dirname(atmo_data))
                            else:
                                dir_files = os.listdir(".")

                            file_list = np.sort([file for file in dir_files if fnmatch.fnmatch(file, os.path.basename(atmo_data))])

                            click.echo('\n' + "Computing localization using atmosphere ensemble:")
                            norms = []
                            for k, file_name in enumerate(file_list):                            
                                print('\t' + str(k + 1) + '/' + str(len(file_list)) + '\t' + file_path + file_name + '\t', end='')
                                temp = tribl.run(det_list, file_path + file_name, temp_path + "-" + str(k), bm_width=back_az_width, rng_max=range_max, grid_resol=grid_resol, ll_corner=ll_corner, ur_corner=ur_corner,
                                                latlon_resol=latlon_resol, tm_lims=tm_lims, tm_resol=tm_resol, alt_lims=alt_lims, alt_resol=alt_resol, grnd_snd_spd=grnd_snd_spd, c0_stdev=c0_stdev, det_time_stdev=det_tm_stdev, verbose=False, show_prog=True, pool=pl) 
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

                        else:              
                            result = tribl.run(det_list, atmo_data, temp_path, bm_width=back_az_width, rng_max=range_max, grid_resol=grid_resol, ll_corner=ll_corner, ur_corner=ur_corner,
                                                latlon_resol=latlon_resol, tm_lims=tm_lims, tm_resol=tm_resol, alt_lims=alt_lims, alt_resol=alt_resol, grnd_snd_spd=grnd_snd_spd, c0_stdev=c0_stdev, det_time_stdev=det_tm_stdev, verbose=True, pool=pl)
            else:
                click.echo('\n' + "Can't run TRIBL methods without infraGA installed for ray tracing")
                return

        # Determine output format for BISL results
        click.echo('\n' + "Localization Summary:")
        click.echo(bisl.summarize(result))

        ev_data['location'] = ev_data['location'] + [{'params' : loc_params, 'result' : result}]

        with gzip.open("test.ev.json.gz", 'wt', encoding='UTF-8') as zipfile:
            json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)
    else:
        click.echo("Localization result already exists in this event file for this parameter set.")
        click.echo('\n' + "BISL Summary:")
        click.echo(bisl.summarize(ev_data['location'][param_index]['result']))

