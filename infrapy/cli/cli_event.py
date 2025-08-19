#!/usr/bin/env python
import sys
import os 
import click
import json
import gzip
import warnings

import configparser as cnfg
import numpy as np

from multiprocessing import Pool

from datetime import datetime
from obspy import UTCDateTime

from ..utils import config
from ..utils import data_io

from ..association import hjl

from ..location import bisl
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
