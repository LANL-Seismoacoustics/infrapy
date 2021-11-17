#!/usr/bin/env python
import sys

import click
import warnings

import configparser as cnfg
import numpy as np

from multiprocessing import Pool

from datetime import datetime
from obspy import UTCDateTime

from ..utils import config
from ..utils import data_io
from ..association import hjl


@click.command('run_assoc', short_help="Associate detections into events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-dets-in", help="Detection path and pattern", default=None)
@click.option("--local-events-out", help="Path for event info output", default=None)
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
@click.option("--multithread", help="Use multithreading (default: " + config.defaults['ASSOC']['multithread'] + ")", default=None, type=bool)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def run_assoc(config_file, local_dets_in, local_events_out, starttime, endtime, back_az_width, range_max, resolution, distance_matrix_max, cluster_linkage, 
                cluster_threshold, trimming_threshold, event_population_min, event_station_min, multithread, cpu_cnt):
    '''
    Associate detections

    More detailed description...
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##       Association Analysis      ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo("Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None


    # Data IO parameters
    local_dets_in = config.set_param(user_config, 'DETECTION IO', 'local_dets_in', local_dets_in, 'string')
    local_events_out = config.set_param(user_config, 'DETECTION IO', 'local_events_out', local_events_out, 'string')
    starttime = config.set_param(user_config, 'DETECTION IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'DETECTION IO', 'endtime', endtime, 'string')

    # Data IO parameters
    click.echo('\n' + "Data summary:")
    click.echo("  local_dets_in: " + str(local_dets_in))
    click.echo("  local_events_out: " + str(local_events_out))
    click.echo("  starttime: " + str(starttime))
    click.echo("  endtime: " + str(endtime))

    if local_dets_in is None or local_events_out is None:
        msg = "Association analysis requires detection input (--local-det-in) and output path (--local-events-out)"
        warnings.warn(msg)
        return 0

    # Algorithm parameters
    back_az_width = config.set_param(user_config, 'ASSOC', 'back_az_width', back_az_width, 'float')
    range_max = config.set_param(user_config, 'ASSOC', 'range_max', range_max, 'float')
    resolution = config.set_param(user_config, 'ASSOC', 'resolution', resolution, 'int')
    distance_matrix_max = config.set_param(user_config, 'ASSOC', 'distance_matrix_max', distance_matrix_max, 'float')
    cluster_linkage = config.set_param(user_config, 'ASSOC', 'cluster_linkage', cluster_linkage, 'string')
    cluster_threshold = config.set_param(user_config, 'ASSOC', 'cluster_threshold', cluster_threshold, 'float')
    trimming_threshold = config.set_param(user_config, 'ASSOC', 'trimming_threshold', trimming_threshold, 'float')
    event_population_min = config.set_param(user_config, 'ASSOC', 'event_population_min', event_population_min, 'float')
    event_station_min = config.set_param(user_config, 'ASSOC', 'event_station_min', event_station_min, 'float')
    multithread = config.set_param(user_config, 'ASSOC', 'multithread', multithread, 'bool')
    cpu_cnt = config.set_param(user_config, 'ASSOC', 'cpu_cnt', cpu_cnt, 'int')

    click.echo('\n' + "Parameter summary:")
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  range_max: " + str(range_max))
    click.echo("  resolution: " + str(resolution))
    click.echo("  distance_matrix_max: " + str(distance_matrix_max))
    click.echo("  cluster_linkage: " + str(cluster_linkage))
    click.echo("  cluster_threshold: " + str(cluster_threshold))
    click.echo("  trimming_threshold: " + str(trimming_threshold))
    if multithread or cpu_cnt is not None:
        click.echo("  cpu_cnt: " + str(cpu_cnt))
        pl = Pool(cpu_cnt)
    else:
        pl = None

    det_list = data_io.set_det_list(local_dets_in, merge=True)
    events, event_qls = hjl.id_events(det_list, cluster_threshold, starttime=starttime, endtime=endtime, dist_max=distance_matrix_max, 
                                    bm_width=back_az_width, rng_max=range_max, rad_min=100.0, rad_max=(range_max / 4.0), 
                                    resol=resolution, linkage_method=cluster_linkage, trimming_thresh=trimming_threshold, 
                                    cluster_det_population=event_population_min, cluster_array_population=event_station_min, pool=pl)
    data_io.write_events(events, event_qls, det_list, local_events_out)    
    if pl is not None:
        pl.terminate()
        pl.close()
