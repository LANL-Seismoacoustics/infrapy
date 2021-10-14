#!/usr/bin/env python
import sys

import click
import configparser as cnfg
import numpy as np

from multiprocessing import Pool
from obspy.core import read as obspy_read

from ..utils import config
from ..association import hjl

@click.command('run_assoc', short_help="Associate detections into events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--back_az_width", help="Width of beam projection (default: " + config.defaults['ASSOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Maximum source-receiver range (default: " + config.defaults['ASSOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['ASSOC']['resolution'] + ")", default=None, type=int)
@click.option("--distance-matrix-max", help="Distance matrix maximum (default: " + config.defaults['ASSOC']['distance_matrix_max'] + ")", default=None, type=float)
@click.option("--cluster-threshold", help="Cluster linkage threshold (default: " + config.defaults['ASSOC']['cluster_threshold'] + ")", default=None, type=float)
@click.option("--trimming-threshold", help="Mishapen cluster threshold (default: " + config.defaults['ASSOC']['trimming_threshold'] + ")", default=None, type=float)
@click.option("--multithread", help="Use multithreading (default: " + config.defaults['FK']['multithread'] + ")", default=None, type=bool)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def run_assoc(config_file, back_az_width, range_max, resolution, distance_matrix_max, cluster_threshold, trimming_threshold, multithread, cpu_cnt):
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
    click.echo('\n' + "Data summary:")
    click.echo("  NOT SET UP YET...")

    # Algorithm parameters
    back_az_width = config.set_param(user_config, 'ASSOC', 'back_az_width', back_az_width, 'float')
    range_max = config.set_param(user_config, 'ASSOC', 'range_max', range_max, 'float')
    resolution = config.set_param(user_config, 'ASSOC', 'resolution', resolution, 'int')
    distance_matrix_max = config.set_param(user_config, 'ASSOC', 'distance_matrix_max', distance_matrix_max, 'float')
    cluster_threshold = config.set_param(user_config, 'ASSOC', 'cluster_threshold', cluster_threshold, 'float')
    trimming_threshold = config.set_param(user_config, 'ASSOC', 'trimming_threshold', trimming_threshold, 'float')
    multithread = config.set_param(user_config, 'ASSOC', 'multithread', multithread, 'bool')
    cpu_cnt = config.set_param(user_config, 'ASSOC', 'cpu_cnt', cpu_cnt, 'int')

    click.echo('\n' + "Parameter summary:")
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  range_max: " + str(range_max))
    click.echo("  resolution: " + str(resolution))
    click.echo("  distance_matrix_max: " + str(distance_matrix_max))
    click.echo("  cluster_threshold: " + str(cluster_threshold))
    click.echo("  trimming_threshold: " + str(trimming_threshold))

    click.echo("  multithread: " + str(multithread))
    if multithread:
        click.echo("  cpu_cnt: " + str(cpu_cnt))

    print('\n' + "Run association analysis...")

    print('\n' + "Write results to specified output..." + '\n')