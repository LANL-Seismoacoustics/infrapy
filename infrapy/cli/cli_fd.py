#!/usr/bin/env python
import sys

import click
import configparser as cnfg
import numpy as np

from multiprocessing import Pool
from obspy.core import read as obspy_read

from ..utils import config
from ..detection import beamforming_new as fkd

@click.command('run_fd', short_help="Identify detections from beamforming results")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-fk", help="Local beamforming (fk) data files", default=None)
@click.option("--local-dets", help="Local detection data files", default=None)
@click.option("--window_len", help="Adaptive window length (default: " + config.defaults['FD']['window_len'] + " [s])", default=None, type=float)
@click.option("--p-value", help="Detection p-value (default: " + config.defaults['FD']['p_value'] + ")", default=None, type=float)
@click.option("--min-duration", help="Minimum detection duration (default: " + config.defaults['FD']['min_duration'] + " [s])", default=None, type=float)
@click.option("--back-az-width", help="Maximum azimuth scatter (default: " + config.defaults['FD']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--fixed-thresh", help="Fixed f-stat threshold (default: None)", default=None, type=float)
@click.option("--return-thresh", help="Return threshold (default: " + config.defaults['FD']['return_thresh'] + ")", default=None, type=bool)
def run_fd(config_file, local_fk, local_dets, window_len, p_value, min_duration, back_az_width, fixed_thresh, return_thresh):
    '''
    Identify detections

    More detailed description...
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##     Detection (fd) Analysis     ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Data IO parameters
    # use local ingestion for initial testing
    local_fk = config.set_param(user_config, 'DETECTION IO', 'local_fk', local_fk, 'string')
    local_dets = config.set_param(user_config, 'DETECTION IO', 'local_dets', local_dets, 'string')

    click.echo('\n' + "Data parameters:")
    click.echo("  local_fk: " + local_fk)
    click.echo("  local_dets: " + local_dets)

    # Algorithm parameters
    window_len = config.set_param(user_config, 'FD', 'window_len', window_len, 'float')
    p_value = config.set_param(user_config, 'FD', 'p_value', p_value, 'float')
    min_duration = config.set_param(user_config, 'FD', 'min_duration', min_duration, 'float')
    back_az_width = config.set_param(user_config, 'FD', 'back_az_width', back_az_width, 'float')
    fixed_thresh = config.set_param(user_config, 'FD', 'fixed_thresh', fixed_thresh, 'float')
    return_thresh = config.set_param(user_config, 'FD', 'return_thresh', return_thresh, 'bool')

    click.echo('\n' + "Algorithm arameters:")
    click.echo("  window_len: " + str(window_len))
    click.echo("  p_value: " + str(p_value))
    click.echo("  min_duration: " + str(min_duration))
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  fixed_thresh: " + str(fixed_thresh))
    click.echo("  return_thresh: " + str(return_thresh))

    print('\n' + "Run fd...")

    if local_fk is not None:
        beam_times = np.load(local_fk + ".fk_times.npy")
        beam_peaks = np.load(local_fk + ".fk_peaks.npy")
        beam_meta = open(local_fk + ".fk_meta.txt")
    else:
        print("Non-local data not yet set up...")
        return 0

    TB_prod = (5.0 - 0.5) * 10.0
    channel_cnt = 4
    min_seq = int(max(2, min_duration / window_len))

    dets, thresh_vals = fkd.run_fd(beam_times, beam_peaks, window_len, TB_prod, channel_cnt, p_value, min_seq, back_az_width, fixed_thresh, True)

    for det in dets:
        print(det)

    if return_thresh:
        print(thresh_vals)