#!/usr/bin/env python
import sys

import os
import imp 
import click

import configparser as cnfg
import numpy as np

from obspy.core import read as obspy_read

from multiprocessing import Pool


from ..detection import beamforming_new as fkd
from ..association import hjl
from ..location import bisl

# Set up default configuation
default_config = cnfg.ConfigParser()
default_config.read(imp.find_module('infrapy')[1] + '/resources/default.config')

def set_param(user_config, section, param, cli_val, format='float'):   
    if cli_val is not None:
        return cli_val
    elif user_config is not None:
        try:
            if user_config[section][param] == "None":
                return None
            else:
                if format == 'float':
                    return float(user_config[section][param])
                elif format == 'int':
                    return int(user_config[section][param])
                elif format == 'bool':
                    return user_config[section].getboolean(param)
                else:
                    return user_config[section][param]
        except:
            try:
                if default_config[section][param] == "None":
                    return None
                else:
                    if format == 'float':
                        return float(default_config[section][param])
                    elif format == 'int':
                        return int(default_config[section][param])
                    elif format == 'bool':
                        return default_config[section].getboolean(param)
                    else:
                        return user_config[section][param]
            except:
                return None
    else:
        try:
            return default_config[section][param]
        except:
            return None


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main(args=None):
    '''
    infrapy - Python-based Infrasound Data Analysis Toolkit

    More detailed description...
    '''
    pass


@main.command('run_fk', short_help="Run beamforming methods on waveform data")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-url", help="Database URL for waveform data files", default=None)
@click.option("--db-site", help="Database site table for waveform data files", default=None)
@click.option("--db-wfdisc", help="Database wfdisc table for waveform data files", default=None)
@click.option("--db-origin", help="Database origin table for waveform data files", default=None)

@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)

@click.option("--local-fk-out", help="Local beamforming (fk) data files", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + default_config['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + default_config['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--back-az-min", help="Minimum back azimuth (default: " + default_config['FK']['back_az_min'] + " [deg])", default=None, type=float)
@click.option("--back-az-max", help="Maximum back azimuth (default: " + default_config['FK']['back_az_max'] + " [deg])", default=None, type=float)
@click.option("--back-az-step", help="Back azimuth resolution (default: " + default_config['FK']['back_az_step'] + " [deg])", default=None, type=float)
@click.option("--trace-vel-min", help="Minimum trace velocity (default: " + default_config['FK']['trace_vel_min'] + " [m/s])", default=None, type=float)
@click.option("--trace-vel-max", help="Maximum trace velocity (default: " + default_config['FK']['trace_vel_max'] + " [m/s])", default=None, type=float)
@click.option("--trace-vel-step", help="Trace velocity resolution (default: " + default_config['FK']['trace_vel_step'] + " [m/s])", default=None, type=float)
@click.option("--method", help="Beamforming method (default: " + default_config['FK']['method'] + ")", default=None)
@click.option("--signal_start", help="Start of analysis window (rel. starttime [s])", default=None, type=float)
@click.option("--signal_end", help="End of analysis window (rel. starttime [s])", default=None, type=float)
@click.option("--noise_start", help="Start of noise sample (see GLS algorithm notes)", default=None, type=float)
@click.option("--noise_end", help="End of noise sample (see GLS algorithm notes)", default=None, type=float)
@click.option("--window_len", help="Analysis window length (default: " + default_config['FK']['window_len'] + " [s])", default=None, type=float)
@click.option("--sub-window_len", help="Analysis sub-window length (default: None [s])", default=None, type=float)
@click.option("--window_step", help="Step between analysis windows (default: " + default_config['FK']['window_step'] + " [s])", default=None, type=float)
@click.option("--multithread", help="Use multithreading (default: " + default_config['FK']['multithread'] + ")", default=None, type=bool)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def run_fk(config_file, local_wvfrms, fdsn, db_url, db_site, db_wfdisc, db_origin, network, station, location, channel, starttime, endtime,
    local_fk_out, freq_min, freq_max, back_az_min, back_az_max, back_az_step, trace_vel_min, trace_vel_max, trace_vel_step, method, 
    signal_start, signal_end, noise_start, noise_end, window_len, sub_window_len, window_step, multithread, cpu_cnt):
    '''
    Run fk beamforming methods

    More detailed description...
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##    Beamforming (fk) Analysis    ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Database and data IO parameters   
    db_url = set_param(user_config, 'database', 'url', db_url, 'string')
    db_site = set_param(user_config, 'database', 'site', db_site, 'string')
    db_wfdisc = set_param(user_config, 'database', 'wfdisc', db_wfdisc, 'string')
    db_origin = set_param(user_config, 'database', 'origin', db_origin, 'string')

    # Waveform IO parameters
    local_wvfrms = set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    fdsn = set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')
   
    network = set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       
    starttime = set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    # Result IO
    local_fk_out = set_param(user_config, 'DETECTION IO', 'local_fk_out', local_fk_out, 'string')


    click.echo('\n' + "Data parameters:")
    if local_wvfrms is not None:
        click.echo("  local_wvfrms: " + str(local_wvfrms))
    elif fdsn is not None:
        click.echo("  fdsn: " + str(fdsn))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    elif db_url is not None:
        click.echo("  db_url: " + str(db_url))
        click.echo("  db_site: " + str(db_site))
        click.echo("  db_wfdisc: " + str(db_wfdisc))
        click.echo("  db_origin: " + str(db_origin))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    else:
        click.echo("Invalid data parameters.  Config file requires 1 of:")
        click.echo("  local_wvfrms")
        click.echo("  fdsn")
        click.echo("  db_url (and other database info)")
        
    click.echo("  local_fk_out: " + str(local_fk_out))

    # Algorithm parameters
    freq_min = set_param(user_config, 'FK', 'freq_min', freq_min, 'float')
    freq_max = set_param(user_config, 'FK', 'freq_max', freq_max, 'float')
    back_az_min = set_param(user_config, 'FK', 'back_az_min', back_az_min, 'float')
    back_az_max = set_param(user_config, 'FK', 'back_az_max', back_az_max, 'float')
    back_az_step = set_param(user_config, 'FK', 'back_az_step', back_az_step, 'float')
    trace_vel_min = set_param(user_config, 'FK', 'trace_vel_min', trace_vel_min, 'float')
    trace_vel_max = set_param(user_config, 'FK', 'trace_vel_max', trace_vel_max, 'float')
    trace_vel_step = set_param(user_config, 'FK', 'trace_vel_step', trace_vel_step, 'float')
    method = set_param(user_config, 'FK', 'method', method, 'string')
    signal_start = set_param(user_config, 'FK', 'signal_start', signal_start, 'float')
    signal_end = set_param(user_config, 'FK', 'signal_end', signal_end, 'float')
    noise_start = set_param(user_config, 'FK', 'noise_start', noise_start, 'float')
    noise_end = set_param(user_config, 'FK', 'noise_end', noise_end, 'float')
    window_len = set_param(user_config, 'FK', 'window_len', window_len, 'float')
    sub_window_len = set_param(user_config, 'FK', 'sub_window_len', sub_window_len, 'float')
    window_step = set_param(user_config, 'FK', 'window_step', window_step, 'float')
    multithread = set_param(user_config, 'FK', 'multithread', multithread, 'bool')
    cpu_cnt = set_param(user_config, 'FK', 'cpu_cnt', cpu_cnt, 'int')

    click.echo('\n' + "Algorithm parameters:")
    click.echo("  freq_min: " + str(freq_min))
    click.echo("  freq_max: " + str(freq_max))
    click.echo("  back_az_min: " + str(back_az_min))
    click.echo("  back_az_max: " + str(back_az_max))
    click.echo("  back_az_step: " + str(back_az_step))
    click.echo("  trace_vel_min: " + str(trace_vel_min))
    click.echo("  trace_vel_max: " + str(trace_vel_max))
    click.echo("  trace_vel_step: " + str(trace_vel_step))
    click.echo("  method: " + str(method))
    click.echo("  signal_start: " + str(signal_start))
    click.echo("  signal_end: " + str(signal_end))
    if method == "GLS":
        click.echo("  noise_start: " + str(noise_start))
        click.echo("  noise_end: " + str(noise_end))
    click.echo("  window_len: " + str(window_len))
    click.echo("  sub_window_len: " + str(sub_window_len))
    click.echo("  window_step: " + str(window_step))
    click.echo("  multithread: " + str(multithread))
    if multithread or cpu_cnt is not None:
        click.echo("  cpu_cnt: " + str(cpu_cnt))
        pl = Pool(cpu_cnt)
    else:
        pl = None

    # Check data option and populate obspy Stream
    # Note: this might become a separate function that returns stream and latlon info
    if local_wvfrms is not None:
        print('\n' + "Loading local data files...")
        stream = obspy_read(local_wvfrms)
    elif fdsn is not None:
        print('\n' + "FDSN methods not set up yet...")

        # client = Client(service)
        # stream = client.get_waveforms(network, station, location, channel, startTime, endTime)
        # inventory = client.get_stations(network=network, station=station)
        
        return 0
    elif db_url is not None:
        print('\n' + "Database methods not set up yet...")
        return 0

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.station + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    # Define DOA values
    back_az_vals = np.arange(back_az_min, back_az_max, back_az_step)
    trc_vel_vals = np.arange(trace_vel_min, trace_vel_max, trace_vel_step)

    # Trim streams if necessary
    # if signal_start is not None:
    #  stream.trim()  

    # run fk analysis
    beam_times, beam_peaks = fkd.run_fk(stream, [freq_min, freq_max], window_len, sub_window_len, window_step, method, back_az_vals, trc_vel_vals, pl)

    print('\n' + "Write results to specified output..." + '\n')
    if local_fk_out is not None:
        np.save(local_fk_out + ".fk_times", beam_times)
        np.save(local_fk_out + ".fk_peaks", beam_peaks)

        # need to add a [...].fk_meta.txt output that summarizes run parameters
        file_out = open(local_fk_out + ".fk_meta.txt", 'w')
        print("InfraPy Beamforming (fk) Analysis", file=file_out)
        print("---------------------------------" + '\n', file=file_out)

        print('\n' + "Data summary:", file=file_out)
        for tr in stream:
            print(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.station + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime), file=file_out)

        print('\n' + "Algorithm parameters:", file=file_out)
        print("  freq_min: " + str(freq_min), file=file_out)
        print("  freq_max: " + str(freq_max), file=file_out)
        print("  back_az_min: " + str(back_az_min), file=file_out)
        print("  back_az_max: " + str(back_az_max), file=file_out)
        print("  back_az_step: " + str(back_az_step), file=file_out)
        print("  trace_vel_min: " + str(trace_vel_min), file=file_out)
        print("  trace_vel_max: " + str(trace_vel_max), file=file_out)
        print("  trace_vel_step: " + str(trace_vel_step), file=file_out)
        print("  method: " + str(method), file=file_out)
        print("  signal_start: " + str(signal_start), file=file_out)
        print("  signal_end: " + str(signal_end), file=file_out)
        if method == "GLS":
            print("  noise_start: " + str(noise_start), file=file_out)
            print("  noise_end: " + str(noise_end), file=file_out)
        print("  window_len: " + str(window_len), file=file_out)
        print("  sub_window_len: " + str(sub_window_len), file=file_out)
        print("  window_step: " + str(window_step), file=file_out)
        print("  multithread: " + str(multithread), file=file_out)
        if multithread:
            print("  cpu_cnt: " + str(cpu_cnt), file=file_out)
        file_out.close()

    pl.close()
    pl.terminate()


@main.command('run_fd', short_help="Identify detections from beamforming results")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-fk", help="Local beamforming (fk) data files", default=None)
@click.option("--local-dets", help="Local detection data files", default=None)
@click.option("--window_len", help="Adaptive window length (default: " + default_config['FD']['window_len'] + " [s])", default=None, type=float)
@click.option("--p-value", help="Detection p-value (default: " + default_config['FD']['p_value'] + ")", default=None, type=float)
@click.option("--min-duration", help="Minimum detection duration (default: " + default_config['FD']['min_duration'] + " [s])", default=None, type=float)
@click.option("--back-az-width", help="Maximum azimuth scatter (default: " + default_config['FD']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--fixed-thresh", help="Fixed f-stat threshold (default: None)", default=None, type=float)
@click.option("--return-thresh", help="Return threshold (default: " + default_config['FD']['return_thresh'] + ")", default=None, type=bool)
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
    local_fk = set_param(user_config, 'DETECTION IO', 'local_fk', local_fk, 'string')
    local_dets = set_param(user_config, 'DETECTION IO', 'local_dets', local_dets, 'string')

    click.echo('\n' + "Data parameters:")
    click.echo("  local_fk: " + local_fk)
    click.echo("  local_dets: " + local_dets)

    # Algorithm parameters
    window_len = set_param(user_config, 'FD', 'window_len', window_len, 'float')
    p_value = set_param(user_config, 'FD', 'p_value', p_value, 'float')
    min_duration = set_param(user_config, 'FD', 'min_duration', min_duration, 'float')
    back_az_width = set_param(user_config, 'FD', 'back_az_width', back_az_width, 'float')
    fixed_thresh = set_param(user_config, 'FD', 'fixed_thresh', fixed_thresh, 'float')
    return_thresh = set_param(user_config, 'FD', 'return_thresh', return_thresh, 'bool')

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





@main.command('run_assoc', short_help="Associate detections into events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--back_az_width", help="Width of beam projection (default: " + default_config['ASSOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Maximum source-receiver range (default: " + default_config['ASSOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + default_config['ASSOC']['resolution'] + ")", default=None, type=int)
@click.option("--distance-matrix-max", help="Distance matrix maximum (default: " + default_config['ASSOC']['distance_matrix_max'] + ")", default=None, type=float)
@click.option("--cluster-threshold", help="Cluster linkage threshold (default: " + default_config['ASSOC']['cluster_threshold'] + ")", default=None, type=float)
@click.option("--trimming-threshold", help="Mishapen cluster threshold (default: " + default_config['ASSOC']['trimming_threshold'] + ")", default=None, type=float)
@click.option("--multithread", help="Use multithreading (default: " + default_config['FK']['multithread'] + ")", default=None, type=bool)
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
    back_az_width = set_param(user_config, 'ASSOC', 'back_az_width', back_az_width, 'float')
    range_max = set_param(user_config, 'ASSOC', 'range_max', range_max, 'float')
    resolution = set_param(user_config, 'ASSOC', 'resolution', resolution, 'int')
    distance_matrix_max = set_param(user_config, 'ASSOC', 'distance_matrix_max', distance_matrix_max, 'float')
    cluster_threshold = set_param(user_config, 'ASSOC', 'cluster_threshold', cluster_threshold, 'float')
    trimming_threshold = set_param(user_config, 'ASSOC', 'trimming_threshold', trimming_threshold, 'float')
    multithread = set_param(user_config, 'ASSOC', 'multithread', multithread, 'bool')
    cpu_cnt = set_param(user_config, 'ASSOC', 'cpu_cnt', cpu_cnt, 'int')

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



@main.command('run_loc', short_help="Estimate source locations and times for events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--back_az_width", help="Width of beam projection (default: " + default_config['LOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Maximum source-receiver range (default: " + default_config['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + default_config['LOC']['resolution'] + ")", default=None, type=int)
@click.option("--pgm-model", help="Path geometry model file (default: None)", default=None)
def run_loc(config_file, back_az_width, range_max, resolution, pgm_model):
    '''
    Estimate source locations and times for events


    More detailed description...
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
        click.echo("Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Data IO parameters
    click.echo('\n' + "Data summary:")
    click.echo("  NOT SET UP YET...")

    # Algorithm parameters
    back_az_width = set_param(user_config, 'LOC', 'back_az_width', back_az_width, 'float')
    range_max = set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    resolution = set_param(user_config, 'LOC', 'resolution', resolution, 'int')
    pgm_model = set_param(user_config, 'LOC', 'pgm_model', pgm_model, 'int')

    click.echo('\n' + "Parameter summary:")
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  range_max: " + str(range_max))
    click.echo("  resolution: " + str(resolution))
    click.echo("  pgm_model: " + str(pgm_model))

    print('\n' + "Run localization analysis...")

    print('\n' + "Write results to specified output..." + '\n')


@main.command('run_loc', short_help="Estimate source locations and times for events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--back_az_width", help="Width of beam projection (default: " + default_config['LOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Maximum source-receiver range (default: " + default_config['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + default_config['LOC']['resolution'] + ")", default=None, type=int)
@click.option("--pgm-model", help="Path geometry model file (default: None)", default=None)
def run_loc(config_file, back_az_width, range_max, resolution, pgm_model):
    '''
    Estimate source locations and times for events


    More detailed description...
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
        click.echo("Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Data IO parameters
    click.echo('\n' + "Data summary:")
    click.echo("  NOT SET UP YET...")

    # Algorithm parameters
    back_az_width = set_param(user_config, 'LOC', 'back_az_width', back_az_width, 'float')
    range_max = set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    resolution = set_param(user_config, 'LOC', 'resolution', resolution, 'int')
    pgm_model = set_param(user_config, 'LOC', 'pgm_model', pgm_model, 'int')

    click.echo('\n' + "Parameter summary:")
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  range_max: " + str(range_max))
    click.echo("  resolution: " + str(resolution))
    click.echo("  pgm_model: " + str(pgm_model))

    print('\n' + "Run localization analysis...")

    print('\n' + "Write results to specified output..." + '\n')

if __name__ == '__main__':
    main()
