#!/usr/bin/env python
import sys

import click
import configparser as cnfg
import numpy as np

from multiprocessing import Pool
from obspy.core import read as obspy_read

from ..utils import config
from ..utils import data_io
from ..detection import beamforming_new as fkd

@click.command('run_fk', short_help="Run beamforming methods on waveform data")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-url", help="Database URL for waveform data files", default=None)
@click.option("--db-site", help="Database site table for waveform data files", default=None)
@click.option("--db-wfdisc", help="Database wfdisc table for waveform data files", default=None)
@click.option("--db-origin", help="Database origin table for waveform data files", default=None)

@click.option("--local_latlon", help="Array location information for local waveforms", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)

@click.option("--local-fk-out", help="Local beamforming (fk) data files", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--back-az-min", help="Minimum back azimuth (default: " + config.defaults['FK']['back_az_min'] + " [deg])", default=None, type=float)
@click.option("--back-az-max", help="Maximum back azimuth (default: " + config.defaults['FK']['back_az_max'] + " [deg])", default=None, type=float)
@click.option("--back-az-step", help="Back azimuth resolution (default: " + config.defaults['FK']['back_az_step'] + " [deg])", default=None, type=float)
@click.option("--trace-vel-min", help="Minimum trace velocity (default: " + config.defaults['FK']['trace_vel_min'] + " [m/s])", default=None, type=float)
@click.option("--trace-vel-max", help="Maximum trace velocity (default: " + config.defaults['FK']['trace_vel_max'] + " [m/s])", default=None, type=float)
@click.option("--trace-vel-step", help="Trace velocity resolution (default: " + config.defaults['FK']['trace_vel_step'] + " [m/s])", default=None, type=float)
@click.option("--method", help="Beamforming method (default: " + config.defaults['FK']['method'] + ")", default=None)
@click.option("--signal_start", help="Start of analysis window (rel. starttime [s])", default=None, type=float)
@click.option("--signal_end", help="End of analysis window (rel. starttime [s])", default=None, type=float)
@click.option("--noise_start", help="Start of noise sample (see GLS algorithm notes)", default=None, type=float)
@click.option("--noise_end", help="End of noise sample (see GLS algorithm notes)", default=None, type=float)
@click.option("--window_len", help="Analysis window length (default: " + config.defaults['FK']['window_len'] + " [s])", default=None, type=float)
@click.option("--sub-window_len", help="Analysis sub-window length (default: None [s])", default=None, type=float)
@click.option("--window_step", help="Step between analysis windows (default: " + config.defaults['FK']['window_step'] + " [s])", default=None, type=float)
@click.option("--multithread", help="Use multithreading (default: " + config.defaults['FK']['multithread'] + ")", default=None, type=bool)
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
    db_url = config.set_param(user_config, 'database', 'url', db_url, 'string')
    db_site = config.set_param(user_config, 'database', 'site', db_site, 'string')
    db_wfdisc = config.set_param(user_config, 'database', 'wfdisc', db_wfdisc, 'string')
    db_origin = config.set_param(user_config, 'database', 'origin', db_origin, 'string')

    # Waveform IO parameters
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')
   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       
    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    # Result IO
    local_fk_out = config.set_param(user_config, 'DETECTION IO', 'local_fk_out', local_fk_out, 'string')

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
    freq_min = config.set_param(user_config, 'FK', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'FK', 'freq_max', freq_max, 'float')
    back_az_min = config.set_param(user_config, 'FK', 'back_az_min', back_az_min, 'float')
    back_az_max = config.set_param(user_config, 'FK', 'back_az_max', back_az_max, 'float')
    back_az_step = config.set_param(user_config, 'FK', 'back_az_step', back_az_step, 'float')
    trace_vel_min = config.set_param(user_config, 'FK', 'trace_vel_min', trace_vel_min, 'float')
    trace_vel_max = config.set_param(user_config, 'FK', 'trace_vel_max', trace_vel_max, 'float')
    trace_vel_step = config.set_param(user_config, 'FK', 'trace_vel_step', trace_vel_step, 'float')
    method = config.set_param(user_config, 'FK', 'method', method, 'string')
    signal_start = config.set_param(user_config, 'FK', 'signal_start', signal_start, 'float')
    signal_end = config.set_param(user_config, 'FK', 'signal_end', signal_end, 'float')
    noise_start = config.set_param(user_config, 'FK', 'noise_start', noise_start, 'float')
    noise_end = config.set_param(user_config, 'FK', 'noise_end', noise_end, 'float')
    window_len = config.set_param(user_config, 'FK', 'window_len', window_len, 'float')
    sub_window_len = config.set_param(user_config, 'FK', 'sub_window_len', sub_window_len, 'float')
    window_step = config.set_param(user_config, 'FK', 'window_step', window_step, 'float')
    multithread = config.set_param(user_config, 'FK', 'multithread', multithread, 'bool')
    cpu_cnt = config.set_param(user_config, 'FK', 'cpu_cnt', cpu_cnt, 'int')

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
    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_url, network, station, location, channel, starttime, endtime)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.station + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    # Define DOA values
    back_az_vals = np.arange(back_az_min, back_az_max, back_az_step)
    trc_vel_vals = np.arange(trace_vel_min, trace_vel_max, trace_vel_step)

    # run fk analysis
    beam_times, beam_peaks = fkd.run_fk(stream, latlon, [freq_min, freq_max], window_len, sub_window_len, window_step, method, back_az_vals, trc_vel_vals, pl)

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