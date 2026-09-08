#!/usr/bin/env python

import os 
import click
import warnings

import json
import gzip

import configparser as cnfg
import numpy as np

from multiprocessing import Pool

from obspy import UTCDateTime 

from ..utils import config, data_io
from ..detection import beam as fkd
from ..detection import spectral

@click.command('beam', short_help="Run beamforming-based detection on an array")
@click.option("--cnfg-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)

@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--det-label", help="Label for detection results", default=None)

@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--back-az-min", help="Minimum back azimuth (default: " + config.defaults['FK']['back_az_min'] + " [deg])", default=None, type=float)
@click.option("--back-az-max", help="Maximum back azimuth (default: " + config.defaults['FK']['back_az_max'] + " [deg])", default=None, type=float)
@click.option("--back-az-step", help="Back azimuth resolution (default: " + config.defaults['FK']['back_az_step'] + " [deg])", default=None, type=float)
@click.option("--trace-vel-min", help="Minimum trace velocity (default: " + config.defaults['FK']['trace_vel_min'] + " [m/s])", default=None, type=float)
@click.option("--trace-vel-max", help="Maximum trace velocity (default: " + config.defaults['FK']['trace_vel_max'] + " [m/s])", default=None, type=float)
@click.option("--trace-vel-step", help="Trace velocity resolution (default: " + config.defaults['FK']['trace_vel_step'] + " [m/s])", default=None, type=float)
@click.option("--method", help="Beamforming method (default: " + config.defaults['FK']['method'] + ")", default=None)
@click.option("--signal-start", help="Start of analysis window", default=None)
@click.option("--signal-end", help="End of analysis window", default=None)
@click.option("--noise-start", help="Start of noise sample", default=None)
@click.option("--noise-end", help="End of noise sample", default=None)
@click.option("--fk-window-len", help="Analysis window length (default: " + config.defaults['FK']['window_len'] + " [s])", default=None, type=float)
@click.option("--fk-sub-window-len", help="Analysis sub-window length (default: None [s])", default=None, type=float)
@click.option("--fk-window-step", help="Step between analysis windows (default: " + config.defaults['FK']['window_step'] + " [s])", default=None, type=float)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)

@click.option("--fd-window-len", help="Adaptive window length (default: " + config.defaults['FD']['window_len'] + " [s])", default=None, type=float)
@click.option("--p-value", help="Detection p-value (default: " + config.defaults['FD']['p_value'] + ")", default=None, type=float)
@click.option("--min-duration", help="Minimum detection duration (default: " + config.defaults['FD']['min_duration'] + " [s])", default=None, type=float)
@click.option("--back-az-width", help="Maximum azimuth scatter (default: " + config.defaults['FD']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--fixed-thresh", help="Fixed f-stat threshold (default: None)", default=None, type=float)
@click.option("--thresh-ceil", help="Hybrid f-stat threshold (default: None)", default=None, type=float)
@click.option("--merge-dets", help="Merge detections (default: " + config.defaults['FD']['merge_dets'] + ")", default=None, type=bool)
@click.option("--auto-overwrite", help="Automatically overwrite existing results", default=None, type=bool)
def run_beam_detect(cnfg_file, local_wvfrms, fdsn, db_config, local_latlon, network, station, location, channel, starttime, endtime, 
    det_label, freq_min, freq_max, back_az_min, back_az_max, back_az_step, trace_vel_min, trace_vel_max, trace_vel_step, method, signal_start, 
    signal_end, noise_start, noise_end, fk_window_len, fk_sub_window_len, fk_window_step, cpu_cnt, fd_window_len, p_value, min_duration, 
    back_az_width, fixed_thresh, thresh_ceil, merge_dets, auto_overwrite):
    '''
    Run combined beamforming (fk) and detection analysis to identify detection in array waveform data.
    
    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy detect beam --local-wvfrms 'data/YJ.BRP*.SAC' --cpu-cnt 4
    \tinfrapy detect beam --cnfg-file config/detection_local.config --cpu-cnt 4
    \tinfrapy detect beam --cnfg-file config/detection_fdsn.config --cpu-cnt 4

    '''
    
    click.echo("")
    click.echo("######################################")
    click.echo("##                                  ##")
    click.echo("##              InfraPy             ##")
    click.echo("##  Beamforming Detection Analyses  ##")
    click.echo("##                                  ##")
    click.echo("######################################")
    click.echo("")    

    if cnfg_file:
        click.echo('\n' + "Loading configuration info from: " + cnfg_file)
        if os.path.isfile(cnfg_file):
            user_config = cnfg.ConfigParser()
            user_config.read(cnfg_file)
        else:
            click.echo('\n' + "Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Database configuration and info   
    db_config = config.set_param(user_config, 'DATA IO', 'db_config', db_config, 'string')
    db_info = None

    # Local DATA IO parameters
    local_wvfrms = config.set_param(user_config, 'DATA IO', 'local_wvfrms', local_wvfrms, 'string')
    local_latlon = config.set_param(user_config, 'DATA IO', 'local_latlon', local_latlon, 'string')

    # FDSN DATA IO parameters
    fdsn = config.set_param(user_config, 'DATA IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'DATA IO', 'network', network, 'string')
    station = config.set_param(user_config, 'DATA IO', 'station', station, 'string')
    location = config.set_param(user_config, 'DATA IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'DATA IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'DATA IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'DATA IO', 'endtime', endtime, 'string')

    # Result IO
    det_label = config.set_param(user_config, 'DATA IO', 'det_label', det_label, 'string')

    click.echo('\n' + "Data parameters:")
    if local_wvfrms is not None:
        click.echo("  local_wvfrms: " + str(local_wvfrms))
        click.echo("  local_latlon: " + str(local_latlon))
    elif fdsn is not None:
        click.echo("  fdsn: " + str(fdsn))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    elif db_config is not None:
        db_info = cnfg.ConfigParser()
        db_info.read(db_config)
        click.echo("  db_config: " + str(db_config))
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
        
    click.echo("  det_label: " + str(det_label))

    # Algorithm parameters
    fk_params = {}
    fk_params['freq_min'] = config.set_param(user_config, 'FK', 'freq_min', freq_min, 'float')
    fk_params['freq_max'] = config.set_param(user_config, 'FK', 'freq_max', freq_max, 'float')
    fk_params['back_az_min'] = config.set_param(user_config, 'FK', 'back_az_min', back_az_min, 'float')
    fk_params['back_az_max'] = config.set_param(user_config, 'FK', 'back_az_max', back_az_max, 'float')
    fk_params['back_az_step'] = config.set_param(user_config, 'FK', 'back_az_step', back_az_step, 'float')
    fk_params['trace_vel_min'] = config.set_param(user_config, 'FK', 'trace_vel_min', trace_vel_min, 'float')
    fk_params['trace_vel_max'] = config.set_param(user_config, 'FK', 'trace_vel_max', trace_vel_max, 'float')
    fk_params['trace_vel_step'] = config.set_param(user_config, 'FK', 'trace_vel_step', trace_vel_step, 'float')
    fk_params['method'] = config.set_param(user_config, 'FK', 'method', method, 'string')
    fk_params['signal_start'] = config.set_param(user_config, 'FK', 'signal_start', signal_start, 'string')
    fk_params['signal_end'] = config.set_param(user_config, 'FK', 'signal_end', signal_end, 'string')
    fk_params['noise_start'] = config.set_param(user_config, 'FK', 'noise_start', noise_start, 'string')
    fk_params['noise_end'] = config.set_param(user_config, 'FK', 'noise_end', noise_end, 'string')
    fk_params['window_len'] = config.set_param(user_config, 'FK', 'window_len', fk_window_len, 'float')
    fk_params['sub_window_len'] = config.set_param(user_config, 'FK', 'sub_window_len', fk_sub_window_len, 'float')
    fk_params['window_step'] = config.set_param(user_config, 'FK', 'window_step', fk_window_step, 'float')
    fk_params['cpu_cnt'] = config.set_param(user_config, 'FK', 'cpu_cnt', cpu_cnt, 'int')

    if fk_params['cpu_cnt'] is not None:
        pl = Pool(fk_params['cpu_cnt'])
    else:
        pl = None

    click.echo('\n' + "fk (beam) parameters:")
    for key in fk_params.keys():
        if fk_params[key] is not None:
            click.echo("  " + key + ": " + str(fk_params[key]))

    det_params = {}
    det_params['window_len'] = config.set_param(user_config, 'FD', 'window_len', fd_window_len, 'float')
    det_params['p_value'] = config.set_param(user_config, 'FD', 'p_value', p_value, 'float')
    det_params['min_duration'] = config.set_param(user_config, 'FD', 'min_duration', min_duration, 'float')
    det_params['back_az_width'] = config.set_param(user_config, 'FD', 'back_az_width', back_az_width, 'float')
    det_params['fixed_thresh'] = config.set_param(user_config, 'FD', 'fixed_thresh', fixed_thresh, 'float')
    det_params['thresh_ceil'] = config.set_param(user_config, 'FD', 'thresh_ceil', thresh_ceil, 'float')
    det_params['merge_dets'] = config.set_param(user_config, 'FD', 'merge_dets', merge_dets, 'bool')

    click.echo('\n' + "detection parameters:")
    for key in det_params.keys():
        if det_params[key] is not None:
            click.echo("  " + key + ": " + str(det_params[key]))

    # Read in data
    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if using a noise window for analysis (only used for GLS analysis)
    ns_covar_inv = None
    if noise_start is not None:
        click.echo('\n' + "Analyzing noise window to compute background covariance...")
        click.echo('\t' + "noise start: " + fk_params['noise_start'])
        click.echo('\t' + "noise end: " + fk_params['noise_end'])

        st_noise = stream.copy()
        st_noise.trim(UTCDateTime(fk_params['noise_start']), UTCDateTime(fk_params['noise_end']))

        # Compute noise covariance
        x, t, _, _ = fkd.stream_to_array_data(st_noise, latlon=latlon)
        _, S, _ = fkd.fft_array_data(x, t, sub_window_len=fk_params['window_len'])

        ns_covar_inv = np.empty_like(S)
        for n in range(S.shape[2]):
            S[:, :, n] += 1.0e-3 * np.mean(np.diag(S[:, :, n])) * np.eye(S.shape[0])
            ns_covar_inv[:, :, n] = np.linalg.inv(S[:, :, n])

    # Check if using a signal window
    if fk_params["signal_start"] is not None or fk_params["signal_end"] is not None:
        if fk_params["signal_start"] is not None:
            t1 = UTCDateTime(fk_params["signal_start"])
        else:
            t1 = stream[0].stats.starttime
    
        if fk_params["signal_end"] is not None:
            t2 = UTCDateTime(fk_params["signal_end"])
        else:
            t2 = stream[0].stats.endtime

        if t1 > t2:
            warning_message = "Specified signal_start after signal_end. Stream won't be trimmed."
            warnings.warn((warning_message))
        else:
            if t1 < stream[0].stats.starttime:
                warning_message = "Specified signal_start before data start time."
                warnings.warn((warning_message))
                t1 = stream[0].stats.starttime 
        
            if t2 > stream[0].stats.endtime:
                warning_message = "Specified signal_end after data end time."
                warnings.warn((warning_message))
                t2 = stream[0].stats.endtime 
        
            click.echo('\n' + "Trimming data to signal analysis window...")
            click.echo('\t' + "start time: " + str(t1))
            click.echo('\t' + "end time: " + str(t2))
            stream.trim(t1, t2)

    wvfrm_info = data_io.wvfrm_info(stream, latlon)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    if local_wvfrms is not None and "/" in local_wvfrms:
        output_id = os.path.dirname(local_wvfrms) + "/"
    else:
        output_id = ""
    output_id = output_id + data_io.stream_label(stream)

    if det_label is None or det_label == "auto":
        det_label = output_id

    # check if results already exist
    if not auto_overwrite and os.path.isfile(det_label + ".dets.json.gz"):
        user_opt = input('\nWARNING!!! Detection results file (' + det_label + '.dets.json.gz) already exists and will be overwritten. \nDo you want to proceed? (y/n): ').lower().strip()
        while True:
            if user_opt in ['y', 'yes']:
                break
            elif user_opt in ['n', 'no']:
                click.echo("")
                return 
            else:
                user_opt = input('Invalid input. Proceed and overwrite detection file? (y/n): ').lower().strip()

    # Run beamforming (fk)
    beam_times, beam_peaks = fkd.run_fk_dict(stream, latlon, fk_params, ns_covar_inv, pl)

    print("Running adaptive f-detector..." + '\n')
    dets, thresh_vals = fkd.run_afd_dict(beam_times, beam_peaks, fk_params, det_params, len(stream))

    # save fk results for the full duration
    dt = np.array([(tn - np.datetime64(stream[0].stats.starttime)).astype('m8[ms]').astype(float) * 1.0e-3 for tn in beam_times])

    fk_out = {}
    fk_out['time'] = dt
    fk_out['back az'] = beam_peaks[:, 0]
    fk_out['tr vel'] = beam_peaks[:, 1]
    fk_out['f-stat'] = beam_peaks[:, 2]
    fk_out['thresh'] = thresh_vals 

    dets_out = [fkd.det2dict(stream, latlon, beam_times, beam_peaks, fk_params, det_info) for det_info in dets] 
    
    click.echo("Writing beamforming and detection results into " + det_label + ".dets.json.gz" + '\n')
    det_output = {'wvfrm_info' : [wvfrm_info], 'fk_params' : [fk_params], 'det_params' : [det_params], 'fk' : fk_out, 'det_info' : dets_out}
    with gzip.open(det_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

    if pl is not None:
        pl.terminate()
        pl.close()
        

@click.command('spectral', short_help="Run spectral detection on a single channel")
@click.option("--cnfg-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)

@click.option("--local-latlon", help="Location information for local waveforms", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)

@click.option("--det-label", help="Label for detection results", default=None)

@click.option("--signal-start", help="Start of analysis window", default=None)
@click.option("--signal-end", help="End of analysis window", default=None)

@click.option("--spectral-option", help="Spectral analysis method ('spectogram', 'stft', or 'cwt'), default: " + config.defaults['SD']['spectral_option'] + ")", default=None)
@click.option("--morlet-omega0", help="Morlet parameter for 'cwt', default: " + config.defaults['SD']['morlet_omega0'] + ")", default=None, type=float)

@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--window-len", help="Adaptive window length (default: " + config.defaults['SD']['window_len'] + " [s])", default=None, type=float)
@click.option("--window-step", help="Adaptive window step (default: " + config.defaults['SD']['window_step'] + " [s])", default=None, type=float)
@click.option("--p-value", help="Detection p-value (default: " + config.defaults['SD']['p_value'] + ")", default=None, type=float)
@click.option("--freq-tm-factor", help="Freq./time scaling (sec/decade) (def.: " + config.defaults['SD']['freq_tm_factor'] + ")", default=None, type=float)

@click.option("--cluster-eps", help="Clustering linkage distance (default: " + config.defaults['SD']['cluster_eps'] + ")", default=None, type=float)
@click.option("--cluster-min-samples", help="Clustering minimum samples (default: " + config.defaults['SD']['cluster_min_samples'] + ")", default=None, type=int)
@click.option("--cluster-window-len", help="Clustering linkage distance (default: " + config.defaults['SD']['cluster_window_len'] + ")", default=None, type=float)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
@click.option("--auto-overwrite", help="Automatically overwrite existing results", default=None, type=bool)
def run_spec_detect(cnfg_file, local_wvfrms, fdsn, db_config, local_latlon, network, station, location, channel, starttime, endtime, 
    det_label, signal_start, signal_end, spectral_option, morlet_omega0, freq_min, freq_max, window_len, window_step, 
    p_value, freq_tm_factor, cluster_eps, cluster_min_samples, cluster_window_len, cpu_cnt, auto_overwrite):
    '''
    Run spectral detection methods on a single channel to identify signals of interest.
    
    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy detect spectral --local-wvfrms 'data/YJ.BRP1..EDF.SAC' --cpu-cnt 4   
    \tinfrapy detect spectral --local-wvfrms 'data/YJ.BRP1..EDF.SAC' --cpu-cnt 4 --spectral-option cwt --cluster-min-samples 500 --cluster-eps 4
    
    '''
    

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##   Spectral Detection Analyses   ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("") 

    if cnfg_file:
        click.echo('\n' + "Loading configuration info from: " + cnfg_file)
        if os.path.isfile(cnfg_file):
            user_config = cnfg.ConfigParser()
            user_config.read(cnfg_file)
        else:
            click.echo('\n' + "Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Database configuration and info   
    db_config = config.set_param(user_config, 'DATA IO', 'db_config', db_config, 'string')
    db_info = None

    # Local DATA IO parameters
    local_wvfrms = config.set_param(user_config, 'DATA IO', 'local_wvfrms', local_wvfrms, 'string')
    local_latlon = config.set_param(user_config, 'DATA IO', 'local_latlon', local_latlon, 'string')

    # FDSN DATA IO parameters
    fdsn = config.set_param(user_config, 'DATA IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'DATA IO', 'network', network, 'string')
    station = config.set_param(user_config, 'DATA IO', 'station', station, 'string')
    location = config.set_param(user_config, 'DATA IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'DATA IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'DATA IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'DATA IO', 'endtime', endtime, 'string')

    # Result IO
    det_label = config.set_param(user_config, 'DATA IO', 'det_label', det_label, 'string')

    click.echo('\n' + "Data parameters:")
    if local_wvfrms is not None:
        click.echo("  local_wvfrms: " + str(local_wvfrms))
        click.echo("  local_latlon: " + str(local_latlon))
    elif fdsn is not None:
        click.echo("  fdsn: " + str(fdsn))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    elif db_config is not None:
        db_info = cnfg.ConfigParser()
        db_info.read(db_config)
        click.echo("  db_config: " + str(db_config))
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
        
    click.echo("  det_label: " + str(det_label))
    if cpu_cnt is not None:
        click.echo("  cpu_cnt: " + str(cpu_cnt))
        pl = Pool(cpu_cnt)
    else:
        pl = None

    # Algorithm parameters
    sd_params = {}
    sd_params["spectral_option"] = config.set_param(user_config, 'SD', 'spectral_option', spectral_option, 'string')
    sd_params["morlet_omega0"] = config.set_param(user_config, 'SD', 'morlet_omega0', morlet_omega0, 'float')    
    sd_params["freq_min"] = config.set_param(user_config, 'SD', 'freq_min', freq_min, 'float')
    sd_params["freq_max"] = config.set_param(user_config, 'SD', 'freq_max', freq_max, 'float')
    sd_params["signal_start"] = config.set_param(user_config, 'SD', 'signal_start', signal_start, 'string')
    sd_params["signal_end"] = config.set_param(user_config, 'SD', 'signal_end', signal_end, 'string')
    sd_params["window_len"] = config.set_param(user_config, 'SD', 'window_len', window_len, 'float')
    sd_params["window_step"] = config.set_param(user_config, 'SD', 'window_step', window_step, 'float')
    sd_params["p_value"] = config.set_param(user_config, 'SD', 'p_value', p_value, 'float')
    sd_params["freq_tm_factor"] = config.set_param(user_config, 'SD', 'freq_tm_factor', freq_tm_factor, 'float')
    sd_params["cluster_eps"] = config.set_param(user_config, 'SD', 'cluster_eps', cluster_eps, 'float')
    sd_params["cluster_min_samples"] = config.set_param(user_config, 'SD', 'cluster_min_samples', cluster_min_samples, 'int')
    sd_params["cluster_window_len"] = config.set_param(user_config, 'SD', 'cluster_window_len', cluster_window_len, 'float')
    sd_params["cpu_cnt"] = config.set_param(user_config, 'SD', 'cpu_cnt', cpu_cnt, 'int')

    if sd_params['cpu_cnt'] is not None:
        pl = Pool(sd_params['cpu_cnt'])
    else:
        pl = None

    if 'cwt' not in sd_params['spectral_option']:
        sd_params['morlet_omega0'] = None

    click.echo('\n' + "sd (spectral detector) parameters:")
    for key in sd_params.keys():
        if sd_params[key] is not None:
            click.echo("  " + key + ": " + str(sd_params[key]))

    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if using a signal window
    if sd_params["signal_start"] is not None or sd_params["signal_end"] is not None:
        if sd_params["signal_start"] is not None:
            t1 = UTCDateTime(sd_params["signal_start"])
        else:
            t1 = stream[0].stats.starttime
    
        if sd_params["signal_end"] is not None:
            t2 = UTCDateTime(sd_params["signal_end"])
        else:
            t2 = stream[0].stats.endtime

        if t1 > t2:
            warning_message = "Specified signal_start after signal_end. Stream won't be trimmed."
            warnings.warn((warning_message))
        else:
            if t1 < stream[0].stats.starttime:
                warning_message = "Specified signal_start before data start time."
                warnings.warn((warning_message))
                t1 = stream[0].stats.starttime 
        
            if t2 > stream[0].stats.endtime:
                warning_message = "Specified signal_end after data end time."
                warnings.warn((warning_message))
                t2 = stream[0].stats.endtime 
        
            click.echo('\n' + "Trimming data to signal analysis window...")
            click.echo('\t' + "start time: " + str(t1))
            click.echo('\t' + "end time: " + str(t2))
            stream.trim(t1, t2)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    wvfrm_info = data_io.wvfrm_info(stream, latlon)

    if local_wvfrms is not None and "/" in local_wvfrms:
        output_id = os.path.dirname(local_wvfrms) + "/"
    else:
        output_id = ""
    output_id = output_id + data_io.stream_label(stream)

    if det_label is None or det_label == "auto":
        det_label = output_id

    # check if results already exist
    if not auto_overwrite and os.path.isfile(det_label + ".dets.json.gz"):
        user_opt = input('\nWARNING!!! Detection results file (' + det_label + '.dets.json.gz) already exists and will be overwritten. \nDo you want to proceed? (y/n): ').lower().strip()
        while True:
            if user_opt in ['y', 'yes']:
                break
            elif user_opt in ['n', 'no']:
                click.echo("")
                return 
            else:
                user_opt = input('Invalid input. Proceed and overwrite detection file? (y/n): ').lower().strip()

    det_list, spectrogram, history = spectral.spec_det_dict(stream[0], sd_params, pl)
    det_output = {'wvfrm_info' : wvfrm_info, 'sd_params' : sd_params, 'spectrogram': spectrogram, 'history': history, 'det_info' : det_list}
    with gzip.open(det_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

    if pl is not None:
        pl.terminate()
        pl.close()

