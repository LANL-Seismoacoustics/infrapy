#!/usr/bin/env python

import os 
import click
import warnings

import configparser as cnfg
import numpy as np

from multiprocessing import Pool

from obspy import UTCDateTime 

from infrapy.utils import config
from infrapy.utils import data_io
from infrapy.detection import beamforming_new as fkd


@click.command('run_fk', short_help="Run beamforming methods on waveform data")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-url", help="Database URL for waveform data files", default=None)
@click.option("--db-site", help="Database site table for waveform data files", default=None)
@click.option("--db-wfdisc", help="Database wfdisc table for waveform data files", default=None)

@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)

@click.option("--local-fk-label", help="Label for local output of fk results", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--back-az-min", help="Minimum back azimuth (default: " + config.defaults['FK']['back_az_min'] + " [deg])", default=None, type=float)
@click.option("--back-az-max", help="Maximum back azimuth (default: " + config.defaults['FK']['back_az_max'] + " [deg])", default=None, type=float)
@click.option("--back-az-step", help="Back azimuth resolution (default: " + config.defaults['FK']['back_az_step'] + " [deg])", default=None, type=float)
@click.option("--trace-vel-min", help="Minimum trace velocity (default: " + config.defaults['FK']['trace_vel_min'] + " [m/s])", default=None, type=float)
@click.option("--trace-vel-max", help="Maximum trace velocity (default: " + config.defaults['FK']['trace_vel_max'] + " [m/s])", default=None, type=float)
@click.option("--trace-vel-step", help="Trace velocity resolution (default: " + config.defaults['FK']['trace_vel_step'] + " [m/s])", default=None, type=float)
@click.option("--method", help="Beamforming method (default: " + config.defaults['FK']['method'] + ")", default=None)
@click.option("--signal-start", help="Start of signal window", default=None)
@click.option("--signal-end", help="End of signal window", default=None)
@click.option("--noise-start", help="Start of noise sample", default=None)
@click.option("--noise-end", help="End of noise sample", default=None)
@click.option("--window-len", help="Analysis window length (default: " + config.defaults['FK']['window_len'] + " [s])", default=None, type=float)
@click.option("--sub-window-len", help="Analysis sub-window length (default: None [s])", default=None, type=float)
@click.option("--window-step", help="Step between analysis windows (default: " + config.defaults['FK']['window_step'] + " [s])", default=None, type=float)
@click.option("--multithread", help="Use multithreading (default: " + config.defaults['FK']['multithread'] + ")", default=None, type=bool)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
@click.option("--write-wvfrms", help="Write waveforms into local files (default: " + config.defaults['FK']['write_wvfrms'] + ")", default=None, type=bool)
def run_fk(config_file, local_wvfrms, fdsn, db_url, db_site, db_wfdisc, local_latlon, network, station, location, channel, starttime, endtime,
    local_fk_label, freq_min, freq_max, back_az_min, back_az_max, back_az_step, trace_vel_min, trace_vel_max, trace_vel_step, method, 
    signal_start, signal_end, noise_start, noise_end, window_len, sub_window_len, window_step, multithread, cpu_cnt, write_wvfrms):
    '''
    Run beamforming (fk) analysis

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_fk --local-wvfrms 'data/YJ.BRP*.SAC' --cpu-cnt 4
    \tinfrapy run_fk --config-file config/detection_local.config --cpu-cnt 4
    \tinfrapy run_fk --config-file config/detection_fdsn.config --cpu-cnt 4

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
    db_url = config.set_param(user_config, 'WAVEFORM IO', 'db_url', db_url, 'string')
    db_site = config.set_param(user_config, 'WAVEFORM IO', 'db_site', db_site, 'string')
    db_wfdisc = config.set_param(user_config, 'WAVEFORM IO', 'db_wfdisc', db_wfdisc, 'string')

    # Local waveform IO parameters
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    local_latlon = config.set_param(user_config, 'WAVEFORM IO', 'local_latlon', local_latlon, 'string')

    # FDSN waveform IO parameters
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    # Result IO
    local_fk_label = config.set_param(user_config, 'DETECTION IO', 'local_fk_label', local_fk_label, 'string')

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
    elif db_url is not None:
        click.echo("  db_url: " + str(db_url))
        click.echo("  db_site: " + str(db_site))
        click.echo("  db_wfdisc: " + str(db_wfdisc))
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
        
    click.echo("  local_fk_label: " + str(local_fk_label))

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
    signal_start = config.set_param(user_config, 'FK', 'signal_start', signal_start, 'string')
    signal_end = config.set_param(user_config, 'FK', 'signal_end', signal_end, 'string')
    noise_start = config.set_param(user_config, 'FK', 'noise_start', noise_start, 'float')
    noise_end = config.set_param(user_config, 'FK', 'noise_end', noise_end, 'float')
    window_len = config.set_param(user_config, 'FK', 'window_len', window_len, 'float')
    sub_window_len = config.set_param(user_config, 'FK', 'sub_window_len', sub_window_len, 'float')
    window_step = config.set_param(user_config, 'FK', 'window_step', window_step, 'float')
    multithread = config.set_param(user_config, 'FK', 'multithread', multithread, 'bool')
    cpu_cnt = config.set_param(user_config, 'FK', 'cpu_cnt', cpu_cnt, 'int')
    write_wvfrms = config.set_param(user_config, 'FK', 'write_wvfrms', write_wvfrms, 'bool')

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
    click.echo("  write_wvfrms: " + str(write_wvfrms))

    # Check data option and populate obspy Stream
    if db_url is not None:
        db_info = {'url': db_url, 'site': db_site, 'wfdisc': db_wfdisc}
    else:
        db_info = None
    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    if write_wvfrms:
        if local_wvfrms is None:
            click.echo('\n' + "Writing waveform data to local SAC files...")
            data_io.write_stream_to_sac(stream, latlon)
        else: 
            click.echo('\n' + "Cannot write waveform data when using local data...")

    if local_fk_label is None or local_fk_label == "auto":
        if local_wvfrms is not None and "/" in local_wvfrms:
            local_fk_label = os.path.dirname(local_wvfrms) + "/"
        else:
            local_fk_label = ""

        local_fk_label = local_fk_label + stream[0].stats.network + "." + os.path.commonprefix([tr.stats.station for tr in stream])
        local_fk_label = local_fk_label + '_' + "%02d" % stream[0].stats.starttime.year + ".%02d" % stream[0].stats.starttime.month + ".%02d" % stream[0].stats.starttime.day
        local_fk_label = local_fk_label + '_' + "%02d" % stream[0].stats.starttime.hour + "." + "%02d" % stream[0].stats.starttime.minute + "." + "%02d" % stream[0].stats.starttime.second
        local_fk_label = local_fk_label + '-' + "%02d" % stream[0].stats.endtime.hour + "." + "%02d" % stream[0].stats.endtime.minute + "." + "%02d" % stream[0].stats.endtime.second

    # Define DOA values
    back_az_vals = np.arange(back_az_min, back_az_max, back_az_step)
    trc_vel_vals = np.arange(trace_vel_min, trace_vel_max, trace_vel_step)

    '''
    # Check if using a noise window
    if noise_start is not None:
        print("Analyzing noise window to compute noise covariance...")
        st_noise = stream.copy()
        st_noise.trim(UTCDateTime(noise_start), UTCDateTime(noise_end))

        # Compute noise covariance
        x, t, _, _ = fkd.stream_to_array_data(st_noise, latlon=latlon)
        _, S, _ = fkd.fft_array_data(x, t, sub_window_len=window_len)

        ns_covar_inv = np.empty_like(S)
        for n in range(S.shape[2]):
            S[:, :, n] += 1.0e-3 * np.mean(np.diag(S[:, :, n])) * np.eye(S.shape[0])
            ns_covar_inv[:, :, n] = np.linalg.inv(S[:, :, n])
    '''
            
    # Check if using a signal window
    if signal_start is not None:
        t1 = UTCDateTime(signal_start)
        t2 = UTCDateTime(signal_end)

        click.echo('\n' + "Trimming data to signal analysis window...")
        click.echo('\t' + "start time: " + str(t1))
        click.echo('\t' + "end time: " + str(t2))


        warning_message = "signal_start and signal_end values poorly defined."
        if t1 > t2:
            warning_message = warning_message + "  signal_start after signal_end."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        elif t1 < stream[0].stats.starttime:
            warning_message = warning_message + "  signal_start before data start time."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        elif t2 > stream[0].stats.endtime:
            warning_message = warning_message + "  signal_end after data end time."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        else:
            stream.trim(t1, t2)

    # run fk analysis
    beam_times, beam_peaks = fkd.run_fk(stream, latlon, [freq_min, freq_max], window_len, sub_window_len, window_step, method, back_az_vals, trc_vel_vals, pl)

    # new save methods
    dt = np.array([(tn - np.datetime64(tr.stats.starttime)).astype('m8[ms]').astype(float) * 1.0e-3 for tn in beam_times])
    fk_results = np.hstack((np.atleast_2d(dt).T, beam_peaks))
    fk_header = data_io.fk_header(stream, latlon, freq_min, freq_max, back_az_min, back_az_max, back_az_step, trace_vel_min, trace_vel_max, trace_vel_step, method, 
        signal_start, signal_end, noise_start, noise_end, window_len, sub_window_len, window_step)

    if not os.path.isfile(local_fk_label + ".fk_results.dat"):
        click.echo('\n' + "Writing results into " + local_fk_label + ".fk_results.dat")
        np.savetxt(local_fk_label + ".fk_results.dat", fk_results, header=fk_header)
    else:
        k = 0
        while os.path.isfile(local_fk_label + "-v" + str(k) + ".fk_results.dat"):
            k += 1
        click.echo('\n' + "WARNING!  fk results file(s) already exist." + '\n' + "Writing a new version: " + local_fk_label + "-v" + str(k) + ".fk_results.dat")
        np.savetxt(local_fk_label + "-v" + str(k) + ".fk_results.dat", fk_results, header=fk_header)

    if multithread or cpu_cnt is not None:
        pl.terminate()
        pl.close()

    click.echo('')


@click.command('run_fd', short_help="Identify detections from beamforming results")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-fk-label", help="Local beamforming (fk) results label", default=None)
@click.option("--local-detect-label", help="Label for local detection (fd) results", default=None)
@click.option("--window-len", help="Adaptive window length (default: " + config.defaults['FD']['window_len'] + " [s])", default=None, type=float)
@click.option("--p-value", help="Detection p-value (default: " + config.defaults['FD']['p_value'] + ")", default=None, type=float)
@click.option("--min-duration", help="Minimum detection duration (default: " + config.defaults['FD']['min_duration'] + " [s])", default=None, type=float)
@click.option("--back-az-width", help="Maximum azimuth scatter (default: " + config.defaults['FD']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--fixed-thresh", help="Fixed f-stat threshold (default: None)", default=None, type=float)
@click.option("--return-thresh", help="Return threshold (default: " + config.defaults['FD']['return_thresh'] + ")", default=None, type=bool)
def run_fd(config_file, local_fk_label, local_detect_label, window_len, p_value, min_duration, back_az_width, fixed_thresh, return_thresh):
    '''
    Run fd analysis to identify detections in beamforming results

    \b
    Example usage (run from infrapy/examples directory after run_fk examples):
    \tinfrapy run_fd --local-fk-label data/YJ.BRP_2012.04.09_18.00.00-18.19.59
    \tinfrapy run_fd --local-fk-label IM.I53H_2018.12.19_01.00.00-03.00.00
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
    local_fk_label = config.set_param(user_config, 'DETECTION IO', 'local_fk_label', local_fk_label, 'string')
    local_detect_label = config.set_param(user_config, 'DETECTION IO', 'local_detect_label', local_detect_label, 'string')

    if ".fk_results.dat" in local_fk_label:
        local_fk_label = local_fk_label[:-15]

    if local_detect_label is None or local_detect_label == "auto":
        local_detect_label = local_fk_label

    click.echo('\n' + "Data parameters:")
    click.echo("  local_fk_label: " + local_fk_label)
    click.echo("  local_detect_label: " + local_detect_label)

    # Algorithm parameters
    window_len = config.set_param(user_config, 'FD', 'window_len', window_len, 'float')
    p_value = config.set_param(user_config, 'FD', 'p_value', p_value, 'float')
    min_duration = config.set_param(user_config, 'FD', 'min_duration', min_duration, 'float')
    back_az_width = config.set_param(user_config, 'FD', 'back_az_width', back_az_width, 'float')
    fixed_thresh = config.set_param(user_config, 'FD', 'fixed_thresh', fixed_thresh, 'float')
    return_thresh = config.set_param(user_config, 'FD', 'return_thresh', return_thresh, 'bool')

    click.echo('\n' + "Algorithm parameters:")
    click.echo("  window_len: " + str(window_len))
    click.echo("  p_value: " + str(p_value))
    click.echo("  min_duration: " + str(min_duration))
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  fixed_thresh: " + str(fixed_thresh))
    click.echo("  return_thresh: " + str(return_thresh))

    print('\n' + "Running fd...")

    if local_fk_label is not None:
        temp = np.loadtxt(local_fk_label + ".fk_results.dat")
        dt, beam_peaks = temp[:, 0], temp[:, 1:]

        temp = open(local_fk_label + ".fk_results.dat", 'r')
        for line in temp:
            if "t0:" in line:
                t0 = np.datetime64(line.split(' ')[-1][:-1])
            elif "freq_min" in line:
                freq_min = float(line.split(' ')[-1])
            elif "freq_max" in line:
                freq_max = float(line.split(' ')[-1])
            elif "window_len" in line and "sub_window_len" not in line:
                fk_window_len = float(line.split(' ')[-1])
            elif "channel_cnt" in line:
                channel_cnt = float(line.split(' ')[-1])
            elif "latitude" in line:
                array_lat = float(line.split(' ')[-1])
            elif "longitude" in line:
                array_lon = float(line.split(' ')[-1])

        beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
    else:
        print("Non-local data not yet set up...")
        return 0

    TB_prod = (freq_max - freq_min) * fk_window_len
    min_seq = max(2, int(min_duration / fk_window_len))

    dets, thresh_vals = fkd.run_fd(beam_times, beam_peaks, window_len, TB_prod, channel_cnt, p_value, min_seq, back_az_width, fixed_thresh, True)

    det_list = []
    for det_info in dets:
        det_list = det_list + [data_io.define_detection(det_info, [array_lat, array_lon], channel_cnt, [freq_min,freq_max], note="InfraPy CLI detection")]
    print("Writing detections to " + local_detect_label + ".dets.json")
    data_io.detection_list_to_json(local_detect_label + ".dets.json", det_list)

    if return_thresh:
        np.save(local_detect_label + ".thresholds", thresh_vals)


@click.command('run_fkd', short_help="Run beamforming and detection methods in sequence")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-url", help="Database URL for waveform data files", default=None)
@click.option("--db-site", help="Database site table for waveform data files", default=None)
@click.option("--db-wfdisc", help="Database wfdisc table for waveform data files", default=None)

@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)

@click.option("--local-fk-label", help="Label for local output of fk results", default=None)
@click.option("--local-detect-label", help="Label for local detection (fd) results", default=None)

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
@click.option("--multithread", help="Use multithreading (default: False)", default=None, type=bool)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)

@click.option("--fd-window-len", help="Adaptive window length (default: " + config.defaults['FD']['window_len'] + " [s])", default=None, type=float)
@click.option("--p-value", help="Detection p-value (default: " + config.defaults['FD']['p_value'] + ")", default=None, type=float)
@click.option("--min-duration", help="Minimum detection duration (default: " + config.defaults['FD']['min_duration'] + " [s])", default=None, type=float)
@click.option("--back-az-width", help="Maximum azimuth scatter (default: " + config.defaults['FD']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--fixed-thresh", help="Fixed f-stat threshold (default: None)", default=None, type=float)
@click.option("--return-thresh", help="Return threshold (default: " + config.defaults['FD']['return_thresh'] + ")", default=None, type=bool)
@click.option("--write-wvfrms", help="Write waveforms into local SAC files (default: " + config.defaults['FK']['write_wvfrms'] + ")", default=None, type=bool)
def run_fkd(config_file, local_wvfrms, fdsn, db_url, db_site, db_wfdisc, local_latlon, network, station, location, channel, starttime, endtime,
    local_fk_label, local_detect_label, freq_min, freq_max, back_az_min, back_az_max, back_az_step, trace_vel_min, trace_vel_max, trace_vel_step, method,  signal_start, 
    signal_end, noise_start, noise_end, fk_window_len, fk_sub_window_len, fk_window_step, multithread, cpu_cnt, fd_window_len, p_value, min_duration, 
    back_az_width, fixed_thresh, return_thresh, write_wvfrms):
    '''
    Run combined beamforming (fk) and detection analysis to identify detection in array waveform data.
    
    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_fkd --local-wvfrms 'data/YJ.BRP*.SAC' --cpu-cnt 4
    \tinfrapy run_fkd --config-file config/detection_local.config --cpu-cnt 4
    \tinfrapy run_fkd --config-file config/detection_fdsn.config --cpu-cnt 4

    '''
    

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##    Combined Beamforming and     ##")
    click.echo("##  Detection (fk + fd) Analyses   ##")
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
    db_url = config.set_param(user_config, 'WAVEFORM IO', 'db_url', db_url, 'string')
    db_site = config.set_param(user_config, 'WAVEFORM IO', 'db_site', db_site, 'string')
    db_wfdisc = config.set_param(user_config, 'WAVEFORM IO', 'db_wfdisc', db_wfdisc, 'string')

    # Local waveform IO parameters
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    local_latlon = config.set_param(user_config, 'WAVEFORM IO', 'local_latlon', local_latlon, 'string')

    # FDSN waveform IO parameters
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    # Result IO
    local_fk_label = config.set_param(user_config, 'DETECTION IO', 'local_fk_label', local_fk_label, 'string')
    local_detect_label = config.set_param(user_config, 'DETECTION IO', 'local_detect_label', local_detect_label, 'string')

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
    elif db_url is not None:
        click.echo("  db_url: " + str(db_url))
        click.echo("  db_site: " + str(db_site))
        click.echo("  db_wfdisc: " + str(db_wfdisc))
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
        
    click.echo("  local_fk_label: " + str(local_fk_label))
    click.echo("  local_detect_label: " + str(local_detect_label))

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
    signal_start = config.set_param(user_config, 'FK', 'signal_start', signal_start, 'string')
    signal_end = config.set_param(user_config, 'FK', 'signal_end', signal_end, 'string')
    noise_start = config.set_param(user_config, 'FK', 'noise_start', noise_start, 'string')
    noise_end = config.set_param(user_config, 'FK', 'noise_end', noise_end, 'string')
    fk_window_len = config.set_param(user_config, 'FK', 'window_len', fk_window_len, 'float')
    fk_sub_window_len = config.set_param(user_config, 'FK', 'sub_window_len', fk_sub_window_len, 'float')
    fk_window_step = config.set_param(user_config, 'FK', 'window_step', fk_window_step, 'float')
    multithread = config.set_param(user_config, 'FK', 'multithread', multithread, 'bool')
    cpu_cnt = config.set_param(user_config, 'FK', 'cpu_cnt', cpu_cnt, 'int')
    write_wvfrms = config.set_param(user_config, 'FK', 'write_wvfrms', write_wvfrms, 'bool')

    fd_window_len = config.set_param(user_config, 'FD', 'window_len', fd_window_len, 'float')
    p_value = config.set_param(user_config, 'FD', 'p_value', p_value, 'float')
    min_duration = config.set_param(user_config, 'FD', 'min_duration', min_duration, 'float')
    back_az_width = config.set_param(user_config, 'FD', 'back_az_width', back_az_width, 'float')
    fixed_thresh = config.set_param(user_config, 'FD', 'fixed_thresh', fixed_thresh, 'float')
    return_thresh = config.set_param(user_config, 'FD', 'return_thresh', return_thresh, 'bool')

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
    click.echo("  window_len (fk): " + str(fk_window_len))
    click.echo("  sub_window_len (fk): " + str(fk_sub_window_len))
    click.echo("  window_step (fk): " + str(fk_window_step))
    click.echo("  multithread: " + str(multithread))
    if multithread or cpu_cnt is not None:
        click.echo("  cpu_cnt: " + str(cpu_cnt))
        pl = Pool(cpu_cnt)
    else:
        pl = None
    click.echo("  write_wvfrms: " + str(write_wvfrms))

    click.echo(" ")
    click.echo("  window_len (fd): " + str(fd_window_len))
    click.echo("  p_value: " + str(p_value))
    click.echo("  min_duration: " + str(min_duration))
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  fixed_thresh: " + str(fixed_thresh))
    click.echo("  return_thresh: " + str(return_thresh))

    # Check data option and populate obspy Stream
    if db_url is not None:
        db_info = {'url': db_url, 'site': db_site, 'wfdisc': db_wfdisc}
    else:
        db_info = None
    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    if write_wvfrms:
        if local_wvfrms is None:
            click.echo('\n' + "Writing waveform data to local SAC files...")
            data_io.write_stream_to_sac(stream, latlon)
        else: 
            click.echo('\n' + "Cannot write waveform data when using local data...")

    if latlon:
        array_loc = latlon[0]
    else:
        array_loc = [stream[0].stats.sac['stla'], stream[0].stats.sac['stlo']]

    if local_wvfrms is not None and "/" in local_wvfrms:
        output_id = os.path.dirname(local_wvfrms) + "/"
    else:
        output_id = ""

    output_id = output_id + tr.stats.network + "." + os.path.commonprefix([tr.stats.station for tr in stream])
    output_id = output_id + '_' + "%02d" % tr.stats.starttime.year + ".%02d" % tr.stats.starttime.month + ".%02d" % tr.stats.starttime.day
    output_id = output_id + '_' + "%02d" % tr.stats.starttime.hour + "." + "%02d" % tr.stats.starttime.minute + "." + "%02d" % tr.stats.starttime.second
    output_id = output_id + '-' + "%02d" % tr.stats.endtime.hour + "." + "%02d" % tr.stats.endtime.minute + "." + "%02d" % tr.stats.endtime.second

    # Check if using a signal window
    if signal_start is not None:
        t1 = UTCDateTime(signal_start)
        t2 = UTCDateTime(signal_end)

        click.echo('\n' + "Trimming data to signal analysis window...")
        click.echo('\t' + "start time: " + str(t1))
        click.echo('\t' + "end time: " + str(t2))


        warning_message = "signal_start and signal_end values poorly defined."
        if t1 > t2:
            warning_message = warning_message + "  signal_start after signal_end."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        elif t1 < stream[0].stats.starttime:
            warning_message = warning_message + "  signal_start before data start time."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        elif t2 > stream[0].stats.endtime:
            warning_message = warning_message + "  signal_end after data end time."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        else:
            stream.trim(t1, t2)


    # Define DOA values
    back_az_vals = np.arange(back_az_min, back_az_max, back_az_step)
    trc_vel_vals = np.arange(trace_vel_min, trace_vel_max, trace_vel_step)

    # run fk analysis
    beam_times, beam_peaks = fkd.run_fk(stream, latlon, [freq_min, freq_max], fk_window_len, fk_sub_window_len, fk_window_step, method, back_az_vals, trc_vel_vals, pl)

    print("Running adaptive f-detector..." + '\n')
    TB_prod = (freq_max - freq_min) * fk_window_len
    min_seq = max(2, int(min_duration / fk_window_len))
    dets, thresh_vals = fkd.run_fd(beam_times, beam_peaks, fd_window_len, TB_prod, len(stream), p_value, min_seq, back_az_width, fixed_thresh, True)

    if local_fk_label is None or local_fk_label == "auto":
        local_fk_label = output_id

    # save fk results
    dt = np.array([(tn - np.datetime64(tr.stats.starttime)).astype('m8[ms]').astype(float) * 1.0e-3 for tn in beam_times])
    fk_results = np.hstack((np.atleast_2d(dt).T, beam_peaks))
    fk_header = data_io.fk_header(stream, latlon, freq_min, freq_max, back_az_min, back_az_max, back_az_step, trace_vel_min, trace_vel_max, trace_vel_step, method, 
        signal_start, signal_end, noise_start, noise_end, fk_window_len, fk_sub_window_len, fk_window_step)

    if not os.path.isfile(local_fk_label + ".fk_results.dat"):
        click.echo('\n' + "Writing results into " + local_fk_label + ".fk_results.dat")
        np.savetxt(local_fk_label + ".fk_results.dat", fk_results, header=fk_header)
    else:
        k = 0
        while os.path.isfile(local_fk_label + "-v" + str(k) + ".fk_results.dat"):
            k += 1
        click.echo('\n' + "WARNING!  fk results file(s) already exist." + '\n' + "Writing a new version: " + local_fk_label + "-v" + str(k) + ".fk_results.dat")
        np.savetxt(local_fk_label + "-v" + str(k) + ".fk_results.dat", fk_results, header=fk_header)

    # save detection results
    det_list = []
    for det_info in dets:
        det_list = det_list + [data_io.define_detection(det_info, array_loc, len(stream), [freq_min, freq_max], note="InfraPy CLI detection")]

    if local_detect_label is None or local_detect_label == "auto":
        local_detect_label = output_id

    click.echo("Writing detection results using label: " + local_detect_label)
    data_io.detection_list_to_json(local_detect_label + ".dets.json", det_list)

    if return_thresh:
        np.save(local_detect_label + ".thresholds", thresh_vals)

    if pl is not None:
        pl.terminate()
        pl.close()





