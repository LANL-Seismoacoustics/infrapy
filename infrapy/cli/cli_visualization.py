#!/usr/bin/env python
import sys
import os
import warnings 

import click
import json
import configparser as cnfg
from matplotlib.pyplot import figure
import numpy as np

from obspy import UTCDateTime

from multiprocessing import Pool

from infrapy.location import bisl

from infrapy.utils import config
from infrapy.utils import data_io
from infrapy.detection import visualization as det_vis
from infrapy.location import visualization as loc_vis


@click.command('fk', short_help="Visualize beamforming (fk) results")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--local-fk-label", help="Local beamforming (fk) data files", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def fk(config_file, local_wvfrms, local_latlon, fdsn, db_config, network, station, location, 
    channel, starttime, endtime, freq_min, freq_max, local_fk_label, figure_out, show_figure):
    '''
    Visualize beamforming (fk) results

    \b
    Example usage (run from infrapy/examples directory after running the run_fk examples):
    \tinfrapy plot fk --local-wvfrms 'data/YJ.BRP*'
    \tinfrapy plot fk --config-file config/detection_local.config
    \tinfrapy plot fk --config-file config/detection_fdsn.config --figure-out FDSN_fk-results.png --show-figure False

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##         Beamforming (fk)        ##")
    click.echo("##          Visualization          ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Database configuration and info   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

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

    # Frequency limits
    freq_min = config.set_param(user_config, 'FK', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'FK', 'freq_max', freq_max, 'float')

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
        
    click.echo("  local_fk_label: " + str(local_fk_label))

    click.echo('\n' + "Visualization parameters:")
    click.echo("  freq_min: " + str(freq_min))
    click.echo("  freq_max: " + str(freq_max))
    if figure_out:
        click.echo("  figure_out: " + figure_out)

    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if waveform data is specified and populate obspy Stream
    if stream is not None:
        stream.filter("bandpass", freqmin=freq_min, freqmax=freq_max)

        if local_fk_label is None or local_fk_label == "auto":
            if local_wvfrms is not None and "/" in local_wvfrms:
                local_fk_label = os.path.dirname(local_wvfrms) + "/"
            else:
                local_fk_label = ""
            local_fk_label = local_fk_label + data_io.stream_label(stream)

        if ".fk_results.dat" not in local_fk_label:
            local_fk_label = local_fk_label + ".fk_results.dat"

        temp = np.loadtxt(local_fk_label)
        dt, beam_peaks = temp[:, 0], temp[:, 1:]

        temp = open(local_fk_label, 'r')
        for line in temp:
            if "t0:" in line:
                t0 = np.datetime64(line.split(' ')[-1][:-1])
            elif "freq_min" in line:
                freq_min = float(line.split(' ')[-1])
            elif "freq_max" in line:
                freq_max = float(line.split(' ')[-1])

        beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
        det_vis.plot_fk1(stream, latlon, beam_times, beam_peaks, title=local_fk_label, output_path=figure_out, show_fig=show_figure)
    else:
        if os.path.isfile(local_fk_label):
            temp = np.loadtxt(local_fk_label)
            dt, beam_peaks = temp[:, 0], temp[:, 1:]

            temp = open(local_fk_label, 'r')
            for line in temp:
                if "t0:" in line:
                    t0 = np.datetime64(line.split(' ')[-1][:-1])
                elif "freq_min" in line:
                    freq_min = float(line.split(' ')[-1])
                elif "freq_max" in line:
                    freq_max = float(line.split(' ')[-1])

            beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
            det_vis.plot_fk2(beam_times, beam_peaks, output_path=figure_out, show_fig=show_figure)
        else:
            msg = "Beamforming (fk) results not found.  No file: " + local_fk_label
            warnings.warn(msg)


@click.command('fd', short_help="Visualize detections from beamforming results")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--local-fk-label", help="Local beamforming (fk) data files", default=None)
@click.option("--local-detect-label", help="Local detection data files", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def fd(config_file, local_wvfrms, local_latlon, fdsn, db_config, network, station, location, channel, starttime, endtime,
    local_fk_label, local_detect_label, figure_out, show_figure):
    '''
    Visualize detection (fd) results

    \b
    Example usage (run from infrapy/examples directory after running fd examples or fkd examples):
    \tinfrapy plot fd --local-wvfrms 'data/YJ.BRP*'
    \tinfrapy plot fd --config-file config/detection_local.config
    \tinfrapy plot fd --config-file config/detection_fdsn.config --figure-out FDSN_fd-results.png --show-figure False

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##          Detection (fd)         ##")
    click.echo("##          Visualization          ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Database configuration and info   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

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
        
    click.echo("  local_fk_label: " + str(local_fk_label))
    click.echo("  local_detect_label: " + str(local_detect_label))

    if figure_out:
        click.echo("  figure_out: " + figure_out)

    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if waveform data is specified and populate obspy Stream
    if stream is not None:
        if local_fk_label is None or local_fk_label == "auto":
            if local_wvfrms is not None and "/" in local_wvfrms:
                local_fk_label = os.path.dirname(local_wvfrms) + "/"
            else:
                local_fk_label = ""
            local_fk_label = local_fk_label + data_io.stream_label(stream)
            
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

        beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
        stream.filter("bandpass", freqmin=freq_min, freqmax=freq_max)
        
        # Read in detection list
        if local_detect_label is None or local_detect_label == 'auto':
            local_detect_label = local_fk_label

        det_list = data_io.set_det_list(local_detect_label + ".dets.json", merge=True)
        if len(det_list) == 0:
            click.echo("Note: no detections found in analysis.")

        if os.path.isfile(local_detect_label + ".fd_thresholds.dat"):
            temp = np.loadtxt(local_detect_label + ".fd_thresholds.dat")
            thresh_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
            det_thresh = [thresh_times, temp[:, 1]]

        else:
            det_thresh = None

        det_vis.plot_fk1(stream, latlon, beam_times, beam_peaks, detections=det_list, title=local_fk_label, output_path=figure_out, det_thresh=det_thresh, show_fig=show_figure)
    else:
        if os.path.isfile(local_fk_label + ".fk_times.npy"):
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

            beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
            stream.filter("bandpass", freqmin=freq_min, freqmax=freq_max)

            # Read in detection list
            det_list = data_io.set_det_list(local_detect_label, merge=True)
            if len(det_list) == 0:
                click.echo("Note: no detections found in analysis.")

            det_vis.plot_fk2(beam_times, beam_peaks, detections=det_list, output_path=figure_out, show_fig=show_figure)
        else:
            msg = "Beamforming (fk) results not found.  No file: " + local_fk_label + ".fk_times.npy"
            warnings.warn(msg)


@click.command('sd', short_help="Visualize detection(s) from spectral analysis")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--local-detect-label", help="Local detection data files", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['SD']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['SD']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--signal-start", help="Start of analysis window", default=None)
@click.option("--signal-end", help="End of analysis window", default=None)
@click.option("--single-det-index", help="Index of a single detection", default=None, type=int)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def sd(config_file, local_wvfrms, local_latlon, fdsn, db_config, network, station, location, channel, starttime, endtime,
    local_detect_label, freq_min, freq_max, signal_start, signal_end, single_det_index, figure_out, show_figure):
    '''
    Visualize spectral detection (sd) results

    \b
    Example usage (run from infrapy/examples directory after running fd examples or fkd examples):
    \tinfrapy plot sd --local-wvfrms 'data/YJ.BRP1..EDF.SAC'

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##     Spectral Detection (sd)     ##")
    click.echo("##          Visualization          ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Database configuration and info   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

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
        
    click.echo("  local_detect_label: " + str(local_detect_label))

    if figure_out:
        click.echo("  figure_out: " + figure_out)

    # Algorithm parameters
    freq_min = config.set_param(user_config, 'SD', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'SD', 'freq_max', freq_max, 'float')
    signal_start = config.set_param(user_config, 'SD', 'signal_start', signal_start, 'string')
    signal_end = config.set_param(user_config, 'SD', 'signal_end', signal_end, 'string')

    stream, _ = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

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

    if local_detect_label is not None:
        if ".dets.json" not in local_detect_label:
            local_detect_label = local_detect_label + ".dets.json"
        det_list = json.load(open(local_detect_label))
    else:
        if local_wvfrms is not None and "/" in local_wvfrms:
            output_id = os.path.dirname(local_wvfrms) + "/"
        else:
            output_id = ""
        output_id = output_id + data_io.stream_label(stream)

        if ".dets.json" not in output_id:
            output_id = output_id + ".dets.json"
        det_list = json.load(open(output_id))

    if len(det_list) == 0:
        click.echo("Note: no detections found in analysis.")

    if single_det_index is not None:
        if single_det_index <= len(det_list) - 1:
            click.echo("Plotting detection info for detection index (" + str(single_det_index) + ")..." + '\n')
            det_vis.plot_sd_single(stream[0], det_list[single_det_index], [freq_min, freq_max], output_path=figure_out, show_fig=show_figure)       
        else:
            click.echo("Invalid detection index (" + str(single_det_index) + "), only " + str(len(det_list)) + " detections in file.")
    else:
        click.echo("Plotting spectrogram with detection info..." + '\n')
        det_vis.plot_sd(stream[0], det_list, [freq_min, freq_max], output_path=figure_out, show_fig=show_figure)


@click.command('dets', short_help="Plot detections on a map")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-detect-label", help="Detection path and pattern", default=None)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
def dets(config_file, range_max, local_detect_label, figure_out, offline_maps_dir):
    '''
    Visualize detections on a map

    \b
    Example usage (run from infrapy/examples directory after running run_assoc example):
    \tinfrapy plot dets --local-detect-label 'data/Blom_etal2020_GJI/*'
    \tinfrapy plot dets --local-detect-label 'GJI_example-ev0.dets.json'  --range-max 1000

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##          Detection List         ##")
    click.echo("##             Mapping             ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    


    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    local_detect_label = config.set_param(user_config, 'DETECTION IO', 'local_detect_label', local_detect_label, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  local_detect_label: " + str(local_detect_label))

    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    offline_maps_dir = config.set_param(user_config, 'VISUALIZATION', 'offline_maps_dir', offline_maps_dir, 'string')

    click.echo('\n' + "Visualization parameters:")
    click.echo("  range_max: " + str(range_max) + '\n')
    if offline_maps_dir:
        click.echo("  offline maps directory: {}".format(offline_maps_dir))
        loc_vis.use_offline_maps(offline_maps_dir)

    det_list = data_io.set_det_list(local_detect_label, merge=True)

    click.echo('\n' + "Drawing map with detection back azimuth projections...")
    loc_vis.plot_dets_on_map(det_list, range_max=range_max, output_path=figure_out)


@click.command('loc', short_help="Plot localization result on a map")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-detect-label", help="Detection path and pattern", default=None)
@click.option("--local-loc-label", help="Localization results path", default=None)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--zoom", help="Option to zoom in on the estimated source region", default=False)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--grnd-truth", help="Ground truth location", default=None)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
def loc(config_file, local_detect_label, local_loc_label, range_max, zoom, figure_out, grnd_truth, offline_maps_dir):
    '''
    Visualize BISL results in with wide or zoomed format

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot loc --local-detect-label GJI_example-ev0 --local-loc-label GJI_example-ev0 --range-max 1200.0
    \tinfrapy plot loc --local-detect-label GJI_example-ev0 --local-loc-label GJI_example-ev0 --zoom true

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##       Localization Mapping      ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")  

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    local_detect_label = config.set_param(user_config, 'DETECTION IO', 'local_detect_label', local_detect_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  local_event_label: " + str(local_detect_label))
    click.echo("  local_loc_label: " + str(local_loc_label))

    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    offline_maps_dir = config.set_param(user_config, 'VISUALIZATION', 'offline_maps_dir', offline_maps_dir, 'string')

    click.echo('\n' + "Visualization parameters:")
    click.echo("  range_max: " + str(range_max))
    click.echo("  zoom: " + str(zoom))

    if grnd_truth is not None:
        grnd_truth = [float(val) for val in grnd_truth.strip(' ()[]').split(',')]

    if offline_maps_dir:
        click.echo("  offline maps directory: {}".format(offline_maps_dir))
        loc_vis.use_offline_maps(offline_maps_dir)

    click.echo('\n' + "Reading in detection list...")
    det_list = data_io.set_det_list(local_detect_label, merge=False)
    if ".loc.json" in local_loc_label:
        bisl_result = json.load(open(local_loc_label))
    else:
        bisl_result = json.load(open(local_loc_label + ".loc.json"))

    click.echo('\n' + "BISL Summary:")
    click.echo(bisl.summarize(bisl_result))

    click.echo("Drawing map with BISL source location estimate...")
    loc_vis.plot_loc(det_list, bisl_result, range_max=range_max, zoom=zoom, title=None, output_path=figure_out, grnd_truth=grnd_truth)
    


@click.command('origin-time', short_help="Plot origin time distribution")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-loc-label", help="Localization results", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
def origin_time(config_file, local_loc_label, figure_out):
    '''
    Visualize the BISL origin time distribution

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot origin-time --local-loc-label GJI_example-ev0
    '''
    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##        Origin Time Plot         ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")  

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_event_label', local_loc_label, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  local_detect_label: " + str(local_loc_label))

    click.echo('\n' + "Reading in BISL results...")
    if ".loc.json" in local_loc_label:
        bisl_result = json.load(open(local_loc_label))
    else:
        bisl_result = json.load(open(local_loc_label + ".loc.json"))

    click.echo("Plotting origin time distribution...")
    loc_vis.plot_origin_time(bisl_result, output_path=figure_out)
    

@click.command('yield', short_help="Plot yield estimate distribution")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-yld-label", help="Yield estimate result", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
def yield_plot(config_file, local_yld_label, figure_out):

    '''
    Visualize the SpYE result

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot yield --local-yld-label HRR-5.yld.json

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##       Yield Estimate Plot       ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")  

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    local_yld_label = config.set_param(user_config, 'DETECTION IO', 'local_yld_label', local_yld_label, 'string')
    figure_out = config.set_param(user_config, 'DETECTION IO', 'figure_out', figure_out, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  local_yld_label: " + str(local_yld_label))

    click.echo('\n' + "Reading in SpYE results...")
    if ".yld.json" in local_yld_label:
        spye_result = json.load(open(local_yld_label))
    else:
        spye_result = json.load(open(local_yld_label + ".yld.json"))

    click.echo('\n' + 'Results Summary (tons eq. TNT):')
    click.echo('\t' + "Maximum a Posteriori Yield: " + str(spye_result['yld_vals'][np.argmax(spye_result['yld_pdf'])]))
    click.echo('\t' + "68% Confidence Bounds: " + str(spye_result['conf_bnds'][0]))
    click.echo('\t' + "95% Confidence Bounds: " + str(spye_result['conf_bnds'][1]))

    click.echo('\n' + "Plotting yield PDF...")
    loc_vis.plot_spye(spye_result, output_path=figure_out)

    click.echo("")
