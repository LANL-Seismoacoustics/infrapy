#!/usr/bin/env python
import sys
import os
import warnings 

import click
import configparser as cnfg
from matplotlib.pyplot import figure
import numpy as np

from multiprocessing import Pool

from ..utils import config
from ..utils import data_io

from ..detection import visualization as det_vis


@click.command('plot_fk', short_help="Visualize beamforming (fk) results")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
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
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--local-fk-out", help="Local beamforming (fk) data files", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
def plot_fk(config_file, local_wvfrms, local_latlon, fdsn, db_url, db_site, db_wfdisc, db_origin, network, station, location, 
    channel, starttime, endtime, freq_min, freq_max, local_fk_out, figure_out):
    '''
    Visualize beamforming (fk) results

    \b
    Example usage (run from infrapy/examples directory after running the run_fk examples):
    \tinfrapy plot_fk --local-wvfrms 'data/YJ.BRP*' --figure-out BRP_beam.png
    \tinfrapy plot_fk --config-file config/fk_example2.config

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

    # Database and data IO parameters   
    db_url = config.set_param(user_config, 'database', 'url', db_url, 'string')
    db_site = config.set_param(user_config, 'database', 'site', db_site, 'string')
    db_wfdisc = config.set_param(user_config, 'database', 'wfdisc', db_wfdisc, 'string')
    db_origin = config.set_param(user_config, 'database', 'origin', db_origin, 'string')

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
    local_fk_out = config.set_param(user_config, 'DETECTION IO', 'local_fk_out', local_fk_out, 'string')
    figure_out = config.set_param(user_config, 'VISUALIZATION', 'figure_out', figure_out, 'string')

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
        click.echo("  db_origin: " + str(db_origin))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
        
    click.echo("  local_fk_out: " + str(local_fk_out))
    if figure_out:
        click.echo("  figure_out: " + figure_out)

    # Extract times and peaks from fk results
    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_url, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if waveform data is specified and populate obspy Stream
    if stream is not None:
        stream.filter("bandpass", freqmin=freq_min, freqmax=freq_max)

        if local_fk_out is None or local_fk_out == "auto":
            local_fk_out = stream[-1].stats.network + "." + stream[-1].stats.station
            local_fk_out = local_fk_out + '_' + "%02d" % stream[-1].stats.starttime.hour + "." + "%02d" % stream[-1].stats.starttime.minute + "." + "%02d" % stream[-1].stats.starttime.second
            local_fk_out = local_fk_out + '-' + "%02d" % stream[-1].stats.endtime.hour + "." + "%02d" % stream[-1].stats.endtime.minute + "." + "%02d" % stream[-1].stats.endtime.second
        times = np.load(local_fk_out + ".fk_times.npy")
        peaks = np.load(local_fk_out + ".fk_peaks.npy")

        det_vis.plot_fk1(stream, latlon, times, peaks, title=local_fk_out, output_path=figure_out)
    else:
        if os.path.isfile(local_fk_out + ".fk_times.npy"):
            times = np.load(local_fk_out + ".fk_times.npy")
            peaks = np.load(local_fk_out + ".fk_peaks.npy")
            det_vis.plot_fk2(times, peaks, output_path=figure_out)
        else:
            msg = "Beamforming (fk) results not found.  No file: " + local_fk_out + ".fk_times.npy"
            warnings.warn(msg)


@click.command('plot_fd', short_help="Visualize detections from beamforming results")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
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
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--local-fk-out", help="Local beamforming (fk) data files", default=None)
@click.option("--local-fd-out", help="Local detection data files", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
def plot_fd(config_file, local_wvfrms, local_latlon, fdsn, db_url, db_site, db_wfdisc, db_origin, network, station, location, channel, starttime, endtime,
    freq_min, freq_max, local_fk_out, local_fd_out, figure_out):
    '''
    Visualize detection (fd) results

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot_fd --local-wvfrms 'data/YJ.BRP*' --figure-out BRP_detections.png
    \tinfrapy plot_fd --config-file config/fk_example2.config

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

    # Database and data IO parameters   
    db_url = config.set_param(user_config, 'database', 'url', db_url, 'string')
    db_site = config.set_param(user_config, 'database', 'site', db_site, 'string')
    db_wfdisc = config.set_param(user_config, 'database', 'wfdisc', db_wfdisc, 'string')
    db_origin = config.set_param(user_config, 'database', 'origin', db_origin, 'string')

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
    local_fk_out = config.set_param(user_config, 'DETECTION IO', 'local_fk_out', local_fk_out, 'string')
    local_fd_out = config.set_param(user_config, 'DETECTION IO', 'local_fd_out', local_fd_out, 'string')
    figure_out = config.set_param(user_config, 'VISUALIZATION', 'figure_out', figure_out, 'string')

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
        click.echo("  db_origin: " + str(db_origin))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
        
    click.echo("  local_fk_out: " + str(local_fk_out))
    click.echo("  local_fd_out: " + str(local_fd_out))
    if figure_out:
        click.echo("  figure_out: " + figure_out)

    # Extract times and peaks from fk results
    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_url, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if waveform data is specified and populate obspy Stream
    if stream is not None:
        stream.filter("bandpass", freqmin=freq_min, freqmax=freq_max)

        if local_fk_out is None or local_fk_out == "auto":
            local_fk_out = stream[-1].stats.network + "." + stream[-1].stats.station
            local_fk_out = local_fk_out + '_' + "%02d" % stream[-1].stats.starttime.hour + "." + "%02d" % stream[-1].stats.starttime.minute + "." + "%02d" % stream[-1].stats.starttime.second
            local_fk_out = local_fk_out + '-' + "%02d" % stream[-1].stats.endtime.hour + "." + "%02d" % stream[-1].stats.endtime.minute + "." + "%02d" % stream[-1].stats.endtime.second
        times = np.load(local_fk_out + ".fk_times.npy")
        peaks = np.load(local_fk_out + ".fk_peaks.npy")

        # Read in detection list
        det_list = data_io.set_det_list(local_fk_out + ".dets.json", merge=True)

        det_vis.plot_fk1(stream, latlon, times, peaks, detections=det_list, title=local_fk_out, output_path=figure_out)
    else:
        if os.path.isfile(local_fk_out + ".fk_times.npy"):
            times = np.load(local_fk_out + ".fk_times.npy")
            peaks = np.load(local_fk_out + ".fk_peaks.npy")

            # Read in detection list
            det_list = data_io.set_det_list(local_fd_out, merge=True)

            det_vis.plot_fk2(times, peaks, detections=det_list, output_path=figure_out)
        else:
            msg = "Beamforming (fk) results not found.  No file: " + local_fk_out + ".fk_times.npy"
            warnings.warn(msg)






