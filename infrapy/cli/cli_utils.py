#!which python
"""
cli_utils.py

Utility methods accessible in the command line interface (CLI) of infrapy

Author: pblom@lanl.gov    
"""

import os
import pickle
import imp 

from threading import local 

import warnings

import configparser as cnfg
from inspect import trace

import matplotlib.pyplot as plt 

from pickletools import read_long1
import click

import numpy as np

from obspy import UTCDateTime 

from pyproj import Geod

from infrapy.detection import beamforming_new
from infrapy.propagation import likelihoods as lklhds
from infrapy.utils import config, data_io

@click.command('arrivals2json', short_help="Convert infraGA/GeoAc arrivals to detection file")
@click.option("--arrivals-file", help="InfraGA/GeoAc arrivals file", default=None)
@click.option("--json-file", help="JSON format detection file", default=None)
@click.option("--grnd-snd-spd", help="Ground sound speed", default=340.0)
@click.option("--src-time", help="Source time", default="2020-01-01T00:00:00")
@click.option("--peakf-value", help="Fixed F-value", default=25.0)
@click.option("--array-dim", help="Array dimension", default=6)
def arrivals2json(arrivals_file, json_file, grnd_snd_spd, src_time, peakf_value, array_dim):
    '''
    Convert infraGA/GeoAc eigenray arrival results into a json detection list usable in InfraPy
    
    \b
    Example usage (requires InfraGA/GeoAc arrival output):
    \tinfrapy arrivals2json --arrivals-file example.arrivals.dat --json-file example.dets.json --grnd-snd-spd 335.0 --src-time "2020-12-25T00:00:00"

    '''
    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##        arrivals2json        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  
    
    
    click.echo("")
    click.echo("  arrivals_file: " + str(arrivals_file))
    click.echo("  json_file: " + str(json_file))

    click.echo("")
    click.echo("  grnd_snd_spd: " + str(grnd_snd_spd))
    click.echo("  src_time: " + str(src_time))
    click.echo("  peakF_value: " + str(peakf_value))
    click.echo("  array_dim: " + str(array_dim))

    arrivals = np.loadtxt(arrivals_file)

    det_list = []
    for line in arrivals:
        det = lklhds.InfrasoundDetection(lat_loc=np.round(line[3], 3), lon_loc=np.round(line[4], 3), time=(UTCDateTime(src_time) + line[5]), azimuth=np.round(line[9], 2), f_stat=peakf_value, array_d=array_dim)
        det.trace_velocity = np.round(grnd_snd_spd / np.cos(np.radians(line[8])), 1)
        det.note = "InfraGA/GeoAc arrival output"
        det_list = det_list + [det]

    data_io.detection_list_to_json(json_file, det_list)


@click.command('arrival-time', short_help="Estimate the arrival time for a source-receiver pair")
@click.option("--src-lat", help="Source latitude", default=None, prompt="Enter source latitude: ")
@click.option("--src-lon", help="Source longitude", default=None, prompt="Enter source longitude: ")
@click.option("--src-time", help="Source time", default=None, prompt="Enter source time: ")
@click.option("--rcvr-lat", help="Receiver latitude", default=None)
@click.option("--rcvr-lon", help="Receiver longitude", default=None)
@click.option("--rcvr", help="Reference IMS station (e.g., 'I53')", default=None)
@click.option("--celerity-min", help="Minimum celerity", default=0.24)
@click.option("--celerity-max", help="Maximum celerity", default=0.35)
def arrival_time(src_lat, src_lon, src_time, rcvr_lat, rcvr_lon, rcvr, celerity_min, celerity_max):
    '''
    Compute the range of possible arrivals times for a source-receiver pair given a range of celerity values.
    Can use a receiver latitude/longitude or reference from a list (currently only IMS stations)
    
    \b
    Example usage (requires InfraGA/GeoAc arrival output):
    \tinfrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr-lat 40.0 --rcvr-lon -110.0
    \tinfrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr I57US

    '''
    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         arrival-time        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  
    
    click.echo("  Source Time: " + src_time)
    click.echo("  Source Location: (" + str(src_lat) + ", " + str(src_lon) + ")")
    if rcvr is not None:
        click.echo('\n' + "  User specified reference receiver: " + str(rcvr))

        IMS_info = pickle.load(open(imp.find_module('infrapy')[1] + '/resources/IMS_infrasound_locs.pkl', 'rb'), encoding='latin1')
        for line in IMS_info:
            if rcvr in line[0]:
                click.echo("  Reference IMS station match: " + line[0])
                rcvr_lat, rcvr_lon = line[1][:2]
                click.echo("  Receiver Location: (" + str(rcvr_lat) + ", " + str(rcvr_lon) + ")")
                break
        if rcvr_lat is None:
            warning_message = "Specified reference receiver (" + rcvr + ") not found in IMS info."
            warnings.warn((warning_message))
            return 0

    elif rcvr_lat is not None and rcvr_lon is not None:
        click.echo("  Receiver Location: (" + str(rcvr_lat) + ", " + str(rcvr_lon) + ")")
    else:
        warning_message = "Method requires either a reference receiver or user defined latitude and longitude"
        warnings.warn((warning_message))
        return 0

    click.echo("")
    click.echo("  Celerity Range: (" + str(celerity_min) + ", " + str(celerity_max) + ")")

    sph_proj = Geod(ellps='sphere')
    temp = sph_proj.inv(src_lon, src_lat, rcvr_lon, rcvr_lat, radians=False)
    az, back_az = temp[0], temp[1]
    rng = temp[2] / 1000.0

    if back_az > 180.0:
        back_az = back_az - 360.0
    elif back_az < -180.0:
        back_az = back_az + 360.0

    click.echo("")
    click.echo("  Propagation range: " + str(np.round(rng,2)) + " km")
    click.echo("  Propagation azimuth: " + str(np.round(az, 2)) + " degrees" + '\n')

    click.echo("  Estimated arrival back azimuth: " + str(np.round(back_az, 2)) + " degrees")
    click.echo("  Estimated arrival time range:")
    click.echo("    " + str(UTCDateTime(src_time) + np.round(rng / celerity_max, 0))[:-8])
    click.echo("    " + str(UTCDateTime(src_time) + np.round(rng / celerity_min, 0))[:-8] + '\n')



@click.command('calc-celerity', short_help="Compute the celerity for an arrival from a known source")
@click.option("--src-lat", help="Source latitude", default=None, prompt="Enter source latitude: ")
@click.option("--src-lon", help="Source longitude", default=None, prompt="Enter source longitude: ")
@click.option("--src-time", help="Source time", default=None, prompt="Enter source time: ")
@click.option("--arrival-lat", help="Arrival latitude", default=None, prompt="Enter arrival latitude: ")
@click.option("--arrival-lon", help="Arrival longitude", default=None, prompt="Enter arrival longitude: ")
@click.option("--arrival-time", help="Arrival time", default=None, prompt="Enter arrival time: ")
def calc_celerity(src_lat, src_lon, src_time, arrival_lat, arrival_lon, arrival_time):
    '''
    Compute the range of possible arrivals times for a source-receiver pair given a range of celerity values
    
    \b
    Example usage (requires InfraGA/GeoAc arrival output):
    \tinfrapy utils calc-celerity --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --arrival-lat 40.0 --arrival-lon -110.0 --arrival-time "2020-12-25T01:03:50"

    '''
    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##        calc-celerity        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  
    
    click.echo("  Source Time: " + src_time)
    click.echo("  Source Location: (" + str(src_lat) + ", " + str(src_lon) + ")")

    click.echo('\n' + "  Arrival Time: " + arrival_time)
    click.echo("  Arrival Location: (" + str(arrival_lat) + ", " + str(arrival_lon) + ")")

    dt = UTCDateTime(arrival_time) - UTCDateTime(src_time)

    sph_proj = Geod(ellps='sphere')
    temp = sph_proj.inv(src_lon, src_lat, arrival_lon, arrival_lat, radians=False)
    az = temp[0]
    rng = temp[2]

    click.echo("")
    click.echo("  Propagation time: " + str(np.round(dt, 2)) + " s")
    click.echo("  Propagation range: " + str(np.round(rng / 1000.0, 2)) + " km")
    click.echo("  Propagation azimuth: " + str(np.round(az, 2)) + " degrees")
    click.echo("  Arrival celerity: " + str(np.round(rng / dt, 1)) + " m/s" + '\n')



@click.command('check-db-wvfrms', short_help="Check waveform pull from database")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--db-config", help="Database configuration file", default=None)

@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)

@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
def check_db_wvfrm(config_file, db_config, network, station, location, channel, starttime, endtime):
    '''
    Test database pull of waveform data for beamforming (fk or fdk) analysis

    \b
    Example usage (detection_db.config will be unique to your database pull):
    \tinfrapy run_fk --config-file config/detection_db.config

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##       check_db_wvfrms       ##")
    click.echo("##                             ##")
    click.echo("#################################")
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

    # Database and data IO parameters   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    click.echo('\n' + "Data parameters:")
    click.echo("  db_config: " + str(db_config))
    click.echo("  network: " + str(network))
    click.echo("  station: " + str(station))
    click.echo("  location: " + str(location))
    click.echo("  channel: " + str(channel))
    click.echo("  starttime: " + str(starttime))
    click.echo("  endtime: " + str(endtime))

    # Check data option and populate obspy Stream
    db_info = cnfg.ConfigParser()
    db_info.read(db_config)

    stream, latlon = data_io.set_stream(None, None, db_info, network, station, location, channel, starttime, endtime, None)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    click.echo('\nLocation info:')    
    for line in latlon:
        click.echo(str(line[0]) + '\t' +  str(line[1]))


@click.command('write-wvfrms', short_help="Save waveforms from FDSN or database")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)

@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)

@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
def write_wvfrms(config_file, db_config, fdsn, network, station, location, channel, starttime, endtime):
    '''
    Write waveform data from an FDSN or database pull into local SAC files

    \b
    Example usage (detection_db.config will be unique to your database pull):
    \tinfrapy utils write-wvfrms --config-file config/detection_fdsn.config

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         write-wvfrms        ##")
    click.echo("##                             ##")
    click.echo("#################################")
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

    # Database and data IO parameters   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

    # FDSN waveform IO parameters
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    click.echo('\n' + "Data parameters:")
    if fdsn is not None:
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

        click.echo("  db_url: " + str(db_config))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    else:
        click.echo("Invalid data parameters.  Requires fdsn or db info.")

    stream, latlon = data_io.set_stream(None, fdsn, db_info, network, station, location, channel, starttime, endtime, None)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    click.echo('\n' + "Writing waveform data to local SAC files...")
    data_io.write_stream_to_sac(stream, latlon)


@click.command('best-beam', short_help="Compute the best beam via shift/stack")
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
@click.option("--back-az", help="Back azimuth of user specified beam (degrees)", default=None, type=float)
@click.option("--trace-vel", help="Trace velocity of user specified beam (m/s))", default=None, type=float)
@click.option("--signal-start", help="Start of signal window", default=None)
@click.option("--signal-end", help="End of signal window", default=None)
@click.option("--hold-figure", help="Hold figure open", default=True)
def best_beam(config_file, local_wvfrms, fdsn, db_url, db_site, db_wfdisc, local_latlon, network, station, location, channel, starttime, endtime, local_fk_label, freq_min, freq_max,
    back_az, trace_vel, signal_start, signal_end, hold_figure):
    '''
    Shift and stack the array data to compute the best beam.  Can be run adaptively using the fk_results.dat file or along a specific beam.

    \b
    Example usage (requires 'infrapy run_fk --config-file config/detection_local.config' run first):
    \tinfrapy utils best-beam --config-file config/detection_local.config
    \tinfrapy utils best-beam --config-file config/detection_local.config --back-az -39.0 --trace-vel 358.0
    \tinfrapy utils best-beam --config-file config/detection_local.config --signal-start '2012-04-09T18:13:00' --signal-end '2012-04-09T18:15:00'

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##          best-beam          ##")
    click.echo("##                             ##")
    click.echo("#################################")
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

    # Local fk file
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
        click.echo("Invalid data parameters.  Requires fdsn or db info.")

    if local_fk_label is not None:
        click.echo("  local_fk_label: " + str(local_fk_label))

    # Algorithm parameters
    freq_min = config.set_param(user_config, 'FK', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'FK', 'freq_max', freq_max, 'float')

    signal_start = config.set_param(user_config, 'FK', 'signal_start', signal_start, 'string')
    signal_end = config.set_param(user_config, 'FK', 'signal_end', signal_end, 'string')

    click.echo('\n' + "Algorithm parameters:")
    click.echo("  freq_min: " + str(freq_min))
    click.echo("  freq_max: " + str(freq_max))
    click.echo("  signal_start: " + str(signal_start))
    click.echo("  signal_end: " + str(signal_end))
    if back_az is not None and trace_vel is not None:
        click.echo("  back_az_step: " + str(back_az))
        click.echo("  trace_vel_min: " + str(trace_vel))

    # Check data option and populate obspy Stream
    if db_url is not None:
        db_info = {'url': db_url, 'site': db_site, 'wfdisc': db_wfdisc}
    else:
        db_info = None
    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    if local_fk_label is None or local_fk_label == "auto":
        local_fk_label = ""
        if local_wvfrms is not None:
            if "/" in local_wvfrms:
                local_fk_label = os.path.dirname(local_wvfrms) + "/"
        
        local_fk_label = local_fk_label + tr.stats.network + "." + os.path.commonprefix([tr.stats.station for tr in stream])
        local_fk_label = local_fk_label + '_' + "%02d" % tr.stats.starttime.year + ".%02d" % tr.stats.starttime.month + ".%02d" % tr.stats.starttime.day
        local_fk_label = local_fk_label + '_' + "%02d" % tr.stats.starttime.hour + "." + "%02d" % tr.stats.starttime.minute + "." + "%02d" % tr.stats.starttime.second
        local_fk_label = local_fk_label + '-' + "%02d" % tr.stats.endtime.hour + "." + "%02d" % tr.stats.endtime.minute + "." + "%02d" % tr.stats.endtime.second
    else:
        if local_fk_label[-15:] == ".fk_results.dat":
            local_fk_label = local_fk_label[:-15]

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

    if back_az is not None and trace_vel is not None:
        click.echo('\n' + "Computing best beam with user specified beam...")
        click.echo('\t' + "Back Azimuth: " + str(back_az))
        click.echo('\t' + "Trace Velocity: " + str(trace_vel))

        stream.filter('bandpass', freqmin=freq_min, freqmax=freq_max)
        x, t, t0, geom = beamforming_new.stream_to_array_data(stream, latlon=latlon)
        X, _, f = beamforming_new.fft_array_data(x, t, fft_window="boxcar")

        sig_est, residual = beamforming_new.extract_signal(X, f, [back_az, trace_vel], geom)
        best_beam = np.fft.irfft(sig_est)[:len(t)] / (t[1] - t[0])
        residuals = np.fft.irfft(residual, axis=1)[:, :len(t)]  / (t[1] - t[0])

    else:
        click.echo('\n' + "Computing adaptive best beam...")
        click.echo('\t' + "fk results file: " + local_fk_label + ".fk_results.dat")

        def _envelope(t0, t1, t2, sigma):
            X1 = np.exp(-(t0 - t1) / sigma)
            X2 = np.exp(-(t0 - t2) / sigma)

            return X2 / ((1.0 + X1) * (1.0 + X2))

        # Read in the fk_results
        fk_t0 = None
        freq_min, freq_max = 0.5, 5.0

        temp = open(local_fk_label + ".fk_results.dat", 'r')
        for line in temp:
            if "t0:" in line:
                fk_t0 = np.datetime64(line.strip('\n').split(' ')[-1][:-1])
            elif "freq_min" in line:
                freq_min = float(line.split(' ')[-1])
            elif "freq_max" in line:
                freq_max = float(line.split(' ')[-1])
        temp.close()

        temp = np.loadtxt(local_fk_label + ".fk_results.dat")
        dt, beam_results = temp[:, 0], temp[:, 1:]
        beam_times = np.array([fk_t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])

        # Filter and extract stream info
        stream.filter('bandpass', freqmin=freq_min, freqmax=freq_max)
        x, t, t0, geom = beamforming_new.stream_to_array_data(stream, latlon=latlon)
        M, _ = x.shape

        best_beam = np.zeros_like(t)
        residuals = np.zeros_like(x)

        window_step = (beam_times[1] - beam_times[0]).astype(float) / 1.0e6
        for n, tn in enumerate((beam_times - t0).astype(float) / 1.0e6):
            if tn >= 0.0 and tn <= (t[-1] - t[0]):
                X, _, f = beamforming_new.fft_array_data(x, t, window=[tn - window_step, tn + window_step], fft_window="boxcar")

                sig_est, residual = beamforming_new.extract_signal(X, f, beam_results[n, :2], geom)
                signal_wvfrm = np.fft.irfft(sig_est) / (t[1] - t[0])
                resid_wvfrms = np.fft.irfft(residual, axis=1) / (t[1] - t[0])

                mask = np.logical_and(tn - window_step <= t, t <= tn + window_step)
                best_beam[mask] = best_beam[mask] + signal_wvfrm[:sum(mask.astype(int))] * _envelope(t[mask], tn - window_step / 2.0, tn + window_step / 2.0, window_step / 20.0)
                for nM in range(M):
                    residuals[nM][mask] = residuals[nM][mask] + resid_wvfrms[nM][:sum(mask.astype(int))] * _envelope(t[mask], tn - window_step / 2.0, tn + window_step / 2.0, window_step / 20.0)

    # add output of waveform data (columns: t : beam : resid_1 : resid_2 : ... : resid_M)
    click.echo('\n' + "Writing results into " + local_fk_label + ".best-beam.dat" + '\n')
    header = "InfraPy Best Beam Results" + '\n'
    header = header + '\n' + "Data summary:" + '\n'
    for tr in stream:
        header = header + "    " + tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime) + '\n'

    header = header + "  t0: " + str(stream[0].stats.starttime) + '\n\n'

    if back_az is not None and trace_vel is not None:
        header = header + "Back Azimuth: " + str(back_az)
        header = header + "Trace Velocity: " + str(trace_vel) + '\n'
    else:
        header = header + "Beamforming (fk) results file: " + local_fk_label + ".fk_results.dat" + '\n'

    header = header + '\n' + "Column summary:" + '\n'
    header = header + "time (rel t0) [s] : beam [Pa] : Resid. 1 [Pa] : Resid. 2[Pa] : ... : Resid. M [Pa]" + '\n'

    output_vals = np.vstack((t, np.vstack((best_beam, residuals))))
    np.savetxt(local_fk_label + ".best-beam.dat", output_vals.T, header=header)
    
    # visualize results with UTC time
    plot_times = np.array([t0 + np.timedelta64(int(tn * 1000.0), 'ms') for tn in t])
    plt.plot(plot_times, best_beam, '-k', linewidth=1.0)
    for nM in range(len(residuals)):
        plt.plot(plot_times, residuals[nM], '-r', linewidth=0.25)

    if hold_figure:
        plt.show()
    else:
        plt.show(block=False)
        plt.pause(5.0)
        plt.close()



