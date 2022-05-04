#!which python
"""
cli_utils.py

Utility methods accessible in the command line interface (CLI) of infrapy

Author: pblom@lanl.gov    
"""

import configparser as cnfg

from pickletools import read_long1
import click

import numpy as np

from obspy import UTCDateTime 

from pyproj import Geod

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
@click.option("--rcvr-lat", help="Receiver latitude", default=None, prompt="Enter receiver latitude: ")
@click.option("--rcvr-lon", help="Receiver longitude", default=None, prompt="Enter receiver longitude: ")
@click.option("--celerity-min", help="Minimum celerity", default=0.24)
@click.option("--celerity-max", help="Maximum celerity", default=0.35)
def arrival_time(src_lat, src_lon, src_time, rcvr_lat, rcvr_lon, celerity_min, celerity_max):
    '''
    Compute the range of possible arrivals times for a source-receiver pair given a range of celerity values
    
    \b
    Example usage (requires InfraGA/GeoAc arrival output):
    \tinfrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr-lat 40.0 --rcvr-lon -110.0

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
    click.echo("  Receiver Location: (" + str(rcvr_lat) + ", " + str(rcvr_lon) + ")")

    click.echo("")
    click.echo("  Celerity Range: (" + str(celerity_min) + ", " + str(celerity_max) + ")")

    sph_proj = Geod(ellps='sphere')
    temp = sph_proj.inv(src_lon, src_lat, rcvr_lon, rcvr_lat, radians=False)
    az = temp[0]
    rng = temp[2] / 1000.0

    click.echo("")
    click.echo("  Propagation range: " + str(np.round(rng,2)) + " km")
    click.echo("  Propagation azimuth: " + str(np.round(az, 2)) + " degrees")
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



@click.command('check_db_wvfrms', short_help="Check waveform pull from database")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--db-url", help="Database URL for waveform data files", default=None)
@click.option("--db-site", help="Database site table for waveform data files", default=None)
@click.option("--db-wfdisc", help="Database wfdisc table for waveform data files", default=None)

@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)

@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
def check_db_wvfrm(config_file, db_url, db_site, db_wfdisc, network, station, location, channel, starttime, endtime):
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
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Database and data IO parameters   
    db_url = config.set_param(user_config, 'WAVEFORM IO', 'db_url', db_url, 'string')
    db_site = config.set_param(user_config, 'WAVEFORM IO', 'db_site', db_site, 'string')
    db_wfdisc = config.set_param(user_config, 'WAVEFORM IO', 'db_wfdisc', db_wfdisc, 'string')

    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    click.echo('\n' + "Data parameters:")
    click.echo("  db_url: " + str(db_url))
    click.echo("  db_site: " + str(db_site))
    click.echo("  db_wfdisc: " + str(db_wfdisc))
    click.echo("  network: " + str(network))
    click.echo("  station: " + str(station))
    click.echo("  location: " + str(location))
    click.echo("  channel: " + str(channel))
    click.echo("  starttime: " + str(starttime))
    click.echo("  endtime: " + str(endtime))


    # Check data option and populate obspy Stream
    db_info = {'url': db_url, 'site': db_site, 'wfdisc': db_wfdisc}
    stream, latlon = data_io.set_stream(None, None, db_info, network, station, location, channel, starttime, endtime, None)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    for line in latlon:
        print(line[0], '\t', line[1])



