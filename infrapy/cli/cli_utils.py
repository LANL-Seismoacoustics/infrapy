#!which python
"""
cli_utils.py

Utility methods accessible in the command line interface (CLI) of infrapy

Author: pblom@lanl.gov    
"""

from pickletools import read_long1
import click

import numpy as np

from obspy import UTCDateTime 

from pyproj import Geod

from infrapy.propagation import likelihoods as lklhds
from infrapy.utils import data_io

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
    click.echo("##           InfraPy           ##")
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
@click.option("--rcvr-lat", help="Source latitude", default=None, prompt="Enter receiver latitude: ")
@click.option("--rcvr-lon", help="Source longitude", default=None, prompt="Enter receiver longitude: ")
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
    click.echo("##           InfraPy           ##")
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
    click.echo("  Propagation Range: " + str(np.round(rng,2)) + " km")
    click.echo("  Propagation Azimuth: " + str(np.round(az, 2)) + " degrees")
    click.echo("  Estimted arrival time range:")
    click.echo("    " + str(UTCDateTime(src_time) + np.round(rng / celerity_max, 0))[:-8])
    click.echo("    " + str(UTCDateTime(src_time) + np.round(rng / celerity_min, 0))[:-8] + '\n')








