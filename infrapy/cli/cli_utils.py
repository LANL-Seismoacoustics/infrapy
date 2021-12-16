#!which python
"""
cli_utils.py

Utility methods accessible in the command line interface (CLI) of infrapy

Author: pblom@lanl.gov    
"""

import click

import numpy as np

from obspy import UTCDateTime 

from ..propagation import likelihoods as lklhds

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
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##          arrivals2json          ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
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

    lklhds.detection_list_to_json(json_file, det_list)




