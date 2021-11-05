#!/usr/bin/env python
import sys

import click
import configparser as cnfg
import numpy as np

from ..location import bisl

from ..utils import config
from ..utils import data_io

@click.command('run_loc', short_help="Estimate source locations and times for events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-dets-in", help="Detection path and pattern", default=None)
@click.option("--back_az_width", help="Width of beam projection (default: " + config.defaults['LOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Maximum source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['LOC']['resolution'] + ")", default=None, type=int)
@click.option("--pgm-model", help="Path geometry model (PGM) file (default: None)", default=None)
def run_loc(config_file, local_dets_in, back_az_width, range_max, resolution, pgm_model):
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
    local_dets_in = config.set_param(user_config, 'DETECTION IO', 'local_dets_in', local_dets_in, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  local_dets_in: " + str(local_dets_in))


    # Algorithm parameters
    back_az_width = config.set_param(user_config, 'LOC', 'back_az_width', back_az_width, 'float')
    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    resolution = config.set_param(user_config, 'LOC', 'resolution', resolution, 'int')
    pgm_model = config.set_param(user_config, 'LOC', 'pgm_model', pgm_model, 'int')

    click.echo('\n' + "Parameter summary:")
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  range_max: " + str(range_max))
    click.echo("  resolution: " + str(resolution))
    click.echo("  pgm_model: " + str(pgm_model))

    events = data_io.set_det_list(local_dets_in, merge=False)
    for det_list in events:
        result = bisl.run(det_list, path_geo_model=pgm_model, resol=resolution, bm_width=back_az_width, rng_max=range_max, rad_min=100.0, rad_max=range_max / 4.0)

        # Determine output format for BISL results
        print('\n' + "BISL Summary:")
        print(bisl.summarize(result))

