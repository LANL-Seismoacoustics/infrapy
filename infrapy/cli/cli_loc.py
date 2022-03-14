#!/usr/bin/env python
import sys

import click
import configparser as cnfg
import numpy as np
from scipy.sparse import data

from ..location import bisl
from ..propagation import infrasound

from ..utils import config
from ..utils import data_io

@click.command('run_loc', short_help="Estimate source locations and times for events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-detect-label", help="Detection path and pattern", default=None)
@click.option("--local-loc-label", help="Localization results path", default=None)
@click.option("--back-az-width", help="Width of beam projection (default: " + config.defaults['LOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Maximum source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['LOC']['resolution'] + ")", default=None, type=int)
@click.option("--src-est", help="Estimated source location and radius of region to consider (default: None)", default=None)
@click.option("--pgm-file", help="Path geometry model (PGM) file (default: None)", default=None)
def run_loc(config_file, local_detect_label, local_loc_label, back_az_width, range_max, resolution, src_est, pgm_file):
    '''
    Run Bayesian Infrasonic Source Localization (BISL) methods to estimate the source location and origin time for an event

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_loc --local-detect-label data/detection_set2.json --local-loc-label data/location2
    \tinfrapy run_loc --local-detect-label data/detection_set2.json --local-loc-label data/location2 --pgm-file ../infrapy/propagation/priors/UTTR_models/UTTR_06_1800UTC.pgm
    \tinfrapy run_loc --config-file config/loc_example.config
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
    local_detect_label = config.set_param(user_config, 'DETECTION IO', 'local_detect_label', local_detect_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  local_detect_label: " + str(local_detect_label))
    click.echo("  local_loc_label: " + str(local_loc_label))

    # Algorithm parameters
    back_az_width = config.set_param(user_config, 'LOC', 'back_az_width', back_az_width, 'float')
    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    resolution = config.set_param(user_config, 'LOC', 'resolution', resolution, 'int')
    src_est = config.set_param(user_config, 'LOC', 'src_est', src_est, 'string')
    pgm_file = config.set_param(user_config, 'LOC', 'pgm_file', pgm_file, 'str')

    if src_est is not None:
        src_est = [float(x.strip('[( )]')) for x in src_est.split(',')]

    click.echo('\n' + "Parameter summary:")
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  range_max: " + str(range_max))
    click.echo("  resolution: " + str(resolution))
    click.echo("  src_est: " + str(src_est))
    click.echo("  pgm_file: " + str(pgm_file))

    if pgm_file is not None:
        click.echo("")
        pgm = infrasound.PathGeometryModel()
        pgm.load(pgm_file)
    else:
        pgm = None

    click.echo("")
    events = data_io.set_det_list(local_detect_label, merge=False)
    if type(events[0]) is list:
        # run localization analysis for multiple detection sets
        for j, det_list in enumerate(events):
            click.echo('\n' + "Running BISL on event " + str(j + 1) + " of " + str(len(events)))
            result = bisl.run(det_list, path_geo_model=pgm, custom_region=src_est, resol=resolution, bm_width=back_az_width, rng_max=range_max, rad_min=100.0, rad_max=range_max/4.0)

            # Determine output format for BISL results
            click.echo('\n' + "BISL Summary:")
            click.echo(bisl.summarize(result))

            click.echo("Writing localization result into " + local_loc_label + "_ev-" + str(j + 1) + "loc.json")
            data_io.write_locs(result, local_loc_label + "_ev-" + str(j + 1) + "loc.json")

    else:
        # run a single localization analysis
        click.echo("")
        result = bisl.run(events, path_geo_model=pgm, custom_region=src_est, resol=resolution, bm_width=back_az_width, rng_max=range_max, rad_min=100.0, rad_max=range_max/4.0)

        # Determine output format for BISL results
        click.echo('\n' + "BISL Summary:")
        click.echo(bisl.summarize(result))

        click.echo("Writing localization result into " + local_loc_label + ".loc.json")
        data_io.write_locs(result, local_loc_label + ".loc.json")

