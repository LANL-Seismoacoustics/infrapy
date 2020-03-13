#!/usr/bin/env python
import sys

import os
import click

from . import run_fk as fk

from . import run_fd as fdet
from . import run_assoc as assoc
from . import run_loc as loc


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main(args=None):
    '''
    infrapy - Python-based Infrasound Data Analysis Toolkit

    More detailed description...
    '''
    pass

@main.command('run_fk', short_help="Run beamforming methods on waveform data")
@click.option("--config_file", help="configuration file (required)", prompt = "Specify configuration file")
@click.option("--array", help="override to array in configuration file", default=None)
def run_fk(config_file, array):
    '''
    Run the fk beamforming methods

    More detailed description...
    '''

    if array:
        click.echo("Running fk beamforming with configuration file: '" + config_file + "' with array override: '" + array + "'...")
    else:
        click.echo("Running fk beamforming with configuration file: '" + config_file + "'...")

    if os.path.isfile(config_file):
        fk.run(config_file, array)
    else:
        click.echo("Invalid configuration file")



@main.command('run_fd', short_help="Identify detections from beamforming results")
@click.option("--config_file", help="configuration file (required)", prompt = "Specify configuration file")
@click.option("--array", help="override to array in configuration file", default=None)
def run_fd(config_file, array):
    '''
    Identify detections

    More detailed description...
    '''

    if array:
        click.echo("Running f-detector with configuration file: '" + config_file + "' with array override: '" + array + "'...")
    else:
        click.echo("Running f-detector with configuration file: '" + config_file + "'...")

    if os.path.isfile(config_file):
        fdet.run(config_file, array)
    else:
        click.echo("Invalid configuration file")


@main.command('run_assoc', short_help="Associate detections into events")
@click.option("--config_file", help="configuration file (required)", prompt = "Specify configuration file")
def run_assoc(config_file):
    '''
    Associate detections

    More detailed description...
    '''

    click.echo("Running association methods with configuration file: '" + config_file + "'...")

    if os.path.isfile(config_file):
        assoc.run(config_file)
    else:
        click.echo("Invalid configuration file")

@main.command('run_loc', short_help="Estimate source locations and times for events")
@click.option("--config_file", help="configuration file (required)", prompt = "Specify configuration file")
def run_loc(config_file):
    '''
    Estiamte source locations and times for events


    More detailed description...
    '''

    click.echo("Running location methods with configuration file: '" + config_file + "'...")

    if os.path.isfile(config_file):
        loc.run(config_file)
    else:
        click.echo("Invalid configuration file")



if __name__ == '__main__':
    main()
