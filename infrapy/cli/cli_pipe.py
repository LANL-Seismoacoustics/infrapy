#!/usr/bin/env python
import os 
import click
import json
import gzip
import fnmatch
import tempfile
import warnings

import configparser as cnfg
import numpy as np

from multiprocessing import Pool
from importlib.util import find_spec

import matplotlib.pyplot as plt 

from ..utils import config
from ..utils import data_io

from ..association import hjl

from ..location import bisl, tribl
from ..propagation import infrasound

from ..characterization import spye


@click.command('pipeline', short_help="Automated infrapy analysis pipeline (prototype)")
@click.option("--cnfg-file", help="Configuration file", default=None)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def pipeline(cnfg_file, cpu_cnt):
    '''
    Run infrapy pipeline (infrapype) analysis 

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapype --cnfg-file config/infrapype_HRR5.cnfg
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##            InfraPype            ##")
    click.echo("##   Automated Analysis Pipeline   ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if cnfg_file:
        click.echo('\n' + "Loading configuration info from: " + cnfg_file)
        if os.path.isfile(cnfg_file):
            user_config = cnfg.ConfigParser()
            user_config.read(cnfg_file)
        else:
            click.echo('\n' + "Invalid configuration file (file not found)")
            return 0
    else:
        click.echo('\n' + "Pipeline methods require configuration file.")
        return 0
    
    pipe_params = dict(user_config["PIPELINE"])

    print('\n' + "pipeline parameters:")
    for key in pipe_params.keys():
        if "," in pipe_params[key]:
            pipe_params[key] = pipe_params[key].replace('\n','').split(',')
        
        if 'cnfgs' in key and not isinstance(pipe_params[key], list):
            pipe_params[key] = [pipe_params[key]]

        print("  " + key + ": " + str(pipe_params[key]))
    print("")




