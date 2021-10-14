#!/usr/bin/env python
import sys

import os
import imp 
import click

import configparser as cnfg
import numpy as np


from . import cli_fk
from . import cli_fd
from . import cli_assoc
from . import cli_loc

@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    '''
    infrapy - Python-based Infrasound Data Analysis Toolkit

    More detailed description...
    '''
    pass

main.add_command(cli_fk.run_fk)
main.add_command(cli_fd.run_fd)
main.add_command(cli_assoc.run_assoc)
main.add_command(cli_loc.run_loc)

if __name__ == '__main__':
    main()
