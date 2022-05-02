#!/usr/bin/env python

import click

from . import cli_detection
from . import cli_assoc
from . import cli_loc
from . import cli_visualization
from . import cli_utils

@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    '''
    infrapy - Python-based Infrasound Signal Analysis Toolkit

    Command line interface (CLI) for running and visualizing infrasound analysis
    '''
    pass


@click.group('plot', short_help="Visualize infrapy analysis results", context_settings={'help_option_names': ['-h', '--help']})
def plot():
    '''
    infrapy plot - Visualization methods for analysis results
    
    '''
    pass 


@click.group('utils', short_help="Various utility functions for infrapy analysis", context_settings={'help_option_names': ['-h', '--help']})
def utils():
    '''
    infrapy utils - various utility functions for infrapy analysis
    
    '''
    pass 


main.add_command(plot)
main.add_command(utils)

# main functions (run analysis)
main.add_command(cli_detection.run_fk)
main.add_command(cli_detection.run_fd)
main.add_command(cli_detection.run_fkd)
main.add_command(cli_assoc.run_assoc)
main.add_command(cli_loc.run_loc)

# visualizations
plot.add_command(cli_visualization.fk)
plot.add_command(cli_visualization.fd)
plot.add_command(cli_visualization.dets)
plot.add_command(cli_visualization.loc)
plot.add_command(cli_visualization.origin_time)

# Utilities
utils.add_command(cli_utils.arrivals2json)
utils.add_command(cli_utils.arrival_time)
utils.add_command(cli_utils.calc_celerity)


if __name__ == '__main__':
    main()
