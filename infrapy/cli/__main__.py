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


@click.group('run_spye', short_help="Estimate explosive yield from a surface explosion", context_settings={'help_option_names': ['-h', '--help']})
def run_spye():
    '''
    infrapy run_spye - explosive yield estimation methods
    
    '''
    pass 


@click.group('plot', short_help="Visualize infrapy analysis results", context_settings={'help_option_names': ['-h', '--help']})
def plot():
    '''
    infrapy plot - visualization methods for analysis results
    
    '''
    pass 


@click.group('utils', short_help="Various utility functions for infrapy analysis", context_settings={'help_option_names': ['-h', '--help']})
def utils():
    '''
    infrapy utils - various utility functions for infrapy usage
    
    '''
    pass 


main.add_command(run_spye)
main.add_command(plot)
main.add_command(utils)

# main functions (run analysis)
main.add_command(cli_detection.run_fk)
main.add_command(cli_detection.run_fd)
main.add_command(cli_detection.run_fkd)
main.add_command(cli_detection.run_sd)
main.add_command(cli_assoc.run_assoc)
main.add_command(cli_loc.run_loc)

# SpYE
run_spye.add_command(cli_loc.regional)
run_spye.add_command(cli_loc.single_station)
run_spye.add_command(cli_loc.combine)

# visualizations
plot.add_command(cli_visualization.fk)
plot.add_command(cli_visualization.fd)
plot.add_command(cli_visualization.sd)

plot.add_command(cli_visualization.dets)
plot.add_command(cli_visualization.loc)
plot.add_command(cli_visualization.origin_time)
plot.add_command(cli_visualization.yield_plot)

# Utilities
utils.add_command(cli_utils.arrivals2json)
utils.add_command(cli_utils.arrival_time)
utils.add_command(cli_utils.calc_celerity)
utils.add_command(cli_utils.check_db_wvfrm)
utils.add_command(cli_utils.write_wvfrms)
utils.add_command(cli_utils.best_beam)


if __name__ == '__main__':
    main()
