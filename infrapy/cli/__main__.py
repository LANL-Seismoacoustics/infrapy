#!/usr/bin/env python

import os 
import click
import webbrowser
import subprocess
import shlex

from importlib.util import find_spec 

from . import cli_detection
from . import cli_event
from . import cli_visualization
from . import cli_utils

from . import cli_assoc
from . import cli_loc


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    '''
    infrapy - Python-based Infrasound Signal Analysis Toolkit

    Command line interface (CLI) for running and visualizing infrasound analysis
    '''
    pass


@click.group('detect', short_help="Detect signatures in infrasound data", context_settings={'help_option_names': ['-h', '--help']})
def detect():
    '''
    infrapy detect - run detection analysis
    
    '''
    pass 


@click.group('event', short_help="Build and analyse events", context_settings={'help_option_names': ['-h', '--help']})
def event():
    '''
    infrapy detect - run detection analysis
    
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


#######################
##    Open Manual    ##
#######################
@click.command('doc', short_help="Open infrapy manual")
def open_doc():

    pkg_loc = find_spec('infrapy').submodule_search_locations[0]
    filename = pkg_loc + '/docs/build/html/index.html'

    if not os.path.isfile(filename):      
        print("Compiling manual...")
        subprocess.run(shlex.split("make html -C " + pkg_loc + "/docs/"), shell=False)

    webbrowser.open('file://' + os.path.realpath(filename), new=2)


main.add_command(open_doc)
main.add_command(detect)
main.add_command(event)
main.add_command(plot)
main.add_command(utils)

# Analysis methods
detect.add_command(cli_detection.run_beam_detect)
detect.add_command(cli_detection.run_spec_detect)

event.add_command(cli_event.build)
event.add_command(cli_event.localize)
event.add_command(cli_event.characterize)

# Visualizations
# plot.add_command(cli_visualization.beam_detect)
# plot.add_command(cli_visualization.spec_detect)
plot.add_command(cli_visualization.detect_combined)

plot.add_command(cli_visualization.wvfrms)
plot.add_command(cli_visualization.map_dets)
plot.add_command(cli_visualization.ev_loc)
plot.add_command(cli_visualization.ev_char)

# Utilities
utils.add_command(cli_utils.check_db_wvfrm)
utils.add_command(cli_utils.write_wvfrms)

utils.add_command(cli_utils.db2dets_json)
utils.add_command(cli_utils.db2ev_json)

utils.add_command(cli_utils.merge_dets)
utils.add_command(cli_utils.convert_dets)
utils.add_command(cli_utils.ev_gt)
utils.add_command(cli_utils.ev_summary)
utils.add_command(cli_utils.ev_loc_reset)

# moving to stochprop
utils.add_command(cli_utils.fit_celerity)






####################
## DEPRECATED CLI ##
####################
main.add_command(cli_detection.run_fk)
main.add_command(cli_detection.run_fd)
main.add_command(cli_detection.run_fkd)
main.add_command(cli_detection.run_sd)
main.add_command(cli_assoc.run_assoc)
main.add_command(cli_loc.run_loc)

# SpYE
@click.group('run_spye', short_help="Spectral yield methods", context_settings={'help_option_names': ['-h', '--help']}, hidden=True)
def run_spye():
    '''
    infrapy run_spye - run spectral yield methods
    '''
    pass 


main.add_command(run_spye)
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
utils.add_command(cli_utils.best_beam)





if __name__ == '__main__':
    main()
