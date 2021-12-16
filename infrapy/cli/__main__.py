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
    infrapy - Python-based Infrasound Data Analysis Toolkit

    Command line interface (CLI) for running and visualizing infrasound analysis
    '''
    pass

main.add_command(cli_detection.run_fk)
main.add_command(cli_detection.run_fd)
main.add_command(cli_detection.run_fkd)
main.add_command(cli_assoc.run_assoc)
main.add_command(cli_loc.run_loc)

main.add_command(cli_visualization.plot_fk)
main.add_command(cli_visualization.plot_fd)
main.add_command(cli_visualization.plot_dets)
main.add_command(cli_visualization.plot_loc)
main.add_command(cli_visualization.plot_origin_time)

main.add_command(cli_utils.arrivals2json)


if __name__ == '__main__':
    main()
