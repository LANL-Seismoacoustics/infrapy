#!/usr/bin/env python

import os
import click
import fnmatch
import tempfile

import numpy as np 
import configparser as cnfg

from obspy import UTCDateTime

from ..utils import config

# Function to find common sub-string in a list of strings
def comm_str(str_list):
    if len(str_list) > 1:                    
        str0 = str_list[0]
        result = "" 

        for i in range(len(str0)):
            for j in range(i + 1, len(str0) + 1):
                k = 1
                for k in range(1, len(str_list)):
                    if str0[i:j] not in str_list[k]:
                        break
                if (k + 1 == len(str_list) and len(result) < (j - i)):
                    result = str0[i:j]
    else:
        result = str_list[0]

    return result


@click.command('infrapype', short_help="Automated infrapy analysis pipeline (prototype)",context_settings={'help_option_names': ['-h', '--help']})
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--out-label", help="Specify a file output prefix (default 'YYYY-MM-DD')", default=None)
@click.option("--starttime", help="Start time of automated analysis", default=None)
@click.option("--endtime", help="End time of automated analysis", default=None)


@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def pipeline(config_file, out_label, starttime, endtime, cpu_cnt):
    '''
    Run infrapy pipeline (infrapype) analysis 

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapype --config-file config/pipeline-HRR5.cnfg
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##            InfraPype            ##")
    click.echo("##   Automated Analysis Pipeline   ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    test_commands = False

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo('\n' + "Invalid configuration file (file not found)")
            return 0
    else:
        click.echo('\n' + "Pipeline methods require configuration file.")
        return 0
    
    pipe_params = dict(user_config["PIPELINE"])

    # Read in any command line options
    pipe_params['out_label'] = config.set_param(user_config, 'PIPELINE', 'out_label', out_label, 'string')
    pipe_params['starttime'] = config.set_param(user_config, 'PIPELINE', 'starttime', starttime, 'string')
    pipe_params['endtime'] = config.set_param(user_config, 'PIPELINE', 'endtime', endtime, 'string')

    if pipe_params['out_label'] is None:
        t0 = UTCDateTime(pipe_params['starttime'])
        pipe_params['out_label'] = str(t0.datetime).replace(" ","T").replace(":",".")

    # summarize pipeline parameters
    click.echo('\n' + "pipeline parameters:")
    for key in pipe_params.keys():       
        if 'trace_ids' in key:
            if '|' in pipe_params[key]:
                pipe_params[key] = pipe_params[key].replace('\n','').split('|')
            else:
                pipe_params[key] = pipe_params[key].replace('\n','').split(',')

        if "," in pipe_params[key]:
            pipe_params[key] = pipe_params[key].replace('\n','').split(',')

        if 'cnfgs' in key and not isinstance(pipe_params[key], list):
            pipe_params[key] = [pipe_params[key]]
        
        if type(pipe_params[key]) is list:
            if len(pipe_params[key]) > 1:
                click.echo("  " + key + ":")
                for val in pipe_params[key]:
                    click.echo("    " + str(val))
            else:
                click.echo("  " + key + ": " + str(pipe_params[key][0]))
        else:
            click.echo("  " + key + ": " + str(pipe_params[key]))
    click.echo("")

    # Build directories (if needed)
    for dir_key in ['det_dir', 'ev_dir', 'figs_dir']:
        dir = pipe_params[dir_key]

        if not os.path.isdir(dir):
            click.echo("Creating directory: " + dir)
            os.makedirs(dir)

    # Define temporary directory (if needed)
    with tempfile.TemporaryDirectory(prefix='infraga_') as temp_path:

        if 'temp_dir' in pipe_params.keys():
            if not os.path.isdir(pipe_params['temp_dir']):
                os.makedirs(pipe_params['temp_dir'])               
            temp_path = pipe_params['temp_dir']

        if temp_path[-1] != "/":
            temp_path = temp_path + "/"
        
        # Run detection, merge, and plot
        for id in pipe_params['trace_ids']:
            if '*' in id:
                net, sta, loc, cha = id.split(".")               
                det_label = pipe_params['out_label'] + "_" + id.replace("*","")
            else:
                net = np.unique([item.split(".")[0] for item in id.replace(" ","").split(",")])
                sta = np.unique([item.split(".")[1] for item in id.replace(" ","").split(",")])
                loc = np.unique([item.split(".")[2] for item in id.replace(" ","").split(",")])
                cha = np.unique([item.split(".")[3] for item in id.replace(" ","").split(",")])

                det_label = pipe_params['out_label'] + "_" + comm_str(net) + "." + comm_str(sta) + "." + comm_str(loc) + "." + comm_str(cha)

                net = ",".join(net)
                sta = ",".join(sta)
                loc = ",".join(loc)
                cha = ",".join(cha)

            if not os.path.isfile(pipe_params["det_dir"] + det_label + ".dets.json.gz"):
                for bm_j, bm_config in enumerate(pipe_params["beam_cnfgs"]):

                    command = "infrapy detect beam"
                    command = command + " --fdsn " + pipe_params["fdsn"]
                    command = command + " --network '" + net + "' --station '" + sta + "'"
                    command = command + " --location '" + loc + "' --channel '" + cha + "'"
                    command = command + " --starttime " + pipe_params["starttime"]
                    command = command + " --endtime " + pipe_params["endtime"]

                    command = command + " --config-file " + pipe_params["config_dir"] + bm_config
                    command = command + " --detect-label " + temp_path + det_label + "-" + str(bm_j)
                    click.echo(command)
                    if not test_commands:
                        os.system(command)

                command = "infrapy utils merge-dets --det-files '" + temp_path + det_label + "-*' --merged-label " + pipe_params["det_dir"] + det_label
                click.echo(command)
                if not test_commands:
                    os.system(command)

            fig_label = det_label.replace(".","_")              
            if not os.path.isfile(pipe_params["figs_dir"] + fig_label + "_det0.png"):
                command = "infrapy plot beam --det-file " + pipe_params["det_dir"] + det_label + ".dets.json.gz"
                command = command + " --figure-out " + pipe_params["figs_dir"] + fig_label
                command = command + " --plot-all-dets true  --show-figure False"
                click.echo(command) 
                if not test_commands:
                    os.system(command)
                click.echo("")


        # Build and analyze events if config file(s) are specified 
        if "ev_build_cnfgs" not in pipe_params.keys(): return 

        for ev_j, ev_config in enumerate(pipe_params["ev_build_cnfgs"]):
            ev_j_id = "_ev-bld-cnfg-" + str(ev_j) if len(pipe_params['ev_build_cnfgs']) > 1 else "" 

            if not os.path.isfile(pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id + ".ev.json.gz"):
                command = "infrapy event build --detect-files '" + pipe_params["det_dir"] + "*.dets.json.gz' --event-label " + pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id
                command = command + " --config-file " + pipe_params["config_dir"] + ev_config

                print(command) 
                if not test_commands:
                    os.system(command)

        
        # Cycle through events, plot the projections and compute localizations
        ev_files = [file for file in np.sort(os.listdir(pipe_params["ev_dir"])) if fnmatch.fnmatch(file,"*.ev.json.gz")]

        for ev_file in ev_files:
            command = "infrapy plot map_dets --event-file " + pipe_params["ev_dir"] + ev_file + " --figure-out " + pipe_params["figs_dir"] + ev_file.split(".")[0] + ".back-proj.png --show-figure False"
            print('\n' + command) 
            if not test_commands:
                os.system(command)

            for k, loc_config in enumerate(pipe_params["ev_loc_cnfgs"]):
                command = "infrapy event localize --event-file " + pipe_params["ev_dir"] + ev_file + " --config-file " + pipe_params["config_dir"] + loc_config
                print('\n' + command) 
                if not test_commands:
                    os.system(command)

                command = "infrapy plot localize --event-file " + pipe_params["ev_dir"] + ev_file + "  --event-loc-index " + str(k) 
                command = command + " --figure-out " + pipe_params["figs_dir"] + ev_file.split(".")[0] + ".loc-" + str(k) + ".png --show-figure False"
                print(command) 
                if not test_commands:
                    os.system(command)

            for l, char_config in enumerate(pipe_params["ev_char_cnfgs"]):
                # run yield estimation
                command = "infrapy event characterize --event-file HRR/events/HRR5_0-0.ev.json.gz --config-file " + pipe_params["config_dir"] + char_config
                print(command)
                if not test_commands:
                    os.system(command)


if __name__ == '__main__':
    pipeline()
