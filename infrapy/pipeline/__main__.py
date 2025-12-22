#!/usr/bin/env python

import os
import click
import fnmatch

import numpy as np 
import configparser as cnfg


@click.command('infrapype', short_help="Automated infrapy analysis pipeline (prototype)",context_settings={'help_option_names': ['-h', '--help']})
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--out-label", help="Specify a file output prefix (default YR-JDAY)", default=None)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def pipeline(config_file, out_label, cpu_cnt):
    '''
    Run infrapy pipeline (infrapype) analysis 

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapype --config-file config/infrapype_HRR5.cnfg
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##            InfraPype            ##")
    click.echo("##   Automated Analysis Pipeline   ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

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

    click.echo('\n' + "pipeline parameters:")
    for key in pipe_params.keys():
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


    ###################
    ## RUN DETECTION ##
    ###################
    test_commands = True
    # ev_label = "test"

    if ev_label is None:
        ev_label = "test2"


    # Run detection, merge, and plot
    for bm_j, bm_det_config in enumerate(pipe_params["beam_cnfgs"]):
        for id in pipe_params['trace_ids']:

            # Need to add some parsing logic here for lists or wildcards
            net, sta, loc, cha = id.split(".")
            output_label = id.replace("*","")

            if not os.path.isfile(pipe_params["temp_dir"] +  output_label  + "-" + str(bm_j) + ".dets.json.gz"):
                command = "infrapy detect beam"
                command = command + " --fdsn " + pipe_params["fdsn"]
                command = command + " --network '" + net + "' --station '" + sta + "'"
                command = command + " --location '" + loc + "' --channel '" + cha + "'"
                command = command + " --starttime " + pipe_params["starttime"]
                command = command + " --endtime " + pipe_params["endtime"]

                command = command + " --config-file " + pipe_params["config_dir"] + bm_det_config
                command = command + " --detect-label " + pipe_params["temp_dir"] + output_label + "-" + str(bm_j)
                print(command)
                if not test_commands:
                    os.system(command)
        

    # Merge detections together and plot
    for id in pipe_params['trace_ids']:
        net, sta, loc, cha = id.split(".")               
        output_label = id.replace("*","")
        fig_label = output_label.replace(".","_")              

        if not os.path.isfile(pipe_params["det_dir"] + output_label + ".dets.json.gz"):
            command = "infrapy utils merge-dets --det-files '"+ pipe_params["temp_dir"] + output_label + "-*' --merged-label " + pipe_params["det_dir"] + output_label
            print(command)
            if not test_commands:
               os.system(command)

        if not os.path.isfile(pipe_params["figs_dir"] + fig_label + "_det0.png"):
            command = "infrapy plot beam --det-file " + pipe_params["det_dir"] + output_label + ".dets.json.gz"
            command = command + " --figure-out " + pipe_params["figs_dir"] + fig_label
            command = command + " --plot-all-dets true  --show-figure False"
            print(command) 
            if not test_commands:
                os.system(command)

    ##################
    ## BUILD EVENTS ##
    ##################

    for ev_j, ev_config in enumerate(pipe_params["ev_build_cnfgs"]):
        if not os.path.isfile(pipe_params["ev_dir"] + ev_label + "-" + str(ev_j) + ".ev.json.gz"):
            command = "infrapy event build --detect-files '" + pipe_params["det_dir"] + "*.dets.json.gz' --event-label " + pipe_params["ev_dir"] + ev_label + "_" + str(ev_j)
            command = command + " --config-file " + pipe_params["config_dir"] + ev_config

            print(command) 
            if not test_commands:
                os.system(command)


    # cycle through events, plot the projections and compute localizations
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
