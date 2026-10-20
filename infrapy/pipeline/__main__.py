#!/usr/bin/env python

import os
import click
import fnmatch
import tempfile
import json
import gzip

import numpy as np
import configparser as cnfg

from obspy import UTCDateTime

from ..utils import config, data_io
from ..detection import beam

from ..detection import visualization as det_vis
from ..location import visualization as loc_vis


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


@click.command('infrapype2', short_help="Automated infrapy analysis pipeline (prototype)",context_settings={'help_option_names': ['-h', '--help']})
@click.option("--cnfg-file", help="Configuration file", default=None)
@click.option("--out-label", help="Specify a file output prefix (default 'YYYY-MM-DD')", default=None)
@click.option("--starttime", help="Start time of automated analysis", default=None)
@click.option("--endtime", help="End time of automated analysis", default=None)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def pipeline_old(cnfg_file, out_label, starttime, endtime, cpu_cnt):
    '''
    Run infrapy pipeline (infrapype) analysis

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapype --cnfg-file config/pipeline-HRR5.cnfg
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

                    command = command + " --cnfg-file " + pipe_params["config_dir"] + bm_config
                    command = command + " --det-label " + temp_path + det_label + "-" + str(bm_j)
                    click.echo(command)
                    if not test_commands:
                        os.system(command)

                command = "infrapy utils merge_dets --dets-files '" + temp_path + det_label + "-*' --merged-label " + pipe_params["det_dir"] + det_label
                click.echo(command)
                if not test_commands:
                    os.system(command)

            fig_label = det_label.replace(".","_")
            if not os.path.isfile(pipe_params["figs_dir"] + fig_label + "_det0.png"):
                command = "infrapy plot detect --dets-file " + pipe_params["det_dir"] + det_label + ".dets.json.gz"
                command = command + " --figure-out " + pipe_params["figs_dir"] + fig_label
                command = command + " --plot-all true  --show-figure False"
                click.echo('\n' + command)
                if not test_commands:
                    os.system(command)
                click.echo("")

        # Build and analyze events if config file(s) are specified
        if "ev_build_cnfgs" not in pipe_params.keys(): return

        for ev_j, ev_config in enumerate(pipe_params["ev_build_cnfgs"]):
            ev_j_id = "_ev-bld-cnfg-" + str(ev_j) if len(pipe_params['ev_build_cnfgs']) > 1 else ""

            if not os.path.isfile(pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id + ".ev.json.gz"):
                command = "infrapy event build --dets-files '" + pipe_params["det_dir"] + pipe_params['out_label'] + "*.dets.json.gz' --ev-label " + pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id
                command = command + " --cnfg-file " + pipe_params["config_dir"] + ev_config

                click.echo('\n' + command)
                if not test_commands:
                    os.system(command)

        # Cycle through events, plot the projections and compute localizations
        if "ev_loc_cnfgs" not in pipe_params.keys(): return
        ev_files = [file for file in np.sort(os.listdir(pipe_params["ev_dir"])) if fnmatch.fnmatch(file,"*.ev.json.gz")]
        for ev_file in ev_files:
            fig_path = pipe_params["figs_dir"] + ev_file.split(".ev")[0]

            command = "infrapy plot event --ev-file " + pipe_params["ev_dir"] + ev_file + " --figure-out " + fig_path + ".back-proj.png --show-figure False"
            click.echo('\n' + command)
            if not test_commands:
                os.system(command)

            for k, loc_config in enumerate(pipe_params["ev_loc_cnfgs"]):
                command = "infrapy event locate --ev-file " + pipe_params["ev_dir"] + ev_file + " --cnfg-file " + pipe_params["config_dir"] + loc_config
                click.echo('\n' + command)
                if not test_commands:
                    os.system(command)

                command = "infrapy plot event --ev-file " + pipe_params["ev_dir"] + ev_file + "  --loc-index " + str(k)
                command = command + " --figure-out " + fig_path + ".loc-" + str(k) + ".png --show-figure False"
                click.echo(command)
                if not test_commands:
                    os.system(command)

            if "ev_char_cnfgs" not in pipe_params.keys(): return
            for l, char_config in enumerate(pipe_params["ev_char_cnfgs"]):
                # run yield estimation
                command = "infrapy event characterize --ev-file " + pipe_params["ev_dir"] + ev_file + " --cnfg-file " + pipe_params["config_dir"] + char_config
                click.echo(command)
                if not test_commands:
                    os.system(command)


def set_bm_config(cnfg_file):

    bm_config = cnfg.ConfigParser()
    bm_config.read(cnfg_file)

    bm_params = {}
    bm_params['freq_min'] = config.set_param(bm_config, 'FK', 'freq_min', None, 'float')
    bm_params['freq_max'] = config.set_param(bm_config, 'FK', 'freq_max', None, 'float')
    bm_params['back_az_min'] = config.set_param(bm_config, 'FK', 'back_az_min', None, 'float')
    bm_params['back_az_max'] = config.set_param(bm_config, 'FK', 'back_az_max', None, 'float')
    bm_params['back_az_step'] = config.set_param(bm_config, 'FK', 'back_az_step', None, 'float')
    bm_params['trace_vel_min'] = config.set_param(bm_config, 'FK', 'trace_vel_min', None, 'float')
    bm_params['trace_vel_max'] = config.set_param(bm_config, 'FK', 'trace_vel_max', None, 'float')
    bm_params['trace_vel_step'] = config.set_param(bm_config, 'FK', 'trace_vel_step', None, 'float')
    bm_params['method'] = config.set_param(bm_config, 'FK', 'method', None, 'string')
    bm_params['signal_start'] = config.set_param(bm_config, 'FK', 'signal_start', None, 'string')
    bm_params['signal_end'] = config.set_param(bm_config, 'FK', 'signal_end', None, 'string')
    bm_params['noise_start'] = config.set_param(bm_config, 'FK', 'noise_start', None, 'string')
    bm_params['noise_end'] = config.set_param(bm_config, 'FK', 'noise_end', None, 'string')
    bm_params['window_len'] = config.set_param(bm_config, 'FK', 'window_len', None, 'float')
    bm_params['sub_window_len'] = config.set_param(bm_config, 'FK', 'sub_window_len', None, 'float')
    bm_params['window_step'] = config.set_param(bm_config, 'FK', 'window_step', None, 'float')
    bm_params['cpu_cnt'] = config.set_param(bm_config, 'FK', 'cpu_cnt', None, 'int')

    det_params = {}
    det_params['window_len'] = config.set_param(bm_config, 'FD', 'window_len', None, 'float')
    det_params['p_value'] = config.set_param(bm_config, 'FD', 'p_value', None, 'float')
    det_params['min_duration'] = config.set_param(bm_config, 'FD', 'min_duration', None, 'float')
    det_params['back_az_width'] = config.set_param(bm_config, 'FD', 'back_az_width', None, 'float')
    det_params['fixed_thresh'] = config.set_param(bm_config, 'FD', 'fixed_thresh', None, 'float')
    det_params['thresh_ceil'] = config.set_param(bm_config, 'FD', 'thresh_ceil', None, 'float')
    det_params['merge_dets'] = config.set_param(bm_config, 'FD', 'merge_dets', None, 'bool')

    return bm_params, det_params


@click.command('infrapype', short_help="Automated infrapy analysis pipeline (prototype)",context_settings={'help_option_names': ['-h', '--help']})
@click.option("--cnfg-file", help="Configuration file", default=None)
@click.option("--out-label", help="Specify a file output prefix (default 'YYYY-MM-DDTHH.MM.SS')", default=None)
@click.option("--starttime", help="Start time of automated analysis", default=None)
@click.option("--endtime", help="End time of automated analysis", default=None)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: 1)", default=1, type=int)
def pipeline(cnfg_file, out_label, starttime, endtime, cpu_cnt):
    '''
    Run infrapy pipeline (infrapype) analysis

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapype2 --cnfg-file config/pipeline-HRR5.cnfg
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

    # Read in any command line options
    pipe_params['out_label'] = config.set_param(user_config, 'PIPELINE', 'out_label', out_label, 'string')
    pipe_params['starttime'] = config.set_param(user_config, 'PIPELINE', 'starttime', starttime, 'string')
    pipe_params['endtime'] = config.set_param(user_config, 'PIPELINE', 'endtime', endtime, 'string')
    pipe_params['cpu_cnt'] = config.set_param(user_config, 'PIPELINE', 'cpu_cnt', cpu_cnt, 'int')

    # Check data source (local files, FDSN, DB)
    pipe_params['local_wvfrm_data'] = config.set_param(user_config, 'PIPELINE', 'local_wvfrm_data', None, 'string')
    pipe_params['fdsn'] = config.set_param(user_config, 'PIPELINE', 'fdsn', None, 'string')
    pipe_params['db_config'] = config.set_param(user_config, 'PIPELINE', 'db_config', None, 'string')

    # summarize pipeline parameters
    click.echo('\n' + "pipeline parameters:")
    for key in pipe_params.keys():
        if pipe_params[key] is not None:
            if type(pipe_params[key]) is str:
                if '|' in pipe_params[key]:
                    pipe_params[key] = pipe_params[key].replace('\n','').replace(' ','').split('|')

                if "," in pipe_params[key]:
                    pipe_params[key] = pipe_params[key].replace('\n','').replace(' ','').split(',')

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

    # Set automatic output label if needed
    if pipe_params['starttime'] is not None:
        t0 = UTCDateTime(pipe_params['starttime'])
        pipe_params['out_label'] = str(t0.datetime).replace(" ","T").replace(":",".")
    else:
        click.echo('\n' + "Warning!! Configuration file does not incluce 'starttime' or 'out_label' parameter.")
        click.echo('\t' + "Using general 'infrapype' output label.")
        pipe_params['out_label'] = 'infrapype'

    # Build directories if needed
    for dir_key in ['det_dir', 'ev_dir', 'figs_dir']:
        if not os.path.isdir(pipe_params[dir_key]):
            click.echo("Creating directory: " + pipe_params[dir_key])
            os.makedirs(pipe_params[dir_key])

    # Define temporary directory (if needed)
    with tempfile.TemporaryDirectory(prefix='infraga_') as temp_path:
        if 'temp_dir' in pipe_params.keys():
            if not os.path.isdir(pipe_params['temp_dir']):
                os.makedirs(pipe_params['temp_dir'])
            temp_path = pipe_params['temp_dir']
        if temp_path[-1] != "/":
            temp_path = temp_path + "/"

        # Extract waveform info and set up for iterating over
        if pipe_params['local_wvfrm_data'] is not None:
            wvfrm_sets = [{'local_wvfrms': wvfrms, "fdsn": None, "db_config" : None,
                           'net' : None, 'sta' : None, 'loc' : None, 'cha' : None} for wvfrms in pipe_params['local_wvfrm_data']]

        elif pipe_params['fdsn'] is not None or pipe_params['db_config'] is not None:
            # Parse network, station, channel, location info and set local_wvfrm to Nones
            wvfrm_sets = []
            for id in pipe_params['trace_ids']:
                wvfrm_sets = wvfrm_sets + [{"local_wvfrms" : None, "fdsn": pipe_params['fdsn'], "db_config" : pipe_params['db_config']}]
                if '*' in id:
                    wvfrm_sets[-1]['net'] = id.split(".")[0]
                    wvfrm_sets[-1]['sta'] = id.split(".")[1]
                    wvfrm_sets[-1]['loc'] = id.split(".")[2]
                    wvfrm_sets[-1]['cha'] = id.split(".")[3]
                else:
                    wvfrm_sets[-1]['net'] = ",".join(np.unique([item.split(".")[0] for item in id.replace(" ","").split(",")]))
                    wvfrm_sets[-1]['sta'] = ",".join(np.unique([item.split(".")[1] for item in id.replace(" ","").split(",")]))
                    wvfrm_sets[-1]['loc'] = ",".join(np.unique([item.split(".")[2] for item in id.replace(" ","").split(",")]))
                    wvfrm_sets[-1]['cha'] = ",".join(np.unique([item.split(".")[3] for item in id.replace(" ","").split(",")]))

        # Run detection, merge, and plot
        for info in wvfrm_sets:
            click.echo('\n' + '=' * 50 + '\n' + '=' * 50)
            stream, latlon = data_io.set_stream(info['local_wvfrms'], info["fdsn"], info["db_config"],
                                                info['net'], info['sta'], info['loc'], info['cha'],
                                                pipe_params['starttime'], pipe_params['endtime'], None)

            wvfrm_info = data_io.wvfrm_info(stream, latlon)
            click.echo("Data summary:")
            for tr in stream:
                click.echo('  ' + tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))
            click.echo('')

            det_label = pipe_params['out_label'] + "_"
            det_label = det_label + comm_str(np.unique([tr.id.split(".")[0] for tr in stream])) + "."
            det_label = det_label + comm_str(np.unique([tr.id.split(".")[1] for tr in stream])) + "."
            det_label = det_label + comm_str(np.unique([tr.id.split(".")[2] for tr in stream])) + "."
            det_label = det_label + comm_str(np.unique([tr.id.split(".")[3] for tr in stream]))

            # Run beamforming and combine results
            if not os.path.isfile(pipe_params["det_dir"] + det_label + ".dets.json.gz"):
                for bm_j, bm_config_file in enumerate(pipe_params["beam_cnfgs"]):
                    click.echo("Running beam analysis using " + bm_config_file + '\t', nl=False)

                    bm_params, det_params = set_bm_config(pipe_params["config_dir"] + bm_config_file)
                    beam_times, beam_peaks = beam.run_fk_dict(stream, latlon, bm_params, None, cpu_cnt, verbose=False)
                    dets, _ = beam.run_afd_dict(beam_times, beam_peaks, bm_params, det_params, len(stream))
                    dets_out = [beam.det2dict(stream, latlon, beam_times, beam_peaks, bm_params, det_info) for det_info in dets]

                    det_output = {'wvfrm_info' : [wvfrm_info],
                                  'fk_params' : [bm_params],
                                  'det_params' : [det_params],
                                  'det_info' : dets_out}

                    with gzip.open(temp_path + det_label + "-" + str(bm_j) + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
                        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

                click.echo('')
                dets_all = data_io._load_dets_json(temp_path + det_label + "-*.dets.json.gz")
                det_combined = []
                for entry in dets_all:
                    for det in entry["det_info"]:
                        det_combined = det_combined + [det]
                        det_combined[-1]["wvfrm_info"] = entry["wvfrm_info"]
                        det_combined[-1]["fk_params"] = entry["fk_params"]
                        det_combined[-1]["det_params"] = entry["det_params"]

                det_output = {'det_info' : beam.merge_dets(det_combined)}
                with gzip.open(pipe_params["det_dir"] + det_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
                    json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

            if not os.path.isfile(pipe_params["figs_dir"] + det_label.replace(".","_") + "_det0.png"):
                det_data = json.load(gzip.open(pipe_params["det_dir"] + det_label + ".dets.json.gz", 'rt'))
                for k in range(len(det_data['det_info'])):
                    det = det_data["det_info"][k]
                    det_vis.plot_det_json(det_data, k, output_path=pipe_params["figs_dir"] + det_label.replace(".","_") + "_det" + str(k), show_fig=False)

        '''
        # Build and analyze events if config file(s) are specified
        if "ev_build_cnfgs" not in pipe_params.keys(): return

        for ev_j, ev_config in enumerate(pipe_params["ev_build_cnfgs"]):
            ev_j_id = "_ev-bld-cnfg-" + str(ev_j) if len(pipe_params['ev_build_cnfgs']) > 1 else ""

            if not os.path.isfile(pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id + ".ev.json.gz"):
                command = "infrapy event build --dets-files '" + pipe_params["det_dir"] + pipe_params['out_label'] + "*.dets.json.gz' --ev-label " + pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id
                command = command + " --cnfg-file " + pipe_params["config_dir"] + ev_config

                click.echo('\n' + command)
                if not test_commands:
                    os.system(command)

        # Cycle through events, plot the projections and compute localizations
        if "ev_loc_cnfgs" not in pipe_params.keys(): return
        ev_files = [file for file in np.sort(os.listdir(pipe_params["ev_dir"])) if fnmatch.fnmatch(file,"*.ev.json.gz")]
        for ev_file in ev_files:
            fig_path = pipe_params["figs_dir"] + ev_file.split(".ev")[0]

            command = "infrapy plot event --ev-file " + pipe_params["ev_dir"] + ev_file + " --figure-out " + fig_path + ".back-proj.png --show-figure False"
            click.echo('\n' + command)
            if not test_commands:
                os.system(command)

            for k, loc_config in enumerate(pipe_params["ev_loc_cnfgs"]):
                command = "infrapy event locate --ev-file " + pipe_params["ev_dir"] + ev_file + " --cnfg-file " + pipe_params["config_dir"] + loc_config
                click.echo('\n' + command)
                if not test_commands:
                    os.system(command)

                command = "infrapy plot event --ev-file " + pipe_params["ev_dir"] + ev_file + "  --loc-index " + str(k)
                command = command + " --figure-out " + fig_path + ".loc-" + str(k) + ".png --show-figure False"
                click.echo(command)
                if not test_commands:
                    os.system(command)

            if "ev_char_cnfgs" not in pipe_params.keys(): return
            for l, char_config in enumerate(pipe_params["ev_char_cnfgs"]):
                # run yield estimation
                command = "infrapy event characterize --ev-file " + pipe_params["ev_dir"] + ev_file + " --cnfg-file " + pipe_params["config_dir"] + char_config
                click.echo(command)
                if not test_commands:
                    os.system(command)
        '''

if __name__ == '__main__':
    pipeline2()
