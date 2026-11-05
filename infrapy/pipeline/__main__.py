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
from ..propagation import infrasound

from ..detection import beam
from ..association import hjl
from ..cli import cli_event
from ..characterization import spye

from importlib.util import find_spec
if find_spec('stochprop'):
    import stochprop.propagation as sp_prop


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

    # Read in command line options
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
    if pipe_params['out_label'] is None:
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

        #################
        ##  DETECTION  ##
        #################
        click.echo('\n' + '*' * 50 + '\n**' + ' ' * 14 + 'Detection Analysis' + ' ' * 14 + '**\n' + '*' * 50)
        for k, info in enumerate(wvfrm_sets):
            if k > 0:
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

                    bm_params, det_params = config.set_bm_config(pipe_params["config_dir"] + bm_config_file)
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
            else:
                click.echo("Existing results found for " + pipe_params["det_dir"] + det_label + ".dets.json.gz")

            '''
            # maybe move visualization out of the pipeline loop?
            # make a separate plotting pipeline method?
            if not os.path.isfile(pipe_params["figs_dir"] + det_label.replace(".","_") + "_det0.png"):
                det_data = json.load(gzip.open(pipe_params["det_dir"] + det_label + ".dets.json.gz", 'rt'))
                for k in range(len(det_data['det_info'])):
                    det = det_data["det_info"][k]
                    det_vis.plot_det_json(det_data, k, output_path=pipe_params["figs_dir"] + det_label.replace(".","_") + "_det" + str(k), show_fig=False)
            '''

        ####################
        ## EVENT BUILDING ##
        ####################
        if "ev_build_cnfgs" not in pipe_params.keys():
            return

        click.echo('\n' + '*' * 50 + '\n**' + ' ' * 16 + 'Event Analysis' + ' ' * 16 + '**\n' + '*' * 50 + '\n')
        for ev_j, ev_config in enumerate(pipe_params["ev_build_cnfgs"]):
            ev_j_id = "_ev-bld-cnfg-" + str(ev_j) if len(pipe_params['ev_build_cnfgs']) > 1 else ""

            if not os.path.isfile(pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id + "-0.ev.json.gz"):
                ev_label =  pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id
                ev_bld_params = config.set_ev_bld_config(pipe_params["config_dir"] + ev_config)
                ev_bld_params['starttime'], ev_bld_params['endtime'] = None, None

                click.echo("Loading detections from " + pipe_params["det_dir"] + pipe_params['out_label'] + '*.dets.json.gz\n')
                det_data = data_io._load_dets_json(pipe_params["det_dir"] + pipe_params['out_label'] + "*.dets.json.gz")

                det_dicts = []
                for entry in det_data:
                    for det in entry["det_info"]:
                        det_dicts = det_dicts + [det]
                        if 'fk_params' in entry.keys():
                            det_dicts[-1]["wvfrm_info"] = entry["wvfrm_info"]
                            det_dicts[-1]["fk_params"] = entry["fk_params"]
                            det_dicts[-1]["det_params"] = entry["det_params"]
                det_list = [data_io._det_dict_to_likelihood(dict) for dict in det_dicts]

                infrasound._load_celerity_model(ev_bld_params['celerity_model'])
                events, ev_qls = hjl.id_events_dict(det_list, ev_bld_params, cpu_cnt, verbose=True)

                click.echo('\nIdentified ' + str(len(events)) + ' event(s):')
                for j, ev in enumerate(events):
                    dist_mat = hjl.build_distance_matrix([det_list[k] for k in ev], bm_width=ev_bld_params['back_az_width'],
                                                         rng_max=ev_bld_params['range_max'], rad_min=100.0, rad_max=(ev_bld_params['range_max'] / 4.0),
                                                         resol=ev_bld_params['resolution'],  cpu_cnt=cpu_cnt, progress=False)

                    ev_output = {'ground truth' : {}, 'det_info' : [det_dicts[k] for k in ev], 'assoc_params' : ev_bld_params, 'dist_matrix' : dist_mat, 'location' : [], 'characterization' : []}
                    click.echo('  Writing event into ' + ev_label + '-' + str(j) +'.ev.json.gz')
                    with gzip.open(ev_label + "-" + str(j) + ".ev.json.gz", 'wt', encoding='UTF-8') as zipfile:
                        json.dump(ev_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)
            else:
                click.echo("Warning!  Existing event file (" + pipe_params["ev_dir"] + pipe_params['out_label'] + ev_j_id + '-0.ev.json.gz) found.  Exiting pipeline analysis.')
                # return 0

        ######################
        ##   LOCATION AND   ##
        ## CHARACTERIZATION ##
        ######################
        if "ev_loc_cnfgs" not in pipe_params.keys():
            return

        click.echo('\n' + '=' * 50 + '\n' + '=' * 50)
        ev_files = [file for file in np.sort(os.listdir(pipe_params["ev_dir"])) if fnmatch.fnmatch(file,"*.ev.json.gz")]
        for ev_j, ev_file in enumerate(ev_files):
            click.echo('\nLoading event from ' + pipe_params["ev_dir"] + ev_file)

            ev_data = data_io._load_dets_json(pipe_params["ev_dir"] + ev_file)[0]
            det_list = [data_io._det_dict_to_likelihood(dict) for dict in ev_data['det_info']]

            loc_info = []
            for loc_config in pipe_params["ev_loc_cnfgs"]:
                loc_params, tm_lims, pgm = config.set_ev_loc_config(pipe_params["config_dir"] + loc_config)
                infrasound._load_celerity_model(loc_params["celerity_model"])
                loc_params['cpu_cnt'] = cpu_cnt

                click.echo('  Running location analysis using ' + loc_config + '\t', nl=False)
                result = cli_event.location_wrapper(det_list, loc_params, tm_lims, pgm, verbose=False)

                if result['norm'] > 0.0:
                    loc_info = loc_info + [{'params' : loc_params, 'result' : result}]

            char_info = []
            for char_config in pipe_params["ev_char_cnfgs"]:
                click.echo('  Running characterization analysis using ' + char_config + '\t', nl=False)
                char_params = config.set_ev_char_config(pipe_params["config_dir"] + char_config)

                if char_params['det_mask'] is None:
                    char_params['det_mask'] = np.ones(len(ev_data['det_info']))
                else:
                    char_params['det_mask'] = [int(val) for val in char_params['det_mask'].strip(' ()[]').split(',')]

                loc_index = np.argmin([loc['result']['NS_stdev'] * loc['result']['EW_stdev'] for loc in loc_info])
                src_loc = [loc_info[loc_index]['result']['lat_mean'],
                           loc_info[loc_index]['result']['lon_mean']]

                # extract spectra and scale to dB
                det_info2 = [dict for k, dict in enumerate(ev_data['det_info']) if bool(char_params['det_mask'][k])]
                det_list2 = [data_io._det_dict_to_likelihood(dict) for dict in det_info2]
                det_specs = [spye.extract_json_spectra(det['spec']) for det in det_info2]
                det_specs = [np.array([spec[0], 10.0 * np.log10(spec[1])]) for spec in det_specs]

                # set up TLMs
                tlm_dir = os.path.dirname(char_params['tlm_label'])
                tlm_pattern = char_params['tlm_label'].split("/")[-1]
                tlm_files = [file_name for file_name in np.sort(os.listdir(tlm_dir)) if fnmatch.fnmatch(file_name, tlm_pattern + "*")]

                models = [0] * 2
                models[0] = [float(file_name.split("Hz")[0][len(tlm_pattern):]) for file_name in tlm_files]
                models[1] = [0] * len(tlm_files)
                for n in range(len(tlm_files)):
                    models[1][n] = sp_prop.TLossModel()
                    models[1][n].load(tlm_dir + "/" + tlm_files[n], verbose=False)

                spye_result = spye.run(det_list2, det_specs, src_loc, np.array([char_params['freq_min'], char_params['freq_max']]), models,
                                    yld_rng=np.array([char_params['yld_min'] * 1.0e3, char_params['yld_max'] * 1.0e3]),
                                    ref_src_rng=char_params['ref_rng'], resol=char_params['resolution'], grnd_brst= char_params['grnd_burst'],
                                    p_amb= char_params['amb_press'], T_amb= char_params['amb_temp'], exp_type= char_params['exp_type'])

                char_info = char_info + [{'params' : char_params, 'result' : spye_result}]

            click.echo('  Writing results into event file...')
            ev_data['location'] = loc_info
            ev_data['characterization'] = char_info
            with gzip.open(pipe_params["ev_dir"] + ev_file, 'wt', encoding='UTF-8') as zipfile:
                json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

if __name__ == '__main__':
    pipeline()
