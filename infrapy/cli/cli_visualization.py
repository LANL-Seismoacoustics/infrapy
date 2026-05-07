#!/usr/bin/env python
import sys
import os
import warnings 

import click

import json
import gzip

import configparser as cnfg
import numpy as np

import matplotlib.pyplot as plt 


from scipy.stats import chi2
from obspy import UTCDateTime

from multiprocessing import Pool


from ..utils import config
from ..utils import data_io
from ..detection import visualization as det_vis
from ..location import visualization as loc_vis
from ..location import bisl


@click.command('beam', short_help="Plot detections from beamforming")
@click.option("--det-file", help="Detection GZIP file", default=None)
@click.option("--det-index", help="Index of a single detection", default=None, type=int)
@click.option("--plot-all-dets", help="Plot all detections", default=False)
@click.option("--param-index", help="Index of a parameter set (merged dets)", default=None, type=int)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def beam_detect(det_file, det_index, plot_all_dets, param_index, figure_out, show_figure):
    '''
    Visualize beam detection results

    Example usage:
    \tinfrapy plot beam --det-file data/YJ.BRP_2012.04.09T18.00.00.dets.json.gz
    \tinfrapy plot beam --det-file data/YJ.BRP_2012.04.09T18.00.00.dets.json.gz  --det-index 2
    '''

    click.echo("")
    click.echo("###############################")
    click.echo("##                           ##")
    click.echo("##          InfraPy          ##")
    click.echo("##       Beam Detection      ##")
    click.echo("##       Visualization       ##")
    click.echo("##                           ##")
    click.echo("###############################")
    click.echo("")    

    if os.path.splitext(det_file)[-1] == ".gz":
        det_data = json.load(gzip.open(det_file, 'rt'))
    else:
        det_data = json.load(open(det_file))
    
    if "wvfrm_info" in det_data:
        click.echo('\n' + "waveform summary:")
        for wvfrm in det_data['wvfrm_info'][0]:
            print('\t' + wvfrm['trace id'], end="")
            print('\t' + wvfrm['starttime'] + ' - ' + wvfrm['endtime'])

        click.echo('\n' + "fk (beam) parameters:")
        for key in det_data['fk_params'][0].keys():
            click.echo("  " + key + ": " + str(det_data['fk_params'][0][key]))

        click.echo('\n' + "detection parameters:")
        for key in det_data['det_params'][0].keys():
            click.echo("  " + key + ": " + str(det_data['det_params'][0][key]))

        if det_index is not None:
            if det_index > len(det_data["det_info"]):
                click.echo('\n' + "detection index (" + str(det_index) + ") doesn't correspond to a detection in this file.")

                click.echo('\n' + "detection summary:")
                for nd, det in enumerate(det_data['det_info']):
                    print("   index: " + str(nd), end='\t')
                    print("   time: " + det['peak f-stat time'], end='\t')
                    print("   f-stat: " + str(np.round(det['f-stat'], 1)), end='\t')
                    print("   back azimuth: " + str(np.round(det['back az'], 1)), end='\t')
                    print("   trace velocity: " + str(np.round(det['tr vel'], 1)), end='\t')
                    print("   duration: " + str(det['start/end'][0][-1] - det['start/end'][0][0]))
                print("")

            else:
                det = det_data['det_info'][det_index]

                click.echo('\n' + "detection summary (index = " + str(det_index) + "):")
                print("   time: " + det['peak f-stat time'])
                print("   f-stat: " + str(np.round(det['f-stat'], 1)))
                print("   back azimuth: " + str(np.round(det['back az'], 1)) + " deg (rel. N)")
                print("   trace velocity: " + str(np.round(det['tr vel'], 1)) + "m/s")
                print("   duration: " + str(det['start/end'][0][-1] - det['start/end'][0][0]) + " sec")

                click.echo('\n' + "Plotting detection index " + str(det_index))
                det_vis.plot_det_json(det_data, det_index, output_path=figure_out, show_fig=show_figure)

        else:
            click.echo('\n' + "detection summary:")
            for nd, det in enumerate(det_data['det_info']):
                print("   index: " + str(nd), end='\t')
                print("   time: " + det['peak f-stat time'], end='\t')
                print("   f-stat: " + str(np.round(det['f-stat'], 1)), end='\t')
                print("   back azimuth: " + str(np.round(det['back az'], 1)), end='\t')
                print("   trace velocity: " + str(np.round(det['tr vel'], 1)), end='\t')                
                print("   duration: " + str(det['start/end'][0][-1] - det['start/end'][0][0]))
            if plot_all_dets:
                click.echo('\n' + "Plotting all detections...")
                for k in range(len(det_data['det_info'])):
                    det_vis.plot_det_json(det_data, k, output_path=figure_out + "_det" + str(k), show_fig=False)
                if show_figure:
                    plt.show()

            else:
                click.echo('\n' + "Plotting full set of fk results and windowed detections...")
                det_vis.plot_fk_json(det_data, output_path=figure_out, show_fig=show_figure)

    else: 
        if det_index is None:
            det_index = 0

        if plot_all_dets:
            click.echo('\n' + "Plotting all detections...")
            for k in range(len(det_data['det_info'])):
                click.echo('\n' + "detection summary (index = " + str(k) + "):")

                det = det_data["det_info"][k]

                click.echo('\nRun index:' + ''.join(['\t\t\t' + str(j) for j in np.arange(len(det["fk_params"]))]))
                click.echo('-' * 32 + '-' * 24 * len(det["fk_params"]) )
                click.echo('Frequency band [Hz]:' + ''.join(['\t\t' + str(fk_j["freq_min"]) + " - " + str(fk_j["freq_max"]) for fk_j in det['fk_params']]))
                click.echo('Windows (len, step, sub) [s]:' + ''.join(['\t' + str(fk_j["window_len"]) + ", " + str(fk_j["window_step"]) + ", " + str(fk_j["sub_window_len"]) + '\t' for fk_j in det['fk_params']]))
                click.echo('Back Azimuth Grid [deg]:' + ''.join(['\t' + str(fk_j["back_az_min"]) + ", " + str(fk_j["back_az_max"]) + ", " + str(fk_j["back_az_step"]) for fk_j in det['fk_params']]))
                click.echo('Trace Vel. Grid [m/s]:\t' + ''.join(['\t' + str(fk_j["trace_vel_min"]) + ", " + str(fk_j["trace_vel_max"]) + ", " + str(fk_j["trace_vel_step"]) for fk_j in det['fk_params']]))
            

                det_vis.plot_det_json(det_data, k, output_path=figure_out + "_det" + str(k), show_fig=False)
            if show_figure:
                plt.show()

        else:
            click.echo('\n' + "detection summary (index = " + str(det_index) + "):")

            det = det_data["det_info"][det_index]

            click.echo('\nRun index:' + ''.join(['\t\t\t' + str(j) for j in np.arange(len(det["fk_params"]))]))
            click.echo('-' * 32 + '-' * 24 * len(det["fk_params"]) )
            click.echo('Frequency band [Hz]:' + ''.join(['\t\t' + str(fk_j["freq_min"]) + " - " + str(fk_j["freq_max"]) for fk_j in det['fk_params']]))
            click.echo('Windows (len, step, sub) [s]:' + ''.join(['\t' + str(fk_j["window_len"]) + ", " + str(fk_j["window_step"]) + ", " + str(fk_j["sub_window_len"]) + '\t' for fk_j in det['fk_params']]))
            click.echo('Back Azimuth Grid [deg]:' + ''.join(['\t' + str(fk_j["back_az_min"]) + ", " + str(fk_j["back_az_max"]) + ", " + str(fk_j["back_az_step"]) for fk_j in det['fk_params']]))
            click.echo('Trace Vel. Grid [m/s]:\t' + ''.join(['\t' + str(fk_j["trace_vel_min"]) + ", " + str(fk_j["trace_vel_max"]) + ", " + str(fk_j["trace_vel_step"]) for fk_j in det['fk_params']]))
            
            click.echo('\n' + "Plotting detection...")
            det_vis.plot_det_json(det_data, det_index, param_index, output_path=figure_out, show_fig=show_figure)



@click.command('wvfrms', short_help="Plot waveform from detection or event")
@click.option("--det-file", help="Detection GZIP file", default=None)
@click.option("--event-file", help="Event GZIP JSON file", default=None)
@click.option("--det-index", help="Index of a single detection", default=0)
@click.option("--plot-all-dets", help="Plot waveforms for all detections", default=False)
@click.option("--use-loc", help="Use location solution to defined ranges", default=False)
@click.option("--use-gt", help="Use ground truth to defined ranges", default=False)
@click.option("--event-loc-index", help="Index of location for event result", default=0)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def wvfrms(det_file, event_file, det_index, plot_all_dets, use_loc, use_gt, event_loc_index, figure_out, show_figure):
    '''
    Summarize the contents of a JSON detections file

    Example usage (requires 'infrapy run_fkd --config-file config/detection_local.config' run first):
    \tinfrapy plot wvfrms --det-file data/YJ.BRP_2012.04.09_18.00.00-18.19.59.dets.json.gz --det-index 0
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##     Waveform Visualization      ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")  

    if det_file is not None:
        if os.path.splitext(det_file)[-1] == ".gz":
            det_data = json.load(gzip.open(det_file, 'rt'))
        else:
            det_data = json.load(open(det_file))

        if plot_all_dets:
            click.echo('\n' + "Plotting all detections...")
            for k in range(len(det_data['det_info'])):
                det = det_data["det_info"][k]

                click.echo('\nRun index:' + ''.join(['\t\t\t' + str(j) for j in np.arange(len(det["fk_params"]))]))
                click.echo('-' * 32 + '-' * 24 * len(det["fk_params"]) )
                click.echo('Frequency band [Hz]:' + ''.join(['\t\t' + str(fk_j["freq_min"]) + " - " + str(fk_j["freq_max"]) for fk_j in det['fk_params']]))
                click.echo('Windows (len, step, sub) [s]:' + ''.join(['\t' + str(fk_j["window_len"]) + ", " + str(fk_j["window_step"]) + ", " + str(fk_j["sub_window_len"]) + '\t' for fk_j in det['fk_params']]))
                click.echo('Back Azimuth Grid [deg]:' + ''.join(['\t' + str(fk_j["back_az_min"]) + ", " + str(fk_j["back_az_max"]) + ", " + str(fk_j["back_az_step"]) for fk_j in det['fk_params']]))
                click.echo('Trace Vel. Grid [m/s]:\t' + ''.join(['\t' + str(fk_j["trace_vel_min"]) + ", " + str(fk_j["trace_vel_max"]) + ", " + str(fk_j["trace_vel_step"]) for fk_j in det['fk_params']]))

                if figure_out is not None:
                    figure_out = figure_out + "_det" + str(k)

                det_vis.plot_wvfrms(det, output_path=figure_out, show_fig=False)

            if show_figure:
                plt.show()
        else:
            det = det_data["det_info"][det_index]

            click.echo('\nRun index:' + ''.join(['\t\t\t' + str(j) for j in np.arange(len(det["fk_params"]))]))
            click.echo('-' * 32 + '-' * 24 * len(det["fk_params"]) )
            click.echo('Frequency band [Hz]:' + ''.join(['\t\t' + str(fk_j["freq_min"]) + " - " + str(fk_j["freq_max"]) for fk_j in det['fk_params']]))
            click.echo('Windows (len, step, sub) [s]:' + ''.join(['\t' + str(fk_j["window_len"]) + ", " + str(fk_j["window_step"]) + ", " + str(fk_j["sub_window_len"]) + '\t' for fk_j in det['fk_params']]))
            click.echo('Back Azimuth Grid [deg]:' + ''.join(['\t' + str(fk_j["back_az_min"]) + ", " + str(fk_j["back_az_max"]) + ", " + str(fk_j["back_az_step"]) for fk_j in det['fk_params']]))
            click.echo('Trace Vel. Grid [m/s]:\t' + ''.join(['\t' + str(fk_j["trace_vel_min"]) + ", " + str(fk_j["trace_vel_max"]) + ", " + str(fk_j["trace_vel_step"]) for fk_j in det['fk_params']]))

            det_vis.plot_wvfrms(det, output_path=figure_out, show_fig=show_figure)

    elif event_file is not None:

        ev_data = json.load(gzip.open(event_file, 'rt'))

        click.echo('\n' + "=" * 17 + '\n' + "Detection Summary" + '\n' + "=" * 17 + '\n')
        for det in ev_data['det_info']:
            click.echo(det['wvfrm_info'][0][0]['trace id'])
            click.echo("  location: " + str(det['wvfrm_info'][0][0]['latitude']) + ", " + str(det['wvfrm_info'][0][0]['longitude']))
            click.echo("  detection time: " + det['peak f-stat time'])
            click.echo("  back azimuth [deg]: " + str(np.round(det["back az"], 2)))
            click.echo("  tface velocity [m/s]: " + str(np.round(det["tr vel"], 2)))
            click.echo("  f-stat: " + str(np.round(det["f-stat"], 2)))
            click.echo("")

        if len(ev_data['location']) > 0:
            click.echo('\n' + "=" * 20 + '\n' + "Localization Summary" + '\n' + "=" * 20)
            for loc_k, loc in enumerate(ev_data['location']):
                click.echo('\n' + "#" * 14)
                click.echo("## " + "index: " + str(loc_k) + " ##")
                click.echo("#" * 14)

                click.echo("parameters" + '\n' + "-" * 10)
                for key in loc['params'].keys():
                    if loc['params'][key] is not None:
                        click.echo("    " + key + ": " + str(loc['params'][key]))

                lat = str(np.round(loc['result']['lat_mean'], 3))
                lon = str(np.round(loc['result']['lon_mean'], 3))
                NS_std = str(np.round(loc['result']['NS_stdev'], 2))
                EW_std = str(np.round(loc['result']['EW_stdev'], 2))
                tm_std = str(np.round(loc['result']['t_stdev'], 1))

                click.echo('\n' + "result" + '\n' + "-" * 6)
                click.echo("    latitude: " + lat + " deg +/- " + NS_std + " km.")
                click.echo("    longitude: " + lon + " deg +/- " + EW_std + " km.")
                click.echo("    origin time: " + loc['result']['t_mean'] + " +/- " + tm_std + " s." + '\n')

        if len(ev_data['ground truth'].keys()) > 0:
            click.echo('\n' + "=" * 20 + '\n' + "Ground Truth Summary"  + '\n' + "=" * 20)
            for key in ev_data['ground truth']:
                click.echo("  " + key + ': ' + str(ev_data['ground truth'][key]))
            click.echo("")
        
        loc_vis.plot_ev_wvfrms(ev_data, use_loc=use_loc, loc_index=event_loc_index, use_gt=use_gt)
           

@click.command('spectral', short_help="Visualize detection(s) from spectral analysis")
@click.option("--det-file", help="Detection GZIP file", default=None)
@click.option("--det-index", help="Index of a single detection", default=None, type=int)
@click.option("--log-scale-freq", help="Visualize frequency in log scaling", default=False)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def spec_detect(det_file, log_scale_freq, det_index, figure_out, show_figure):
    '''
    Visualize spectral detection (sd) results

    \b
    Example usage (run from infrapy/examples directory after running fd examples or fkd examples):
    \tinfrapy plot spectral --det-file 'data/YJ.BRP1_2012.04.09T18.00.00.dets.json.gz'
    \tinfrapy plot spectral --det-file 'data/YJ.BRP1_2012.04.09T18.00.00.dets.json.gz' --det-index 3

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##     Spectral Detection (sd)     ##")
    click.echo("##          Visualization          ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if os.path.splitext(det_file)[-1] == ".gz":
        det_data = json.load(gzip.open(det_file, 'rt'))
    else:
        det_data = json.load(open(det_file))

    click.echo('\n' + "waveform summary:")
    for wvfrm in det_data['wvfrm_info']:
        print('\t' + wvfrm['trace id'], end="")
        print('\t' + wvfrm['starttime'] + ' - ' + wvfrm['endtime'])

    click.echo('\n' + "sd (spectral detector) parameters:")
    for key in det_data['sd_params'].keys():
        click.echo("  " + key + ": " + str(det_data['sd_params'][key]))

    if det_index is not None:
        if det_index > len(det_data["det_info"]):
            click.echo('\n' + "detection index (" + str(det_index) + ") doesn't correspond to a detection in this file.")

            click.echo('\n' + "detection summary:")
            for nd, det in enumerate(det_data['det_info']):
                print("   index: " + str(nd), end='\t')
                print("   time: " + det['peak f-stat time'], end='\t')
                print("   f-stat: " + str(np.round(det['f-stat'], 1)), end='\t')
                print("   back azimuth: " + str(np.round(det['back az'], 1)), end='\t')
                print("   trace velocity: " + str(np.round(det['tr vel'], 1)), end='\t')
                print("   duration: " + str(det['time vals'][-1] - det['time vals'][0]))
            print("")
        else:

            det = det_data['det_info'][det_index]

            spec_pnts = np.array(det['spec pnts'])
            t1, t2 = min(spec_pnts[:, 0]), max(spec_pnts[:, 0])
            f1, f2 = min(spec_pnts[:, 1]), max(spec_pnts[:, 1])

            print('\n' + "detection summary (index = " + str(det_index) + "):")
            print("   time: " + det['peak f-stat time'])
            print("   duration [s]: " + str(np.round(t2 - t1, 2)))
            print("   frequency range [Hz]: " + str(f1) + " - " + str(f2))
        
            click.echo('\n' + "Plotting detection index " + str(det_index))
            det_vis.plot_sd_single_json(det_data, det_index, log_scale_freq=log_scale_freq, output_path=figure_out, show_fig=show_figure)

    else:
        click.echo('\n' + "detection summary:")
        for nd, det in enumerate(det_data['det_info']):
            spec_pnts = np.array(det['spec pnts'])
            t1, t2 = min(spec_pnts[:, 0]), max(spec_pnts[:, 0])
            f1, f2 = min(spec_pnts[:, 1]), max(spec_pnts[:, 1])

            print("   index: " + str(nd), end='\t')
            print("   time: " + det['peak f-stat time'], end='\t')
            print("   duration [s]: " + str(np.round(t2 - t1, 2)), end='\t')
            print("   frequency range [Hz]: " + str(f1) + ", " + str(f2))

        click.echo('\n' + "Plotting spectrogram detection results...")
        det_vis.plot_sd_json(det_data, log_scale_freq=log_scale_freq, output_path=figure_out, show_fig=show_figure)



@click.command('map_dets', short_help="Plot detections on a map")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--detect-files", help="Detection path and pattern (option 1)", default=None)
@click.option("--event-file", help="EVent file (option 2)", default=None)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def map_dets(config_file, detect_files, event_file, range_max, figure_out, offline_maps_dir, show_figure):
    '''
    Visualize detections on a map

    \b
    Example usage (run from infrapy/examples directory after running run_assoc example):
    \tinfrapy plot map_dets --detect-files 'data/Blom_etal2020_GJI/SY*' --range-max 1500
    \tinfrapy plot map_dets --event-file data/Blom_etal2020_GJI/Blom_etal2020_GJI-0.ev.json.gz 

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##          Detection List         ##")
    click.echo("##             Mapping             ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    


    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    detect_files = config.set_param(user_config, 'DETECTION IO', 'detect_files', detect_files, 'string')
    event_file = config.set_param(user_config, 'EVENT IO', 'event_file', event_file, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  detect_files: " + str(detect_files))
    click.echo("  event_file: " + str(event_file))

    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    offline_maps_dir = config.set_param(user_config, 'VISUALIZATION', 'offline_maps_dir', offline_maps_dir, 'string')

    click.echo('\n' + "Visualization parameters:")
    click.echo("  range_max: " + str(range_max) + '\n')
    if offline_maps_dir:
        click.echo("  offline maps directory: {}".format(offline_maps_dir))
        loc_vis.use_offline_maps(offline_maps_dir)

    if detect_files is not None:
        det_data = data_io._load_dets_json(detect_files)
        det_dicts = []
        for entry in det_data:
            for det in entry["det_info"]:
                det_dicts = det_dicts + [det]
                if 'fk_params' in entry.keys():
                    det_dicts[-1]["wvfrm_info"] = entry["wvfrm_info"]
                    det_dicts[-1]["fk_params"] = entry["fk_params"]
                    det_dicts[-1]["det_params"] = entry["det_params"]

    elif event_file is not None:
        ev_info = data_io._load_dets_json(event_file)[0]
        det_dicts = ev_info["det_info"]
        range_max = ev_info["assoc_params"]["range_max"]
        click.echo("Updating range max from event building parameters: " + str(range_max) + " km")
        
    else:
        click.echo("Requires either detection file(s) or event file to plot detection projections)")
        return 

    det_list = [data_io._det_dict_to_likelihood(dict) for dict in det_dicts]

    click.echo('\n' + "Drawing map with detection back azimuth projections...")
    loc_vis.plot_dets_on_map(det_list, range_max=range_max, output_path=figure_out, show_fig=show_figure)


@click.command('localize', short_help="Plot localization result on a map")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--event-file", help="Detection path and pattern", default=None)
@click.option("--event-loc-index", help="Index of location results", default=None, type=int)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--confidence-level", help="Confidence level (default 90%)", default=90.0)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Generate figure on screeen", default=True)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
def ev_loc(config_file, event_file, event_loc_index, range_max, confidence_level, figure_out, show_figure, offline_maps_dir):
    '''
    Visualize BISL results in with wide or zoomed format

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot loc --local-detect-label GJI_example-ev0 --local-loc-label GJI_example-ev0 --range-max 1200.0
    \tinfrapy plot loc --local-detect-label GJI_example-ev0 --local-loc-label GJI_example-ev0 --zoom true

    '''

    click.echo("")
    click.echo("####################################")
    click.echo("##                                ##")
    click.echo("##             InfraPy            ##")
    click.echo("##   Localization visualization   ##")
    click.echo("##                                ##")
    click.echo("####################################")
    click.echo("")  

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    if offline_maps_dir is not None:
        loc_vis.use_offline_maps(offline_maps_dir)

    event_file = config.set_param(user_config, 'DATA IO', 'event_file', event_file, 'string')
    click.echo("Loading event information from event_file: " + str(event_file))
    ev_data = data_io._load_dets_json(event_file)[0]
    loc_info = ev_data['location']

    if range_max is None:
        range_max = ev_data['assoc_params']['range_max']
    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')

    if event_loc_index is not None:
        loc = ev_data['location'][event_loc_index]
        click.echo("Visualizing with event_loc_index: " + str(event_loc_index) + '\n')
        click.echo("localization parameters:")
        for key in loc['params'].keys():
            if loc['params'][key] is not None:
                click.echo("  " + key + ": " + str(loc['params'][key]))
        click.echo("")

    else:
        click.echo("  " + str(len(loc_info)) + " localization results in file")
        for loc_k, loc in enumerate(ev_data['location']):
            click.echo('\n' + "#" * 29)
            click.echo("##  " + "localization index: " + str(loc_k) + "  ##")
            click.echo("#" * 29)

            click.echo("parameters" + '\n' + "-" * 10)
            for key in loc['params'].keys():
                if loc['params'][key] is not None:
                    click.echo("  " + key + ": " + str(loc['params'][key]))

            lat = str(np.round(loc['result']['lat_mean'], 3))
            lon = str(np.round(loc['result']['lon_mean'], 3))
            NS_std = str(np.round(loc['result']['NS_stdev'], 2))
            EW_std = str(np.round(loc['result']['EW_stdev'], 2))
            tm_std = str(np.round(loc['result']['t_stdev'], 1))

            click.echo('\n' + "result" + '\n' + "-" * 6)
            click.echo("  Latitude: " + lat + " deg +/- " + NS_std + " km.")
            click.echo("  Longitude: " + lon + " deg +/- " + EW_std + " km.")
            click.echo("  90% confidence area: " + str(np.round(np.pi * loc['result']['NS_stdev'] * loc['result']['EW_stdev'] * chi2(2).ppf(0.9), 1)) + " sqr km" )

            click.echo("  Origin time: " + loc['result']['t_mean'] + " +/- " + tm_std + " s.")

        click.echo('\n' + "Visualizing index 0 result" + '\n')
        loc = ev_data['location'][0]

    click.echo("Localization Result Summary")
    click.echo("-" * 27)
    click.echo(bisl.summarize(loc['result'], confidence_level=float(confidence_level)))

    det_list = [data_io._det_dict_to_likelihood(dict) for dict in ev_data["det_info"]]
    loc_vis.plot_localization(det_list, loc, ev_data["ground truth"], range_max=range_max, confidence_level=confidence_level, output_path=figure_out, show_fig=show_figure)


@click.command('characterize', short_help="Plot characterization result for an event")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--event-file", help="Detection path and pattern", default=None)
@click.option("--loc-index", help="Index of location results", default=None, type=int)
@click.option("--char-index", help="Index of characterization results", default=None, type=int)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--confidence-level", help="Confidence level (default 90%)", default=90.0, type=float)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Generate figure on screeen", default=True)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
def ev_char(config_file, event_file, loc_index, char_index, range_max, confidence_level, figure_out, show_figure, offline_maps_dir):
    '''
    Visualize characterization results for an event

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot characterize --event-file HRR-6.ev.json.gz

    '''

    click.echo("")
    click.echo("##################################")
    click.echo("##                              ##")
    click.echo("##            InfraPy           ##")
    click.echo("##       Characterization       ##")
    click.echo("##        Visualization         ##")
    click.echo("##                              ##")
    click.echo("##################################")
    click.echo("")  

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    if offline_maps_dir is not None:
        loc_vis.use_offline_maps(offline_maps_dir)

    event_file = config.set_param(user_config, 'DATA IO', 'event_file', event_file, 'string')
    click.echo("Loading event information from event_file: " + str(event_file))
    ev_data = data_io._load_dets_json(event_file)[0]
    loc_info = ev_data['location']

    if range_max is None:
        range_max = ev_data['assoc_params']['range_max']
    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')

    if loc_index is not None:
        loc = ev_data['location'][loc_index]
        click.echo("Visualizing with event_loc_index: " + str(loc_index) + '\n')
        click.echo("localization parameters:")
        for key in loc['params'].keys():
            if loc['params'][key] is not None:
                click.echo("  " + key + ": " + str(loc['params'][key]))
        click.echo("")

    else:
        click.echo("  " + str(len(loc_info)) + " localization results in file")
        for loc_k, loc in enumerate(ev_data['location']):
            click.echo('\n' + "#" * 29)
            click.echo("##  " + "localization index: " + str(loc_k) + "  ##")
            click.echo("#" * 29)

            click.echo("parameters" + '\n' + "-" * 10)
            for key in loc['params'].keys():
                if loc['params'][key] is not None:
                    click.echo("  " + key + ": " + str(loc['params'][key]))

            lat = str(np.round(loc['result']['lat_mean'], 3))
            lon = str(np.round(loc['result']['lon_mean'], 3))
            NS_std = str(np.round(loc['result']['NS_stdev'], 2))
            EW_std = str(np.round(loc['result']['EW_stdev'], 2))
            tm_std = str(np.round(loc['result']['t_stdev'], 1))

            click.echo('\n' + "result" + '\n' + "-" * 6)
            click.echo("  Latitude: " + lat + " deg +/- " + NS_std + " km.")
            click.echo("  Longitude: " + lon + " deg +/- " + EW_std + " km.")
            click.echo("  90% confidence area: " + str(np.round(np.pi * loc['result']['NS_stdev'] * loc['result']['EW_stdev'] * chi2(2).ppf(0.9), 1)) + " sqr km" )

            click.echo("  Origin time: " + loc['result']['t_mean'] + " +/- " + tm_std + " s.")

        click.echo('\n' + "Visualizing index 0 result" + '\n')
        loc = ev_data['location'][0]



    if char_index is not None:
        char = ev_data['characterization'][char_index]
        click.echo("Visualizing with event_loc_index: " + str(loc_index) + '\n')
        click.echo("characterization parameters:")
        for key in char['params'].keys():
            if char['params'][key] is not None:
                click.echo("  " + key + ": " + str(char['params'][key]))
        click.echo("")

    else:
        click.echo("  " + str(len(loc_info)) + " characterization results in file")
        for char_k, char in enumerate(ev_data['characterization']):
            click.echo('\n' + "#" * 34)
            click.echo("##  " + "characterization index: " + str(char_k) + "  ##")
            click.echo("#" * 34)

            click.echo("parameters" + '\n' + "-" * 10)
            for key in char['params'].keys():
                if char['params'][key] is not None:
                    click.echo("  " + key + ": " + str(char['params'][key]))

            lat = str(np.round(loc['result']['lat_mean'], 3))
            lon = str(np.round(loc['result']['lon_mean'], 3))
            NS_std = str(np.round(loc['result']['NS_stdev'], 2))
            EW_std = str(np.round(loc['result']['EW_stdev'], 2))
            tm_std = str(np.round(loc['result']['t_stdev'], 1))

            click.echo('\n' + "result" + '\n' + "-" * 6)
            click.echo("  Latitude: " + lat + " deg +/- " + NS_std + " km.")
            click.echo("  Longitude: " + lon + " deg +/- " + EW_std + " km.")
            click.echo("  90% confidence area: " + str(np.round(np.pi * loc['result']['NS_stdev'] * loc['result']['EW_stdev'] * chi2(2).ppf(0.9), 1)) + " sqr km" )

            click.echo("  Origin time: " + loc['result']['t_mean'] + " +/- " + tm_std + " s.")

        click.echo('\n' + "Visualizing index 0 result" + '\n')
        char = ev_data['characterization'][0]

    click.echo("Localization Result Summary")
    click.echo("-" * 27)
    click.echo(bisl.summarize(loc['result'], confidence_level=float(confidence_level)))
    
    loc_vis.plot_characterization(ev_data["det_info"], loc, char, ev_data["ground truth"], range_max=range_max, confidence_level=float(confidence_level), output_path=figure_out, show_fig=show_figure)















##########################################
## THE REST OF THESE ARE DEPRECATED AND ## 
##  WILL BE REMOVED IN A FUTURE UPDATE  ##
##########################################





@click.command('fk', short_help="Visualize beamforming (fk) results", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--local-fk-label", help="Local beamforming (fk) data files", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def fk(config_file, local_wvfrms, local_latlon, fdsn, db_config, network, station, location, 
    channel, starttime, endtime, freq_min, freq_max, local_fk_label, figure_out, show_figure):
    '''
    Visualize beamforming (fk) results

    \b
    Example usage (run from infrapy/examples directory after running the run_fk examples):
    \tinfrapy plot fk --local-wvfrms 'data/YJ.BRP*'
    \tinfrapy plot fk --config-file config/detection_local.config
    \tinfrapy plot fk --config-file config/detection_fdsn.config --figure-out FDSN_fk-results.png --show-figure False

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##         Beamforming (fk)        ##")
    click.echo("##          Visualization          ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Database configuration and info   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

    # Local waveform IO parameters
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    local_latlon = config.set_param(user_config, 'WAVEFORM IO', 'local_latlon', local_latlon, 'string')

    # FDSN waveform IO parameters
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    # Frequency limits
    freq_min = config.set_param(user_config, 'FK', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'FK', 'freq_max', freq_max, 'float')

    # Result IO
    local_fk_label = config.set_param(user_config, 'DETECTION IO', 'local_fk_label', local_fk_label, 'string')

    click.echo('\n' + "Data parameters:")
    if local_wvfrms is not None:
        click.echo("  local_wvfrms: " + str(local_wvfrms))
        click.echo("  local_latlon: " + str(local_latlon))
    elif fdsn is not None:
        click.echo("  fdsn: " + str(fdsn))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    elif db_config is not None:
        db_info = cnfg.ConfigParser()
        db_info.read(db_config)
        click.echo("  db_config: " + str(db_config))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
        
    click.echo("  local_fk_label: " + str(local_fk_label))

    click.echo('\n' + "Visualization parameters:")
    click.echo("  freq_min: " + str(freq_min))
    click.echo("  freq_max: " + str(freq_max))
    if figure_out:
        click.echo("  figure_out: " + figure_out)

    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if waveform data is specified and populate obspy Stream
    if stream is not None:
        stream.filter("bandpass", freqmin=freq_min, freqmax=freq_max)

        if local_fk_label is None or local_fk_label == "auto":
            if local_wvfrms is not None and "/" in local_wvfrms:
                local_fk_label = os.path.dirname(local_wvfrms) + "/"
            else:
                local_fk_label = ""
            local_fk_label = local_fk_label + data_io.stream_label(stream)

        if ".fk_results.dat" not in local_fk_label:
            local_fk_label = local_fk_label + ".fk_results.dat"

        temp = np.loadtxt(local_fk_label)
        dt, beam_peaks = temp[:, 0], temp[:, 1:]

        temp = open(local_fk_label, 'r')
        for line in temp:
            if "t0:" in line:
                t0 = np.datetime64(line.split(' ')[-1][:-1])
            elif "freq_min" in line:
                freq_min = float(line.split(' ')[-1])
            elif "freq_max" in line:
                freq_max = float(line.split(' ')[-1])

        beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
        det_vis.plot_fk1(stream, latlon, beam_times, beam_peaks, title=local_fk_label, output_path=figure_out, show_fig=show_figure)
    else:
        if os.path.isfile(local_fk_label):
            temp = np.loadtxt(local_fk_label)
            dt, beam_peaks = temp[:, 0], temp[:, 1:]

            temp = open(local_fk_label, 'r')
            for line in temp:
                if "t0:" in line:
                    t0 = np.datetime64(line.split(' ')[-1][:-1])
                elif "freq_min" in line:
                    freq_min = float(line.split(' ')[-1])
                elif "freq_max" in line:
                    freq_max = float(line.split(' ')[-1])

            beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
            det_vis.plot_fk2(beam_times, beam_peaks, output_path=figure_out, show_fig=show_figure)
        else:
            msg = "Beamforming (fk) results not found.  No file: " + local_fk_label
            warnings.warn(msg)


@click.command('fd', short_help="Visualize detections from beamforming results", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--local-fk-label", help="Local beamforming (fk) data files", default=None)
@click.option("--local-detect-label", help="Local detection data files", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def fd(config_file, local_wvfrms, local_latlon, fdsn, db_config, network, station, location, channel, starttime, endtime,
    local_fk_label, detect_label, figure_out, show_figure):
    '''
    Visualize detection (fd) results

    \b
    Example usage (run from infrapy/examples directory after running fd examples or fkd examples):
    \tinfrapy plot fd --local-wvfrms 'data/YJ.BRP*'
    \tinfrapy plot fd --config-file config/detection_local.config
    \tinfrapy plot fd --config-file config/detection_fdsn.config --figure-out FDSN_fd-results.png --show-figure False

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##          Detection (fd)         ##")
    click.echo("##          Visualization          ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Database configuration and info   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

    # Local waveform IO parameters
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    local_latlon = config.set_param(user_config, 'WAVEFORM IO', 'local_latlon', local_latlon, 'string')

    # FDSN waveform IO parameters
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    # Result IO
    local_fk_label = config.set_param(user_config, 'DETECTION IO', 'local_fk_label', local_fk_label, 'string')
    detect_label = config.set_param(user_config, 'DETECTION IO', 'detect_label', detect_label, 'string')

    click.echo('\n' + "Data parameters:")
    if local_wvfrms is not None:
        click.echo("  local_wvfrms: " + str(local_wvfrms))
        click.echo("  local_latlon: " + str(local_latlon))
    elif fdsn is not None:
        click.echo("  fdsn: " + str(fdsn))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    elif db_config is not None:
        db_info = cnfg.ConfigParser()
        db_info.read(db_config)
        click.echo("  db_config: " + str(db_config))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
        
    click.echo("  local_fk_label: " + str(local_fk_label))
    click.echo("  detect_label: " + str(detect_label))

    if figure_out:
        click.echo("  figure_out: " + figure_out)

    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if waveform data is specified and populate obspy Stream
    if stream is not None:
        if local_fk_label is None or local_fk_label == "auto":
            if local_wvfrms is not None and "/" in local_wvfrms:
                local_fk_label = os.path.dirname(local_wvfrms) + "/"
            else:
                local_fk_label = ""
            local_fk_label = local_fk_label + data_io.stream_label(stream)
            
        temp = np.loadtxt(local_fk_label + ".fk_results.dat")
        dt, beam_peaks = temp[:, 0], temp[:, 1:]

        temp = open(local_fk_label + ".fk_results.dat", 'r')
        for line in temp:
            if "t0:" in line:
                t0 = np.datetime64(line.split(' ')[-1][:-1])
            elif "freq_min" in line:
                freq_min = float(line.split(' ')[-1])
            elif "freq_max" in line:
                freq_max = float(line.split(' ')[-1])

        beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
        stream.filter("bandpass", freqmin=freq_min, freqmax=freq_max)
        
        # Read in detection list
        if detect_label is None or detect_label == 'auto':
            detect_label = local_fk_label

        det_list = data_io.set_det_list(detect_label + ".dets.json", merge=True)
        if len(det_list) == 0:
            click.echo("Note: no detections found in analysis.")

        if os.path.isfile(detect_label + ".fd_thresholds.dat"):
            temp = np.loadtxt(detect_label + ".fd_thresholds.dat")
            thresh_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
            det_thresh = [thresh_times, temp[:, 1]]

        else:
            det_thresh = None

        det_vis.plot_fk1(stream, latlon, beam_times, beam_peaks, detections=det_list, title=local_fk_label, output_path=figure_out, det_thresh=det_thresh, show_fig=show_figure)
    else:
        if os.path.isfile(local_fk_label + ".fk_times.npy"):
            temp = np.loadtxt(local_fk_label + ".fk_results.dat")
            dt, beam_peaks = temp[:, 0], temp[:, 1:]

            temp = open(local_fk_label + ".fk_results.dat", 'r')
            for line in temp:
                if "t0:" in line:
                    t0 = np.datetime64(line.split(' ')[-1][:-1])
                elif "freq_min" in line:
                    freq_min = float(line.split(' ')[-1])
                elif "freq_max" in line:
                    freq_max = float(line.split(' ')[-1])

            beam_times = np.array([t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])
            stream.filter("bandpass", freqmin=freq_min, freqmax=freq_max)

            # Read in detection list
            det_list = data_io.set_det_list(detect_label, merge=True)
            if len(det_list) == 0:
                click.echo("Note: no detections found in analysis.")

            det_vis.plot_fk2(beam_times, beam_peaks, detections=det_list, output_path=figure_out, show_fig=show_figure)
        else:
            msg = "Beamforming (fk) results not found.  No file: " + local_fk_label + ".fk_times.npy"
            warnings.warn(msg)


@click.command('sd', short_help="Visualize detection(s) from spectral analysis", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--local-detect-label", help="Local detection data files", default=None)
@click.option("--spectral-option", help="Spectrogram method ('spectogram', 'stft', or 'cwt'), default: " + config.defaults['SD']['spectral_option'] + ")", default=None)
@click.option("--morlet-omega0", help="Frequency scaling for Morlet wavelet in 'cwt', default: " + config.defaults['SD']['morlet_omega0'] + ")", default=None, type=float)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['SD']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['SD']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--signal-start", help="Start of analysis window", default=None)
@click.option("--signal-end", help="End of analysis window", default=None)
@click.option("--det-index", help="Index of a single detection", default=None, type=int)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def sd(config_file, local_wvfrms, local_latlon, fdsn, db_config, network, station, location, channel, starttime, endtime,
    detect_label, spectral_option, morlet_omega0, freq_min, freq_max, signal_start, signal_end, det_index, figure_out, show_figure):
    '''
    Visualize spectral detection (sd) results

    \b
    Example usage (run from infrapy/examples directory after running fd examples or fkd examples):
    \tinfrapy plot sd --local-wvfrms 'data/YJ.BRP1..EDF.SAC'
    \tinfrapy plot sd --local-wvfrms 'data/YJ.BRP1..EDF.SAC' --spectral-option cwt --morlet-omega0 12.0
    \tinfrapy plot sd --local-wvfrms 'data/YJ.BRP1..EDF.SAC' --det-index 2

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##     Spectral Detection (sd)     ##")
    click.echo("##          Visualization          ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    # Database configuration and info   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

    # Local waveform IO parameters
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    local_latlon = config.set_param(user_config, 'WAVEFORM IO', 'local_latlon', local_latlon, 'string')

    # FDSN waveform IO parameters
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    # Result IO
    detect_label = config.set_param(user_config, 'DETECTION IO', 'detect_label', detect_label, 'string')

    click.echo('\n' + "Data parameters:")
    if local_wvfrms is not None:
        click.echo("  local_wvfrms: " + str(local_wvfrms))
        click.echo("  local_latlon: " + str(local_latlon))
    elif fdsn is not None:
        click.echo("  fdsn: " + str(fdsn))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    elif db_config is not None:
        db_info = cnfg.ConfigParser()
        db_info.read(db_config)
        click.echo("  db_config: " + str(db_config))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
        
    click.echo("  detect_label: " + str(detect_label))

    if figure_out:
        click.echo("  figure_out: " + figure_out)

    # Algorithm parameters
    spectral_option = config.set_param(user_config, 'SD', 'spectral_option', spectral_option, 'string')
    morlet_omega0 = config.set_param(user_config, 'SD', 'morlet_omega0', morlet_omega0, 'float')    
    freq_min = config.set_param(user_config, 'SD', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'SD', 'freq_max', freq_max, 'float')
    signal_start = config.set_param(user_config, 'SD', 'signal_start', signal_start, 'string')
    signal_end = config.set_param(user_config, 'SD', 'signal_end', signal_end, 'string')

    stream, _ = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    # Check if using a signal window
    if signal_start is not None:
        t1 = UTCDateTime(signal_start)
        t2 = UTCDateTime(signal_end)

        click.echo('\n' + "Trimming data to signal analysis window...")
        click.echo('\t' + "start time: " + str(t1))
        click.echo('\t' + "end time: " + str(t2))

        warning_message = "signal_start and signal_end values poorly defined."
        if t1 > t2:
            warning_message = warning_message + "  signal_start after signal_end."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        elif t1 < stream[0].stats.starttime:
            warning_message = warning_message + "  signal_start before data start time."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        elif t2 > stream[0].stats.endtime:
            warning_message = warning_message + "  signal_end after data end time."
            warning_message = warning_message + "  Stream won't be trimmed."
            warnings.warn((warning_message))
        else:
            stream.trim(t1, t2)

    if detect_label is not None:
        if ".dets.json" not in detect_label:
            detect_label = detect_label + ".dets.json"
        det_list = json.load(open(detect_label))
    else:
        if local_wvfrms is not None and "/" in local_wvfrms:
            output_id = os.path.dirname(local_wvfrms) + "/"
        else:
            output_id = ""
        output_id = output_id + data_io.stream_label(stream)

        if ".dets.json" not in output_id:
            output_id = output_id + ".dets.json"
        det_list = json.load(open(output_id))

    if len(det_list) == 0:
        click.echo("Note: no detections found in analysis.")

    if det_index is not None:
        if det_index <= len(det_list) - 1:
            click.echo("Plotting detection info for detection index (" + str(det_index) + ")..." + '\n')
            det_vis.plot_sd_single(stream[0], det_list[det_index], [freq_min, freq_max], output_path=figure_out, show_fig=show_figure)       
        else:
            click.echo("Invalid detection index (" + str(det_index) + "), only " + str(len(det_list)) + " detections in file.")
    else:
        click.echo("Plotting spectrogram with detection info..." + '\n')
        det_vis.plot_sd(stream[0], det_list, [freq_min, freq_max], spec_option=spectral_option, morlet_omega0=morlet_omega0, output_path=figure_out, show_fig=show_figure)


@click.command('dets', short_help="Plot detections on a map", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-detect-label", help="Detection path and pattern", default=None)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
def dets(config_file, range_max, detect_label, figure_out, offline_maps_dir):
    '''
    Visualize detections on a map

    \b
    Example usage (run from infrapy/examples directory after running run_assoc example):
    \tinfrapy plot dets --local-detect-label 'data/Blom_etal2020_GJI/*'
    \tinfrapy plot dets --local-detect-label 'GJI_example-ev0.dets.json'  --range-max 1000

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##          Detection List         ##")
    click.echo("##             Mapping             ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    


    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    detect_label = config.set_param(user_config, 'DETECTION IO', 'detect_label', detect_label, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  detect_label: " + str(detect_label))

    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    offline_maps_dir = config.set_param(user_config, 'VISUALIZATION', 'offline_maps_dir', offline_maps_dir, 'string')

    click.echo('\n' + "Visualization parameters:")
    click.echo("  range_max: " + str(range_max) + '\n')
    if offline_maps_dir:
        click.echo("  offline maps directory: {}".format(offline_maps_dir))
        loc_vis.use_offline_maps(offline_maps_dir)

    det_list = data_io.set_det_list(detect_label, merge=True)

    click.echo('\n' + "Drawing map with detection back azimuth projections...")
    loc_vis.plot_dets_on_map(det_list, range_max=range_max, output_path=figure_out)


@click.command('loc', short_help="Plot localization result on a map", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-detect-label", help="Detection path and pattern", default=None)
@click.option("--local-loc-label", help="Localization results path", default=None)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--zoom", help="Option to zoom in on the estimated source region", default=False)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--grnd-truth", help="Ground truth location", default=None)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
def loc(config_file, detect_label, local_loc_label, range_max, zoom, figure_out, grnd_truth, offline_maps_dir):
    '''
    Visualize BISL results in with wide or zoomed format

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot loc --local-detect-label GJI_example-ev0 --local-loc-label GJI_example-ev0 --range-max 1200.0
    \tinfrapy plot loc --local-detect-label GJI_example-ev0 --local-loc-label GJI_example-ev0 --zoom true

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##       Localization Mapping      ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")  

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    detect_label = config.set_param(user_config, 'DETECTION IO', 'detect_label', detect_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  local_event_label: " + str(detect_label))
    click.echo("  local_loc_label: " + str(local_loc_label))

    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    offline_maps_dir = config.set_param(user_config, 'VISUALIZATION', 'offline_maps_dir', offline_maps_dir, 'string')

    click.echo('\n' + "Visualization parameters:")
    click.echo("  range_max: " + str(range_max))
    click.echo("  zoom: " + str(zoom))

    if grnd_truth is not None:
        grnd_truth = [float(val) for val in grnd_truth.strip(' ()[]').split(',')]

    if offline_maps_dir:
        click.echo("  offline maps directory: {}".format(offline_maps_dir))
        loc_vis.use_offline_maps(offline_maps_dir)

    click.echo('\n' + "Reading in detection list...")
    det_list = data_io.set_det_list(detect_label, merge=False)
    if ".loc.json" in local_loc_label:
        bisl_result = json.load(open(local_loc_label))
    else:
        bisl_result = json.load(open(local_loc_label + ".loc.json"))

    click.echo('\n' + "BISL Summary:")
    click.echo(bisl.summarize(bisl_result))

    click.echo("Drawing map with BISL source location estimate...")
    loc_vis.plot_loc(det_list, bisl_result, range_max=range_max, zoom=zoom, title=None, output_path=figure_out, grnd_truth=grnd_truth)
    


@click.command('origin-time', short_help="Plot origin time distribution", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-loc-label", help="Localization results", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--grnd-truth", help="Ground truth origin time for comparison", default=None)
def origin_time(config_file, local_loc_label, figure_out, grnd_truth):
    '''
    Visualize the BISL origin time distribution

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot origin-time --local-loc-label GJI_example-ev0
    '''
    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##        Origin Time Plot         ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")  

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_event_label', local_loc_label, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  detect_label: " + str(local_loc_label))

    click.echo('\n' + "Reading in BISL results...")
    if ".loc.json" in local_loc_label:
        bisl_result = json.load(open(local_loc_label))
    else:
        bisl_result = json.load(open(local_loc_label + ".loc.json"))

    click.echo("Plotting origin time distribution...")
    loc_vis.plot_origin_time(bisl_result, output_path=figure_out, grnd_truth=grnd_truth)
    

@click.command('yield', short_help="Plot yield estimate distribution", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-yld-label", help="Yield estimate result", default=None)
@click.option("--figure-out", help="Destination for figure", default=None)
def yield_plot(config_file, local_yld_label, figure_out):

    '''
    Visualize the SpYE result

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy plot yield --local-yld-label HRR-5.yld.json

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##       Yield Estimate Plot       ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")  

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        user_config = cnfg.ConfigParser()
        user_config.read(config_file)
    else:
        user_config = None

    local_yld_label = config.set_param(user_config, 'DETECTION IO', 'local_yld_label', local_yld_label, 'string')
    figure_out = config.set_param(user_config, 'DETECTION IO', 'figure_out', figure_out, 'string')

    click.echo('\n' + "Data summary:")
    click.echo("  local_yld_label: " + str(local_yld_label))

    click.echo('\n' + "Reading in SpYE results...")
    if ".yld.json" in local_yld_label:
        spye_result = json.load(open(local_yld_label))
    else:
        spye_result = json.load(open(local_yld_label + ".yld.json"))

    click.echo('\n' + 'Results Summary (tons eq. TNT):')
    click.echo('\t' + "Maximum a Posteriori Yield: " + str(spye_result['yld_vals'][np.argmax(spye_result['yld_pdf'])]))
    click.echo('\t' + "68% Confidence Bounds: " + str(spye_result['conf_bnds'][0]))
    click.echo('\t' + "95% Confidence Bounds: " + str(spye_result['conf_bnds'][1]))

    click.echo('\n' + "Plotting yield PDF...")
    loc_vis.plot_spye(spye_result, output_path=figure_out)

    click.echo("")
