#!which python
"""
cli_utils.py

Utility methods accessible in the command line interface (CLI) of infrapy

Author: pblom@lanl.gov    
"""

import os
import pickle
import click
import json
import gzip

import warnings

import configparser as cnfg
from inspect import trace

import matplotlib.pyplot as plt 

from scipy.stats import gaussian_kde, norm
from scipy.optimize import curve_fit
from scipy.signal import hilbert 

from pathlib import Path

import numpy as np

from obspy import UTCDateTime 

from pyproj import Geod

from infrapy.detection import beamforming_new
from infrapy.propagation import likelihoods as lklhds
from infrapy.utils import config, data_io





@click.command('check-db-wvfrms', short_help="Check waveform pull from database")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--db-config", help="Database configuration file", default=None)

@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)

@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
def check_db_wvfrm(config_file, db_config, network, station, location, channel, starttime, endtime):
    '''
    Test database pull of waveform data for beamforming (fk or fdk) analysis

    \b
    Example usage (detection_db.config will be unique to your database pull):
    \tinfrapy run_fk --config-file config/detection_db.config

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##       check_db_wvfrms       ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Database and data IO parameters   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    click.echo('\n' + "Data parameters:")
    click.echo("  db_config: " + str(db_config))
    click.echo("  network: " + str(network))
    click.echo("  station: " + str(station))
    click.echo("  location: " + str(location))
    click.echo("  channel: " + str(channel))
    click.echo("  starttime: " + str(starttime))
    click.echo("  endtime: " + str(endtime))

    # Check data option and populate obspy Stream
    db_info = cnfg.ConfigParser()
    db_info.read(db_config)

    stream, latlon = data_io.set_stream(None, None, db_info, network, station, location, channel, starttime, endtime, None)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    click.echo('\nLocation info:')    
    for line in latlon:
        click.echo(str(line[0]) + '\t' +  str(line[1]))


@click.command('write-wvfrms', short_help="Save waveforms from FDSN or database")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)

@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)

@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
def write_wvfrms(config_file, db_config, fdsn, network, station, location, channel, starttime, endtime):
    '''
    Write waveform data from an FDSN or database pull into local SAC files

    \b
    Example usage (detection_db.config will be unique to your database pull):
    \tinfrapy utils write-wvfrms --config-file config/detection_fdsn.config

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         write-wvfrms        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")   

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Database and data IO parameters   
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

    # FDSN waveform IO parameters
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    # Trimming times
    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    click.echo('\n' + "Data parameters:")
    if fdsn is not None:
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

        click.echo("  db_url: " + str(db_config))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    else:
        click.echo("Invalid data parameters.  Requires fdsn or db info.")

    stream, latlon = data_io.set_stream(None, fdsn, db_info, network, station, location, channel, starttime, endtime, None)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    click.echo('\n' + "Writing waveform data to local SAC files...")
    data_io.write_stream_to_sac(stream, latlon)



# NOTE: THIS FUNCTION IS MOVING TO STOCHPROP 

@click.command('fit-celerity', short_help="Generate a GMM celerity model")
@click.option("--data-file", help="File containing celerity information", default=None)
@click.option("--cel-index", help="Column index of celerity values", default=6)
@click.option("--atten-index", help="Column index of attenuation values", default=11)
@click.option("--atten-lim", help="Attenuation limit", default=None, type=float)
def fit_celerity(data_file, cel_index, atten_index, atten_lim):
    '''
    Compute a KDE of celerity values and generate parameters for a reciprocal celerity model

    \b
    Example usage (requires a data file with celerities):
    \tinfrapy utils fit-celerity --data-file ToyAtmo.arrivals.dat

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         fit-celerity        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")   


    click.echo("  Loading data from " + data_file)
    data = np.loadtxt(data_file)
    cel_data = data[:, cel_index]

    if atten_lim is not None:
        click.echo("  Building KDE with limited arrivals (" + str(atten_lim) + " dB Sutherland & Bass attenuation limit)")
        atten_data = data[:, atten_index]
        cel_kernel = gaussian_kde(1.0 / cel_data[atten_data > atten_lim])
    else:
        click.echo("  Building KDE for all arrival celerities")
        cel_kernel = gaussian_kde(1.0 / cel_data)

    cel_vals = np.linspace(0.38, 0.18, 200)
    rcel_pdf = cel_kernel(1.0 / cel_vals)

    click.echo("  Generating fit to KDE...")
    def rcel_func(rcel, wt1, wt2, wt3, mn1, mn2, mn3, std1, std2, std3):
        result = (wt1 / std1) * norm.pdf((rcel - mn1) / std1)
        result = result + (wt2 / std2) * norm.pdf((rcel - mn2) / std2)
        result = result + (wt3 / std3) * norm.pdf((rcel - mn3) / std3)

        return result
    
    popt, _ = curve_fit(rcel_func, 1.0 / cel_vals, rcel_pdf,
                         p0=[0.0539, 0.0899, 0.8562, 
                             1.0 / 0.327, 1.0 / 0.293, 1.0 / 0.26,
                             0.066, 0.08, 0.33])
    popt = np.round(popt, 3)

    click.echo('\n' + "  Reciprocal celerity model parameters (CLI and config file formats):")
    click.echo("    --rcel-wts '" + str(popt[0]) + ", " + str(popt[1]) + ", " + str(popt[2]) + "' --rcel-mns '" + str(popt[3]) + ", " + str(popt[4]) + ", " + str(popt[5]) + "' --rcel-sds '" + str(popt[6]) + ", " + str(popt[7]) + ", " + str(popt[8]) + "'" + '\n')

    click.echo("    rcel_wts = '" + str(popt[0]) + ", " + str(popt[1]) + ", " + str(popt[2]) + "'")
    click.echo("    rcel_mns = '" + str(popt[3]) + ", " + str(popt[4]) + ", " + str(popt[5]) + "'")
    click.echo("    rcel_sds = '" + str(popt[6]) + ", " + str(popt[7]) + ", " + str(popt[8]) + "'" + '\n')

    click.echo("    Note: mean reciprocal celerities: 1.0/" + str(np.round(1.0 / popt[3], 3)) + ", 1.0/" + str(np.round(1.0 / popt[4], 3)) + ", 1.0/" + str(np.round(1.0 / popt[5], 3)) + '\n')

    plt.figure(figsize=(7, 4))
    plt.plot(cel_vals, rcel_pdf, '-k', linewidth=4.0, label="Data KDE")
    plt.plot(cel_vals, rcel_func(1.0 / cel_vals, popt[0], popt[1], popt[2],popt[3], popt[4], popt[5], 
                                 popt[6], popt[7], popt[8]), '--r', linewidth=2.0, label="GMM Fit")
    plt.xlabel("Celerity [km/s]")
    plt.ylabel("Probability")
    plt.legend()
    plt.show()




@click.command('merge-dets', short_help="Check waveform pull from database")
@click.option("--det-files", help="Detection GZIP files", default=None)
@click.option("--merged-label", help="Output detection file label", default=None)
def merge_dets(det_files, merged_label):


    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         merge_dets          ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  

    dets_data = data_io._load_dets_json(det_files)

    click.echo('\n' + "Unique fk (beam) parameters:")
    for key in dets_data[0]['fk_params'][0].keys():
        vals = [det['fk_params'][0][key] for det in dets_data]
        if None not in vals:
            vals = np.unique(vals)
            if len(vals) == 1:
                vals = vals[0]
        else:
            if all(val is None for val in vals):
                vals = None 

        click.echo("  " + key + ": " + str(vals))

    click.echo('\n' + "Unique detection parameters:")
    for key in dets_data[0]['det_params'][0].keys():
        vals = [det['det_params'][0][key] for det in dets_data]
        if None not in vals:
            vals = np.unique(vals)
            if len(vals) == 1:
                vals = vals[0]
        else:
            if all(val is None for val in vals):
                vals = None 
        click.echo("  " + key + ": " + str(vals))

    click.echo('\n' + "Merging detections...")
    det_list = []
    for entry in dets_data:
        for det in entry["det_info"]:
            det_list = det_list + [det]
            det_list[-1]["wvfrm_info"] = entry["wvfrm_info"]
            det_list[-1]["fk_params"] = entry["fk_params"]
            det_list[-1]["det_params"] = entry["det_params"]

    dets_out = []
    while len(det_list) > 0:
        merge_indices = [0]
        print("")

        for k, det_k in enumerate(det_list[1:]):
            # check at least one station ID matches
            ids_0 = [ch['trace id'] for ch in det_list[0]['wvfrm_info'][0]]
            ids_k = [ch['trace id'] for ch in det_k['wvfrm_info'][0]]

            if any(id in ids_0 for id in ids_k):
                # compute detection time overlap
                dt = abs(UTCDateTime(det_list[0]["peak f-stat time"]) - UTCDateTime(det_k["peak f-stat time"]))

                dur1 = max(60.0, det_list[0]["start/end"][0][1] - det_list[0]["start/end"][0][0])
                dur2 = max(60.0, det_k["start/end"][0][1] - det_k["start/end"][0][0])
                dt = dt / (2.0 * max(dur1, dur2))

                # check back azimuths are within tolerance 
                daz = abs(det_list[0]["back az"] - det_k["back az"])
                if daz > 360.0:
                    daz = daz - 360.0
                daz = daz / 30.0

                # print("   ", det_list[0]["peak f-stat time"], '\t', det_k["peak f-stat time"], '\t', dt, '\t', daz, '\t', np.sqrt(dt**2 + daz**2))

                if np.sqrt(dt**2 + daz**2) < 0.75:
                    merge_indices = merge_indices + [k + 1]

        dets_to_merge = [det_list[j] for j in merge_indices]
        click.echo('\n' + "Detections to merge:")
        for det in dets_to_merge:
            click.echo("  " + det["peak f-stat time"] + ", " + str(det["back az"]))

        f_stat_vals = [det["f-stat"] for det in dets_to_merge]
        tm_vals = [det["peak f-stat time"] for det in dets_to_merge]
        t0 = tm_vals[np.argmax(f_stat_vals)]

        dets_out = dets_out + [det_list[0]]
        dets_out[-1]["f-stat"] = np.max(f_stat_vals)
        dets_out[-1]["peak f-stat time"] = t0

        dt = UTCDateTime(dets_to_merge[0]["peak f-stat time"]) - UTCDateTime(t0)
        dets_out[-1]["fk"][0]["time"] = np.array(dets_out[-1]["fk"][0]["time"]) + dt
        dets_out[-1]["beam"][0]["time"] = np.array(dets_out[-1]["beam"][0]["time"]) + dt
        dets_out[-1]["start/end"][0] = np.array(dets_out[-1]["start/end"][0]) + dt

        for det in dets_to_merge[1:]:
            for key in ["wvfrm_info", "fk_params", "det_params", "start/end", "fk", "beam", "spec"]:
                dets_out[-1][key] = dets_out[-1][key] + det[key]
                
            dt = UTCDateTime(det["peak f-stat time"]) - UTCDateTime(t0)
            dets_out[-1]["start/end"][-1] = np.array(det["start/end"][-1]) + dt
            dets_out[-1]["fk"][-1]["time"] = np.array(det["fk"][-1]["time"]) + dt
            dets_out[-1]["beam"][-1]["time"] = np.array(det["beam"][-1]["time"]) + dt

        # update back azimuth and trace velocity using weighted mean...        
        az_all, tr_all, fs_all = [np.array([])] * 3
        for det in dets_to_merge:
            tm_mask = np.logical_and(det["start/end"][0][0] <= det["fk"][-1]["time"], det["fk"][-1]["time"] <= det["start/end"][0][1])
            az_all = np.append(az_all, np.array(det["fk"][-1]["back az"])[tm_mask])
            tr_all = np.append(tr_all, np.array(det["fk"][-1]["tr vel"])[tm_mask])
            fs_all = np.append(fs_all, np.array(det["fk"][-1]["f-stat"])[tm_mask])

        dets_out[-1]["back az"] = np.average(az_all, weights=fs_all)
        dets_out[-1]["tr vel"] = np.average(tr_all, weights=fs_all)

        # remove merged detections from the original list and continue
        det_list = [det_list[j] for j in range(len(det_list)) if j not in merge_indices]

    det_output = {'det_info' : dets_out}
    with gzip.open(merged_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)


@click.command('convert-dets', short_help="Convert legacy detection output to new JSON")
@click.option("--det-file", help="Detection GZIP files", default=None)
@click.option("--fk-file", help="Detection GZIP files", default=None)

@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)

@click.option("--local-latlon", help="Location information for local waveforms", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)

@click.option("--output-label", help="Output detection file label", default=None)

def convert_dets(det_file, fk_file, config_file, local_wvfrms, fdsn, db_config, local_latlon, network, station, location, channel, starttime, endtime, output_label):


    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         convert_dets        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  

    if config_file:
        if os.path.isfile(config_file):
            click.echo('\n' + "Loading configuration info from: " + config_file)
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo('\n' + "Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None

    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    local_latlon = config.set_param(user_config, 'WAVEFORM IO', 'local_latlon', local_latlon, 'string')

    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')   
    network = config.set_param(user_config, 'WAVEFORM IO', 'network', network, 'string')
    station = config.set_param(user_config, 'WAVEFORM IO', 'station', station, 'string')
    location = config.set_param(user_config, 'WAVEFORM IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'WAVEFORM IO', 'channel', channel, 'string')       

    starttime = config.set_param(user_config, 'WAVEFORM IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'WAVEFORM IO', 'endtime', endtime, 'string')

    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)
    wvfrm_info = data_io.wvfrm_info(stream, latlon)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    # Load fk parameters from the config file if one is provided (otherwise use defaults)
    fk_params = {}
    fk_params['freq_min'] = config.set_param(user_config, 'FK', 'freq_min', None, 'float')
    fk_params['freq_max'] = config.set_param(user_config, 'FK', 'freq_max', None, 'float')
    fk_params['back_az_min'] = config.set_param(user_config, 'FK', 'back_az_min', None, 'float')
    fk_params['back_az_max'] = config.set_param(user_config, 'FK', 'back_az_max', None, 'float')
    fk_params['back_az_step'] = config.set_param(user_config, 'FK', 'back_az_step', None, 'float')
    fk_params['trace_vel_min'] = config.set_param(user_config, 'FK', 'trace_vel_min', None, 'float')
    fk_params['trace_vel_max'] = config.set_param(user_config, 'FK', 'trace_vel_max', None, 'float')
    fk_params['trace_vel_step'] = config.set_param(user_config, 'FK', 'trace_vel_step', None, 'float')
    fk_params['method'] = config.set_param(user_config, 'FK', 'method', None, 'string')
    fk_params['signal_start'] = config.set_param(user_config, 'FK', 'signal_start', None, 'string')
    fk_params['signal_end'] = config.set_param(user_config, 'FK', 'signal_end', None, 'string')
    fk_params['noise_start'] = config.set_param(user_config, 'FK', 'noise_start', None, 'string')
    fk_params['noise_end'] = config.set_param(user_config, 'FK', 'noise_end', None, 'string')
    fk_params['window_len'] = config.set_param(user_config, 'FK', 'window_len', None, 'float')
    fk_params['sub_window_len'] = config.set_param(user_config, 'FK', 'sub_window_len', None, 'float')
    fk_params['window_step'] = config.set_param(user_config, 'FK', 'window_step', None, 'float')
    fk_params['cpu_cnt'] = config.set_param(user_config, 'FK', 'cpu_cnt', None, 'int')

    # Update fk values from fk_results file



    click.echo('\n' + "fk (beam) parameters:")
    for key in fk_params.keys():
        click.echo("  " + key + ": " + str(fk_params[key]))

    # Set detection parameters from config file if one is provided (otherwise use defaults)
    det_params = {}
    det_params['window_len'] = config.set_param(user_config, 'FD', 'window_len', None, 'float')
    det_params['p_value'] = config.set_param(user_config, 'FD', 'p_value', None, 'float')
    det_params['min_duration'] = config.set_param(user_config, 'FD', 'min_duration', None, 'float')
    det_params['back_az_width'] = config.set_param(user_config, 'FD', 'back_az_width', None, 'float')
    det_params['fixed_thresh'] = config.set_param(user_config, 'FD', 'fixed_thresh', None, 'float')
    det_params['thresh_ceil'] = config.set_param(user_config, 'FD', 'thresh_ceil', None, 'float')
    det_params['return_thresh'] = config.set_param(user_config, 'FD', 'return_thresh', None, 'bool')
    det_params['merge_dets'] = config.set_param(user_config, 'FD', 'merge_dets', None, 'bool')

    click.echo('\n' + "detection parameters:")
    for key in det_params.keys():
        click.echo("  " + key + ": " + str(det_params[key]))

    # Extract fk results
    fk_vals = np.loadtxt(fk_file)

    fk_out = {}
    fk_out['time'] = fk_vals[:, 0]
    fk_out['back az'] = fk_vals[:, 1]
    fk_out['tr vel'] = fk_vals[:, 2]
    fk_out['f-stat'] = fk_vals[:, 3]
    fk_out['thresh'] = np.zeros_like(fk_vals[:, 0])
    
    dets_orig = data_io._load_dets_json(det_file)

    dets_out = []
    for det in dets_orig[0]:
        click.echo('')
        click.echo(det)

        dets_out = dets_out + [{}]
        dets_out[-1]['peak f-stat time'] = det['Time (UTC)']
        dets_out[-1]['start/end'] = [[det['Start'], det['End']]]
        dets_out[-1]['f-stat'] = det['F Stat.']

        dt_ref = UTCDateTime(str(dets_out[-1]['peak f-stat time'])) - UTCDateTime(stream[0].stats.starttime)
        det_mask = np.logical_and(det['Start'] <= fk_out['time'] - dt_ref, fk_out['time'] - dt_ref <= det['End'])
        dets_out[-1]['back az'] = np.average(fk_out['back az'][det_mask], weights=fk_out['f-stat'][det_mask])
        dets_out[-1]['tr vel'] = np.average(fk_out['tr vel'][det_mask], weights=fk_out['f-stat'][det_mask])

        # Extract the fk results
        det_buffer = (det['Start'] - det['End']) * 0.15
        det_buffer = max(min(det_buffer, 60.0), 15.0)
        det_buffer = fk_params['window_step'] * np.round(det_buffer/fk_params['window_step'])
        
        det_mask = np.logical_and(det['Start'] - det_buffer <= fk_out['time'] - dt_ref,
                                    fk_out['time'] - dt_ref <= det['End'] + det_buffer)

        dets_out[-1]['fk'] = [{}]
        dets_out[-1]['fk'][0]['time'] = fk_out['time'][det_mask] - dt_ref
        dets_out[-1]['fk'][0]['back az'] = fk_out['back az'][det_mask]
        dets_out[-1]['fk'][0]['tr vel'] = fk_out['tr vel'][det_mask]
        dets_out[-1]['fk'][0]['f-stat'] = fk_out['f-stat'][det_mask]

        # Extract the beam and spectral information
        st_bm = stream.copy()

        t_ref = UTCDateTime(str(dets_out[-1]['peak f-stat time']))
        t1 = t_ref + det['Start'] - det_buffer
        t2 = t_ref + det['End'] + det_buffer

        st_bm.detrend().filter('bandpass', freqmin=fk_params['freq_min'], freqmax=fk_params['freq_max'])
        st_bm.trim(t1, t2)
        x_bm, t_bm, _, geom_bm = beamforming_new.stream_to_array_data(st_bm, latlon=latlon)
        X_bm, _, f_bm = beamforming_new.fft_array_data(x_bm, t_bm, fft_window="boxcar")

        sig_est, residual = beamforming_new.extract_signal(X_bm, f_bm, [dets_out[-1]['back az'], dets_out[-1]['tr vel']], geom_bm)

        sig_wvfrm = np.fft.irfft(sig_est)[:len(t_bm)] / (t_bm[1] - t_bm[0])
        resid_wvfrms = np.fft.irfft(residual, axis=1)[:, :len(t_bm)]  / (t_bm[1] - t_bm[0])
        resid_env = np.mean([np.abs(hilbert(resid_wvfrms[nM])) for nM in range(len(resid_wvfrms))], axis=0)
 
        dets_out[-1]['beam'] = [{}]
        dets_out[-1]['beam'][0]['time'] = t_bm + det['Start'] - det_buffer
        dets_out[-1]['beam'][0]['signal'] = sig_wvfrm
        dets_out[-1]['beam'][0]['resid'] = resid_env

        # repeat without the bandpass filter for the spectra
        st_bm2 = stream.copy()
        st_bm2.trim(t1, t2)
        st_bm2.detrend()

        x_bm2, t_bm2, _, geom_bm2 = beamforming_new.stream_to_array_data(st_bm2, latlon=latlon)
        X_bm2, _, f_bm2 = beamforming_new.fft_array_data(x_bm2, t_bm2, fft_window="boxcar")
        sig_est2, residual2 = beamforming_new.extract_signal(X_bm2, f_bm2, [dets_out[-1]['back az'], dets_out[-1]['tr vel']], geom_bm2)

        dets_out[-1]['spec'] = [{}]
        dets_out[-1]['spec'][0]['freq'] = f_bm2
        dets_out[-1]['spec'][0]['signal'] = np.abs(sig_est2)
        dets_out[-1]['spec'][0]['resid'] = np.mean(np.abs(residual2), axis=0)

    det_output = {'wvfrm_info' : [wvfrm_info], 'fk_params' : [fk_params], 'det_params' : [det_params], 'fk' : fk_out, 'det_info' : dets_out}
    with gzip.open(output_label + "-update.dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)



@click.command('event-gt', short_help="Populate event ground truth information")
@click.option("--event-file", help="Event GZIP JSON files", default=None)

@click.option("--latitude", help="event latitude (deg)", default=None)
@click.option("--longitude", help="Event longitude (deg)", default=None)
@click.option("--orig-tm", help="Origin datetime", default=None)
@click.option("--eq-tnt", help="Explosive yield (eq. TNT) [kg]", default=None)

@click.option("--user-entry", help="User specified info (comma sep. key/val)", default=None, multiple=True)
@click.option("--entry-mode", help="Add values or overwite ('append' or 'replace')", default='append')

def event_gt(event_file, latitude, longitude, orig_tm, eq_tnt, user_entry, entry_mode):


    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##     Write Event GT Info     ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  
    
    click.echo("Loading event information from event_file: " + str(event_file))
    ev_data = data_io._load_dets_json(event_file)[0]
    gt_dict = ev_data['ground truth']

    def entry_check(key, val):
        if key not in gt_dict.keys():
            gt_dict[key] = val
        else:
            if gt_dict[key] is None:
                gt_dict[key] = val
            elif entry_mode == 'replace':
                click.echo("** replacing existing entry for '" + key + "': " + gt_dict[key])            
                gt_dict[key] = val
            else:
                click.echo("Skipping '" + str(key) + "' that already has an entry.")


    base_keys = ['latitude', 'longitude', 'orig_tm', 'eq_tnt']
    base_vals = [latitude, longitude, orig_tm, eq_tnt]
    for k in range(4):
        entry_check(base_keys[k], base_vals[k])

    if gt_dict['orig_tm'] is not None:
        gt_dict['orig_tm'] = UTCDateTime(gt_dict['orig_tm'])

    for e in user_entry:
        key, val = e.split(":")
        entry_check(key, val)

    click.echo('\n' + "Ground truth summary:")
    for key in gt_dict:
        click.echo("  " + key + ': ' + str(gt_dict[key]))
    click.echo("")

    ev_data['Event GT Summary'] = gt_dict
    with gzip.open(event_file, 'wt', encoding='UTF-8') as zipfile:
        json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)


@click.command('event-summary', short_help="Summarize information in an event file")
@click.option("--event-file", help="Event GZIP JSON files", default=None)
def event_summary(event_file):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##     Summarize Event File    ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  

    click.echo("Loading information from event_file: " + str(event_file))
    ev_data = data_io._load_dets_json(event_file)[0]

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
            click.echo("    origin time: " + loc['result']['t_mean'] + " +/- " + tm_std + " s.")

    if len(ev_data['characterization']) > 0:
        click.echo('\n' + "=" * 24 + '\n' + "Characterization Summary" + '\n' + "=" * 24)
        for char_k, char in enumerate(ev_data['characterization']):
            click.echo('\n' + "#" * 14)
            click.echo("## " + "index: " + str(char_k) + " ##")
            click.echo("#" * 14)

            click.echo("parameters" + '\n' + "-" * 10)
            for key in char['params'].keys():
                if char['params'][key] is not None:
                    click.echo("    " + key + ": " + str(char['params'][key]))

            click.echo('\n' + "result" + '\n' + "-" * 6)
            click.echo("    maximum likelihood yield: " + str(np.round(char['result']['yld_vals'][np.argmax(char['result']['yld_pdf'])], 2)) + " tons eq. TNT")
            click.echo("    68% confidence bounds: " + str(char['result']['conf_bnds'][0]))
            click.echo("    95% confidence bounds: " + str(char['result']['conf_bnds'][1]))

    if len(ev_data['ground truth'].keys()) > 0:
        click.echo('\n' + "=" * 20 + '\n' + "Ground Truth Summary"  + '\n' + "=" * 20)
        for key in ev_data['ground truth']:
            click.echo("  " + key + ': ' + str(ev_data['ground truth'][key]))
        click.echo("")



##########################################
## THE REST OF THESE ARE DEPRECATED AND ## 
##  WILL BE REMOVED IN A FUTURE UPDATE  ##
##########################################


@click.command('arrivals2json', short_help="Convert infraGA/GeoAc arrivals to detection file", hidden=True)
@click.option("--arrivals-file", help="InfraGA/GeoAc arrivals file", default=None)
@click.option("--json-file", help="JSON format detection file", default=None)
@click.option("--grnd-snd-spd", help="Ground sound speed", default=340.0)
@click.option("--src-time", help="Source time", default="2020-01-01T00:00:00")
@click.option("--peakf-value", help="Fixed F-value", default=25.0)
@click.option("--array-dim", help="Array dimension", default=6)
def arrivals2json(arrivals_file, json_file, grnd_snd_spd, src_time, peakf_value, array_dim):
    '''
    Convert infraGA/GeoAc eigenray arrival results into a json detection list usable in InfraPy
    
    \b
    Example usage (requires InfraGA/GeoAc arrival output):
    \tinfrapy arrivals2json --arrivals-file example.arrivals.dat --json-file example.dets.json --grnd-snd-spd 335.0 --src-time "2020-12-25T00:00:00"

    '''
    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##        arrivals2json        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  
    
    
    click.echo("")
    click.echo("  arrivals_file: " + str(arrivals_file))
    click.echo("  json_file: " + str(json_file))

    click.echo("")
    click.echo("  grnd_snd_spd: " + str(grnd_snd_spd))
    click.echo("  src_time: " + str(src_time))
    click.echo("  peakF_value: " + str(peakf_value))
    click.echo("  array_dim: " + str(array_dim))

    arrivals = np.loadtxt(arrivals_file)

    det_list = []
    for line in arrivals:
        det = lklhds.InfrasoundDetection(lat_loc=np.round(line[3], 3), lon_loc=np.round(line[4], 3), time=(UTCDateTime(src_time) + line[5]), azimuth=np.round(line[9], 2), f_stat=peakf_value, array_d=array_dim)
        det.trace_velocity = np.round(grnd_snd_spd / np.cos(np.radians(line[8])), 1)
        det.note = "InfraGA/GeoAc arrival output"
        det_list = det_list + [det]

    data_io.detection_list_to_json(json_file, det_list)


@click.command('arrival-time', short_help="Estimate the arrival time for a source-receiver pair", hidden=True)
@click.option("--src-lat", help="Source latitude", default=None, prompt="Enter source latitude: ")
@click.option("--src-lon", help="Source longitude", default=None, prompt="Enter source longitude: ")
@click.option("--src-time", help="Source time", default=None, prompt="Enter source time: ")
@click.option("--rcvr-lat", help="Receiver latitude", default=None)
@click.option("--rcvr-lon", help="Receiver longitude", default=None)
@click.option("--rcvr", help="Reference IMS station (e.g., 'I53')", default=None)
@click.option("--celerity-min", help="Minimum celerity", default=0.24)
@click.option("--celerity-max", help="Maximum celerity", default=0.35)
def arrival_time(src_lat, src_lon, src_time, rcvr_lat, rcvr_lon, rcvr, celerity_min, celerity_max):
    '''
    Compute the range of possible arrivals times for a source-receiver pair given a range of celerity values.
    Can use a receiver latitude/longitude or reference from a list (currently only IMS stations)
    
    \b
    Example usage (requires InfraGA/GeoAc arrival output):
    \tinfrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr-lat 40.0 --rcvr-lon -110.0
    \tinfrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr I57US

    '''
    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         arrival-time        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  
    
    click.echo("  Source Time: " + src_time)
    click.echo("  Source Location: (" + str(src_lat) + ", " + str(src_lon) + ")")
    if rcvr is not None:
        click.echo('\n' + "  User specified reference receiver: " + str(rcvr))

        try:
            ims_locs_file = str(Path(__file__).parent.parent / "resources" / "IMS_infrasound_locs.pkl")
            
            with open(ims_locs_file, 'rb') as infile:
                IMS_info = pickle.load(infile, encoding='latin1')
        except FileNotFoundError:
            warning_message = " IMS locations file (IMS_infrasound_locs.pkl) is not found"
            warnings.warn(warning_message)

        for line in IMS_info:
            if rcvr in line[0]:
                click.echo("  Reference IMS station match: " + line[0])
                rcvr_lat, rcvr_lon = line[1][:2]
                click.echo("  Receiver Location: (" + str(rcvr_lat) + ", " + str(rcvr_lon) + ")")
                break
        if rcvr_lat is None:
            warning_message = "Specified reference receiver (" + rcvr + ") not found in IMS info."
            warnings.warn((warning_message))
            return 0

    elif rcvr_lat is not None and rcvr_lon is not None:
        click.echo("  Receiver Location: (" + str(rcvr_lat) + ", " + str(rcvr_lon) + ")")
    else:
        warning_message = "Method requires either a reference receiver or user defined latitude and longitude"
        warnings.warn((warning_message))
        return 0

    click.echo("")
    click.echo("  Celerity Range: (" + str(celerity_min) + ", " + str(celerity_max) + ")")

    sph_proj = Geod(ellps='sphere')
    temp = sph_proj.inv(src_lon, src_lat, rcvr_lon, rcvr_lat, radians=False)
    az, back_az = temp[0], temp[1]
    rng = temp[2] / 1000.0

    if back_az > 180.0:
        back_az = back_az - 360.0
    elif back_az < -180.0:
        back_az = back_az + 360.0

    click.echo("")
    click.echo("  Propagation range: " + str(np.round(rng,2)) + " km")
    click.echo("  Propagation azimuth: " + str(np.round(az, 2)) + " degrees" + '\n')

    click.echo("  Estimated arrival back azimuth: " + str(np.round(back_az, 2)) + " degrees")
    click.echo("  Estimated arrival time range:")
    click.echo("    " + str(UTCDateTime(src_time) + np.round(rng / celerity_max, 0))[:-8])
    click.echo("    " + str(UTCDateTime(src_time) + np.round(rng / celerity_min, 0))[:-8] + '\n')


@click.command('calc-celerity', short_help="Compute the celerity for an arrival from a known source", hidden=True)
@click.option("--src-lat", help="Source latitude", default=None, prompt="Enter source latitude: ")
@click.option("--src-lon", help="Source longitude", default=None, prompt="Enter source longitude: ")
@click.option("--src-time", help="Source time", default=None, prompt="Enter source time: ")
@click.option("--arrival-lat", help="Arrival latitude", default=None, prompt="Enter arrival latitude: ")
@click.option("--arrival-lon", help="Arrival longitude", default=None, prompt="Enter arrival longitude: ")
@click.option("--arrival-time", help="Arrival time", default=None, prompt="Enter arrival time: ")
def calc_celerity(src_lat, src_lon, src_time, arrival_lat, arrival_lon, arrival_time):
    '''
    Compute the range of possible arrivals times for a source-receiver pair given a range of celerity values
    
    \b
    Example usage (requires InfraGA/GeoAc arrival output):
    \tinfrapy utils calc-celerity --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --arrival-lat 40.0 --arrival-lon -110.0 --arrival-time "2020-12-25T01:03:50"

    '''
    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##        calc-celerity        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")  
    
    click.echo("  Source Time: " + src_time)
    click.echo("  Source Location: (" + str(src_lat) + ", " + str(src_lon) + ")")

    click.echo('\n' + "  Arrival Time: " + arrival_time)
    click.echo("  Arrival Location: (" + str(arrival_lat) + ", " + str(arrival_lon) + ")")

    dt = UTCDateTime(arrival_time) - UTCDateTime(src_time)

    sph_proj = Geod(ellps='sphere')
    temp = sph_proj.inv(src_lon, src_lat, arrival_lon, arrival_lat, radians=False)
    az = temp[0]
    rng = temp[2]

    click.echo("")
    click.echo("  Propagation time: " + str(np.round(dt, 2)) + " s")
    click.echo("  Propagation range: " + str(np.round(rng / 1000.0, 2)) + " km")
    click.echo("  Propagation azimuth: " + str(np.round(az, 2)) + " degrees")
    click.echo("  Arrival celerity: " + str(np.round(rng / dt, 1)) + " m/s" + '\n')


@click.command('best-beam', short_help="Compute the best beam via shift/stack", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-url", help="Database URL for waveform data files", default=None)
@click.option("--db-site", help="Database site table for waveform data files", default=None)
@click.option("--db-wfdisc", help="Database wfdisc table for waveform data files", default=None)
@click.option("--local-latlon", help="Array location information for local waveforms", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
@click.option("--local-fk-label", help="Label for local output of fk results", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['FK']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['FK']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--back-az", help="Back azimuth of user specified beam (degrees)", default=None, type=float)
@click.option("--trace-vel", help="Trace velocity of user specified beam (m/s))", default=None, type=float)
@click.option("--signal-start", help="Start of signal window", default=None)
@click.option("--signal-end", help="End of signal window", default=None)
@click.option("--hold-figure", help="Hold figure open", default=True)
def best_beam(config_file, local_wvfrms, fdsn, db_url, db_site, db_wfdisc, local_latlon, network, station, location, channel, starttime, endtime, local_fk_label, freq_min, freq_max,
    back_az, trace_vel, signal_start, signal_end, hold_figure):
    '''
    Shift and stack the array data to compute the best beam.  Can be run adaptively using the fk_results.dat file or along a specific beam.

    \b
    Example usage (requires 'infrapy run_fk --config-file config/detection_local.config' run first):
    \tinfrapy utils best-beam --config-file config/detection_local.config
    \tinfrapy utils best-beam --config-file config/detection_local.config --back-az -39.0 --trace-vel 358.0
    \tinfrapy utils best-beam --config-file config/detection_local.config --signal-start '2012-04-09T18:13:00' --signal-end '2012-04-09T18:15:00'

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##          best-beam          ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")   

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Database and data IO parameters   
    db_url = config.set_param(user_config, 'WAVEFORM IO', 'db_url', db_url, 'string')
    db_site = config.set_param(user_config, 'WAVEFORM IO', 'db_site', db_site, 'string')
    db_wfdisc = config.set_param(user_config, 'WAVEFORM IO', 'db_wfdisc', db_wfdisc, 'string')

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

    # Local fk file
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
    elif db_url is not None:
        click.echo("  db_url: " + str(db_url))
        click.echo("  db_site: " + str(db_site))
        click.echo("  db_wfdisc: " + str(db_wfdisc))
        click.echo("  network: " + str(network))
        click.echo("  station: " + str(station))
        click.echo("  location: " + str(location))
        click.echo("  channel: " + str(channel))
        click.echo("  starttime: " + str(starttime))
        click.echo("  endtime: " + str(endtime))
    else:
        click.echo("Invalid data parameters.  Requires fdsn or db info.")

    if local_fk_label is not None:
        click.echo("  local_fk_label: " + str(local_fk_label))

    # Algorithm parameters
    freq_min = config.set_param(user_config, 'FK', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'FK', 'freq_max', freq_max, 'float')

    signal_start = config.set_param(user_config, 'FK', 'signal_start', signal_start, 'string')
    signal_end = config.set_param(user_config, 'FK', 'signal_end', signal_end, 'string')

    click.echo('\n' + "Algorithm parameters:")
    click.echo("  freq_min: " + str(freq_min))
    click.echo("  freq_max: " + str(freq_max))
    click.echo("  signal_start: " + str(signal_start))
    click.echo("  signal_end: " + str(signal_end))
    if back_az is not None and trace_vel is not None:
        click.echo("  back_az_step: " + str(back_az))
        click.echo("  trace_vel_min: " + str(trace_vel))

    # Check data option and populate obspy Stream
    if db_url is not None:
        db_info = {'url': db_url, 'site': db_site, 'wfdisc': db_wfdisc}
    else:
        db_info = None
    stream, latlon = data_io.set_stream(local_wvfrms, fdsn, db_info, network, station, location, channel, starttime, endtime, local_latlon)

    click.echo('\n' + "Data summary:")
    for tr in stream:
        click.echo(tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime))

    if local_fk_label is None or local_fk_label == "auto":
        local_fk_label = ""
        if local_wvfrms is not None:
            if "/" in local_wvfrms:
                local_fk_label = os.path.dirname(local_wvfrms) + "/"
        
        local_fk_label = local_fk_label + tr.stats.network + "." + os.path.commonprefix([tr.stats.station for tr in stream])
        local_fk_label = local_fk_label + '_' + "%02d" % tr.stats.starttime.year + ".%02d" % tr.stats.starttime.month + ".%02d" % tr.stats.starttime.day
        local_fk_label = local_fk_label + '_' + "%02d" % tr.stats.starttime.hour + "." + "%02d" % tr.stats.starttime.minute + "." + "%02d" % tr.stats.starttime.second
        local_fk_label = local_fk_label + '-' + "%02d" % tr.stats.endtime.hour + "." + "%02d" % tr.stats.endtime.minute + "." + "%02d" % tr.stats.endtime.second
    else:
        if local_fk_label[-15:] == ".fk_results.dat":
            local_fk_label = local_fk_label[:-15]

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

    if back_az is not None and trace_vel is not None:
        click.echo('\n' + "Computing best beam with user specified beam...")
        click.echo('\t' + "Back Azimuth: " + str(back_az))
        click.echo('\t' + "Trace Velocity: " + str(trace_vel))

        stream.filter('bandpass', freqmin=freq_min, freqmax=freq_max)
        x, t, t0, geom = beamforming_new.stream_to_array_data(stream, latlon=latlon)
        X, _, f = beamforming_new.fft_array_data(x, t, fft_window="boxcar")

        sig_est, residual = beamforming_new.extract_signal(X, f, [back_az, trace_vel], geom)
        best_beam = np.fft.irfft(sig_est)[:len(t)] / (t[1] - t[0])
        residuals = np.fft.irfft(residual, axis=1)[:, :len(t)]  / (t[1] - t[0])

    else:
        click.echo('\n' + "Computing adaptive best beam...")
        click.echo('\t' + "fk results file: " + local_fk_label + ".fk_results.dat")

        def _envelope(t0, t1, t2, sigma):
            X1 = np.exp(-(t0 - t1) / sigma)
            X2 = np.exp(-(t0 - t2) / sigma)

            return X2 / ((1.0 + X1) * (1.0 + X2))

        # Read in the fk_results
        fk_t0 = None
        freq_min, freq_max = 0.5, 5.0

        temp = open(local_fk_label + ".fk_results.dat", 'r')
        for line in temp:
            if "t0:" in line:
                fk_t0 = np.datetime64(line.strip('\n').split(' ')[-1][:-1])
            elif "freq_min" in line:
                freq_min = float(line.split(' ')[-1])
            elif "freq_max" in line:
                freq_max = float(line.split(' ')[-1])
        temp.close()

        temp = np.loadtxt(local_fk_label + ".fk_results.dat")
        dt, beam_results = temp[:, 0], temp[:, 1:]
        beam_times = np.array([fk_t0 + np.timedelta64(int(dt_n * 1e3), 'ms') for dt_n in dt])

        # Filter and extract stream info
        stream.filter('bandpass', freqmin=freq_min, freqmax=freq_max)
        x, t, t0, geom = beamforming_new.stream_to_array_data(stream, latlon=latlon)
        M, _ = x.shape

        best_beam = np.zeros_like(t)
        residuals = np.zeros_like(x)

        window_step = (beam_times[1] - beam_times[0]).astype(float) / 1.0e6
        for n, tn in enumerate((beam_times - t0).astype(float) / 1.0e6):
            if tn >= 0.0 and tn <= (t[-1] - t[0]):
                X, _, f = beamforming_new.fft_array_data(x, t, window=[tn - window_step, tn + window_step], fft_window="boxcar")

                sig_est, residual = beamforming_new.extract_signal(X, f, beam_results[n, :2], geom)
                signal_wvfrm = np.fft.irfft(sig_est) / (t[1] - t[0])
                resid_wvfrms = np.fft.irfft(residual, axis=1) / (t[1] - t[0])

                mask = np.logical_and(tn - window_step <= t, t <= tn + window_step)
                best_beam[mask] = best_beam[mask] + signal_wvfrm[:sum(mask.astype(int))] * _envelope(t[mask], tn - window_step / 2.0, tn + window_step / 2.0, window_step / 20.0)
                for nM in range(M):
                    residuals[nM][mask] = residuals[nM][mask] + resid_wvfrms[nM][:sum(mask.astype(int))] * _envelope(t[mask], tn - window_step / 2.0, tn + window_step / 2.0, window_step / 20.0)

    # add output of waveform data (columns: t : beam : resid_1 : resid_2 : ... : resid_M)
    click.echo('\n' + "Writing results into " + local_fk_label + ".best-beam.dat" + '\n')
    header = "InfraPy Best Beam Results" + '\n'
    header = header + '\n' + "Data summary:" + '\n'
    for tr in stream:
        header = header + "    " + tr.id + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime) + '\n'

    header = header + "  t0: " + str(stream[0].stats.starttime) + '\n\n'

    if back_az is not None and trace_vel is not None:
        header = header + "Back Azimuth: " + str(back_az)
        header = header + "Trace Velocity: " + str(trace_vel) + '\n'
    else:
        header = header + "Beamforming (fk) results file: " + local_fk_label + ".fk_results.dat" + '\n'

    header = header + '\n' + "Column summary:" + '\n'
    header = header + "time (rel t0) [s] : beam [Pa] : Resid. 1 [Pa] : Resid. 2[Pa] : ... : Resid. M [Pa]" + '\n'

    output_vals = np.vstack((t, np.vstack((best_beam, residuals))))
    np.savetxt(local_fk_label + ".best-beam.dat", output_vals.T, header=header)
    
    # visualize results with UTC time
    plot_times = np.array([t0 + np.timedelta64(int(tn * 1000.0), 'ms') for tn in t])
    plt.plot(plot_times, best_beam, '-k', linewidth=1.0)
    for nM in range(len(residuals)):
        plt.plot(plot_times, residuals[nM], '-r', linewidth=0.25)

    if hold_figure:
        plt.show()
    else:
        plt.show(block=False)
        plt.pause(5.0)
        plt.close()
