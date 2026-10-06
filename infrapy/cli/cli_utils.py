#!which python
"""
cli_utils.py

Utility methods accessible in the command line interface (CLI) of infrapy

Author: pblom@lanl.gov
"""

import os
import click
import json
import gzip

import configparser as cnfg

import numpy as np

from obspy import UTCDateTime

from infrapy.characterization import spye
from infrapy.detection import beam
from infrapy.propagation import likelihoods as lklhds
from infrapy.utils import config, data_io, database

# Rewrite this as "db2wvfrms" with option to summarize or write to SAC
@click.command('check_db_wvfrms', short_help="Check waveform pull from database")
@click.option("--cnfg-file", help="Configuration file", default=None)
@click.option("--db-config", help="Database configuration file", default=None)

@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)

@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
def check_db_wvfrm(cnfg_file, db_config, network, station, location, channel, starttime, endtime):
    '''
    Test database pull of waveform data for beamforming (fk or fdk) analysis

    \b
    Example usage (detection_db.config will be unique to your database pull):
    \tinfrapy run_fk --cnfg-file config/detection_db.config

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##       check_db_wvfrms       ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    if cnfg_file:
        click.echo('\n' + "Loading configuration info from: " + cnfg_file)
        if os.path.isfile(cnfg_file):
            user_config = cnfg.ConfigParser()
            user_config.read(cnfg_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Database and data IO parameters
    db_config = config.set_param(user_config, 'DATA IO', 'db_config', db_config, 'string')
    db_info = None

    network = config.set_param(user_config, 'DATA IO', 'network', network, 'string')
    station = config.set_param(user_config, 'DATA IO', 'station', station, 'string')
    location = config.set_param(user_config, 'DATA IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'DATA IO', 'channel', channel, 'string')

    starttime = config.set_param(user_config, 'DATA IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'DATA IO', 'endtime', endtime, 'string')

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


@click.command('db2dets', short_help="Write from arrivals table to dets.json.gz file")
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--lat-bnds", help="Latitude bounds", default=None, prompt="Latitude bounds (comma separated):")
@click.option("--lon-bnds", help="Longitude bounds", default=None, prompt="Longitude bounds (comma separated):")
@click.option("--starttime", help="Start time of analysis window", default=None, prompt="Window start time:")
@click.option("--endtime", help="End time of analysis window", default=None, prompt="Window end time:")
@click.option("--phase-list", help="Phases to include in output (default: 'I')", default="I")
@click.option("--output-label", help="Output label for [...].ev.json.gz file", default=None)
@click.option("--verbose", help="Print retrieved event info to screen (default: True)", default=True)
def db2dets(db_config, lat_bnds, lon_bnds, starttime, endtime, phase_list, output_label, verbose):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##           db2dets           ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    lat_lims = [float(val) for val in lat_bnds.split(",")]
    lon_lims = [float(val) for val in lon_bnds.split(",")]
    starttime = UTCDateTime(starttime)
    endtime = UTCDateTime(endtime)

    print("Pulling arrivals for criterion:")
    print("  Lat/Lon bounds: [" + str(lat_lims[0]) + ", " + str(lon_lims[0]) + "] - [" + str(lat_lims[1]) + ", " + str(lon_lims[1]) + "]")
    print("  Time bounds:", starttime, ",", endtime)
    print("  Phase list:", phase_list)

    db_info = cnfg.ConfigParser()
    db_info.read(db_config)

    print("Setting up database configuration...")
    # set up the session and check connection
    if 'url' in db_info['DATABASE'].keys():
        print("  Connecting to database through url: " + db_info['DATABASE']['url'])
        db_session = database.db_connect_url( db_info['DATABASE']['url'])
    else:
        # clean up the above to simplify this or just require a url?
        db_session = database.db_connect2(db_info)

    # check the session works
    try:
        db_session.get_bind().connect()
    except Exception as e:
        print("Database connection failed")
        return

    det_dicts = database.db2dets(db_session, db_info['DBTABLES'], lat_lims, lon_lims, starttime, endtime, phase_list=phase_list, db_schema="kbcore")

    if len(det_dicts) > 0:
        if verbose:
            print('\n' + str(len(det_dicts)) + " phase(s) found for criterion...")
            for det in det_dicts:
                print("  " + det['wvfrm_info'][0]['trace id'] + ' ' * (16 - len(det['wvfrm_info'][0]['trace id'])), end='\t')
                print(det['phase id'], end='\t')
                print(det['peak f-stat time'], end='\t')
                print(np.round(np.array(det['back az'], dtype=float), 1), end='\t')
                print(np.round(np.array(det['tr vel'], dtype=float), 3), end='\t')
                print(np.round(np.array(det['f-stat'], dtype=float), 1))

        if output_label is not None:
            print("Writing " + str(len(det_dicts)) + " arrival entries into " + output_label + ".dets.json.gz")
            with gzip.open(output_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
                json.dump({"det_info" : det_dicts}, zipfile, indent=4, cls=data_io.Infrapy_Encoder)
    else:
        print('\n' + "No arrivals matched criterion.")


@click.command('db2ev', short_help="Write from arrivals table to ev.json.gz file")
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--evid", help="Event ID to pull", default=0)
@click.option("--phase-list", help="Phases to include in output (default: 'I')", default="I")
@click.option("--output-label", help="Output label for [...].ev.json.gz file", default=None)
@click.option("--verbose", help="Print retrieved event info to screen (default: True)", default=True)
def db2ev(db_config, evid, phase_list, output_label, verbose):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##            db2ev            ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    db_info = cnfg.ConfigParser()
    db_info.read(db_config)

    print("Setting up database configuration...")
    # set up the session and check connection
    if 'url' in db_info['DATABASE'].keys():
        print("  Connecting to database through url: " + db_info['DATABASE']['url'])
        db_session = database.db_connect_url( db_info['DATABASE']['url'])
    else:
        # clean up the above to simplify this or just require a url?
        db_session = database.db_connect2(db_info)

    # check the session works
    try:
        db_session.get_bind().connect()
    except:
        print("Database connection failed")
        return

    ev_output = database.db2ev(db_session, db_info['DBTABLES'], evid, phase_list=phase_list, db_schema="kbcore")

    if len(ev_output['det_info']) > 0:
        if verbose:
            print('\n\n' + str(len(ev_output['det_info'])) + " included phase(s) found for evid: " + str(evid))
            print('Preferred origin info:\n  Location: ' + str(ev_output['ground truth']['latitude']) + ', ' + str(ev_output['ground truth']['longitude']))
            print('  Origin time: ' + str(ev_output['ground truth']['origin time']))
            print('  Name: ' + str(ev_output['ground truth']['name'] + '\n\nDetections list:'))

            for det in ev_output['det_info']:
                print("  " + det['wvfrm_info'][0]['trace id'] + ' ' * (16 - len(det['wvfrm_info'][0]['trace id'])), end='\t')
                print(det['phase id'], end='\t')
                if "I" in det['phase id']:
                    print(det['peak f-stat time'], end='\t')
                    print(np.round(np.array(det['back az'], dtype=float), 1), end='\t')
                    print(np.round(np.array(det['tr vel'], dtype=float), 3), end='\t')
                    print(np.round(np.array(det['f-stat'], dtype=float), 1))
                else:
                    print(det['peak f-stat time'], end='\t')
                    print(np.round(np.array(det['azimuth'], dtype=float), 1), end='\t')
                    print(np.round(np.array(det['slow'], dtype=float), 3), end='\t')
                    print(np.round(np.array(det['snr'], dtype=float), 1))

        if output_label is not None:
            print("Writing event info including " + str(len(ev_output['det_info'])) + " arrival entries into " + output_label + ".ev.json.gz")
            with gzip.open(output_label + ".ev.json.gz", 'wt', encoding='UTF-8') as zipfile:
                json.dump(ev_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)

# Rewrite this as "db2wvfrms" with option to summarize or write to SAC

@click.command('write_wvfrms', short_help="Save waveforms from FDSN or database")
@click.option("--cnfg-file", help="Configuration file", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--network", help="Network code for FDSN and database", default=None)
@click.option("--station", help="Station code for FDSN and database", default=None)
@click.option("--location", help="Location code for FDSN and database", default=None)
@click.option("--channel", help="Channel code for FDSN and database", default=None)
@click.option("--starttime", help="Start time of analysis window", default=None)
@click.option("--endtime", help="End time of analysis window", default=None)
def write_wvfrms(cnfg_file, db_config, fdsn, network, station, location, channel, starttime, endtime):
    '''
    Write waveform data from an FDSN or database pull into local SAC files

    \b
    Example usage (detection_db.config will be unique to your database pull):
    \tinfrapy utils write-wvfrms --cnfg-file config/detection_fdsn.config

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         write-wvfrms        ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    if cnfg_file:
        click.echo('\n' + "Loading configuration info from: " + cnfg_file)
        if os.path.isfile(cnfg_file):
            user_config = cnfg.ConfigParser()
            user_config.read(cnfg_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Database and data IO parameters
    db_config = config.set_param(user_config, 'DATA IO', 'db_config', db_config, 'string')
    db_info = None

    # FDSN DATA IO parameters
    fdsn = config.set_param(user_config, 'DATA IO', 'fdsn', fdsn, 'string')
    network = config.set_param(user_config, 'DATA IO', 'network', network, 'string')
    station = config.set_param(user_config, 'DATA IO', 'station', station, 'string')
    location = config.set_param(user_config, 'DATA IO', 'location', location, 'string')
    channel = config.set_param(user_config, 'DATA IO', 'channel', channel, 'string')

    # Trimming times
    starttime = config.set_param(user_config, 'DATA IO', 'starttime', starttime, 'string')
    endtime = config.set_param(user_config, 'DATA IO', 'endtime', endtime, 'string')

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


@click.command('merge_dets', short_help="Check waveform pull from database")
@click.option("--dets-files", help="Detection GZIP files", default=None)
@click.option("--merged-label", help="Output detection file label", default=None)
def merge_dets(dets_files, merged_label):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##         merge_dets          ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    dets_data = data_io._load_dets_json(dets_files)

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

    # Copy waveform and parameter info into individual detections
    det_list = []
    for entry in dets_data:
        for det in entry["det_info"]:
            det_list = det_list + [det]
            det_list[-1]["wvfrm_info"] = entry["wvfrm_info"]
            det_list[-1]["fk_params"] = entry["fk_params"]
            det_list[-1]["det_params"] = entry["det_params"]

    click.echo('\n' + "Merging detections...")
    det_output = {'det_info' : beam.merge_dets(det_list)}
    with gzip.open(merged_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)


@click.command('ev_gt', short_help="Populate event ground truth information")
@click.option("--ev-file", help="Event GZIP JSON files", default=None)

@click.option("--latitude", help="event latitude (deg)", default=None)
@click.option("--longitude", help="Event longitude (deg)", default=None)
@click.option("--orig-tm", help="Origin datetime", default=None)
@click.option("--eq-tnt", help="Explosive yield (eq. TNT) [kg]", default=None)

@click.option("--user-entry", help="User specified info (comma sep. key/val)", default=None, multiple=True)
@click.option("--entry-mode", help="Add values or overwite ('append' or 'replace')", default='append')

def ev_gt(ev_file, latitude, longitude, orig_tm, eq_tnt, user_entry, entry_mode):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##     Write Event GT Info     ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    click.echo("Loading event information from ev_file: " + str(ev_file))
    ev_data = data_io._load_dets_json(ev_file)[0]
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
    with gzip.open(ev_file, 'wt', encoding='UTF-8') as zipfile:
        json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)


@click.command('ev_summary', short_help="Summarize information in an event file")
@click.option("--ev-file", help="Event GZIP JSON files", default=None)
def ev_summary(ev_file):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##     Summarize Event File    ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    click.echo("Loading information from ev_file: " + str(ev_file))
    ev_data = data_io._load_dets_json(ev_file)[0]

    click.echo('\n' + "=" * 17 + '\n' + "Detection Summary" + '\n' + "=" * 17)
    try:
        det_specs = [spye.extract_json_spectra(det['spec']) for det in ev_data['det_info']]
    except:
        det_specs = None

    for k, det in enumerate(ev_data['det_info']):

        click.echo(det['wvfrm_info'][0][0]['trace id'] + '\t', nl=False)
        click.echo("  loc: " + f"{det['wvfrm_info'][0][0]['latitude']:>7.4f}" + ", " + f"{det['wvfrm_info'][0][0]['longitude']:>8.4f}" + '\t', nl=False)
        click.echo("  time: " + det['peak f-stat time'].split(".")[0] + '\t', nl=False)
        click.echo("  back az [deg]: " + f"{det["back az"]:>8.2f}" + '\t', nl=False)
        click.echo("  tr vel [m/s]: " + f"{det["tr vel"]:.1f}" + '\t', nl=False)
        click.echo("  f-stat: " + f"{det["f-stat"]:>5.1f}" + '\t', nl=False)
        if det_specs is not None:
            freq = np.array(det_specs[k][0])
            mask = np.array((np.array(det_specs[k][1]) / np.array(det_specs[k][2])) > 2.0)
            click.echo("  high snr band: " + f"{np.round(freq[mask][0], 2):.2f}" + " - " + f"{np.round(freq[mask][-1], 2):.2f}" + ' Hz')
        else:
            click.echo('')

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


@click.command('ev_loc_reset', short_help="Reset localization in an event file")
@click.option("--ev-file", help="Event GZIP JSON files", default=None)
def ev_loc_reset(ev_file):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##       Reset Event File      ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    click.echo("Loading information from ev_file: " + str(ev_file))
    ev_data = data_io._load_dets_json(ev_file)[0]

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

        click.echo('\n' + '#' * 40 + '\n' + '#' * 40 + '\n')

    user_opt = input('WARNING!!! This action will remove existing localization result(s) in this event file. \nDo you want to proceed? (y/n): ').lower().strip()
    if user_opt in ['y', 'yes']:
        confirm = True
    else:
        confirm = False

    if confirm:
        click.echo('\nRemoving localization results from event file...')
        ev_data['location'] = []
        with gzip.open(ev_file, 'wt', encoding='UTF-8') as zipfile:
            json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)


@click.command('ev_reset', short_help="Reset location and/or characterization in an event file")
@click.option("--ev-file", help="Event GZIP JSON files", default=None)
@click.option("--reset-loc", help="Flag to reset location info", default=False)
@click.option("--reset-char", help="Flag to reset characterization info", default=False)
def ev_reset(ev_file, reset_loc, reset_char):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##       Reset Event File      ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    click.echo("Loading information from ev_file: " + str(ev_file))
    ev_data = data_io._load_dets_json(ev_file)[0]
    click.echo("  " + str(len(ev_data['location'])) + " location result(s) in file")
    click.echo("  " + str(len(ev_data['characterization'])) + " characterization result(s) in file" + '\n')

    if not reset_loc and not reset_char:
        click.echo("Set '--reset-loc' and/or '--reset-char' to 'True' to remove existing solutions from file." + '\n')
        return

    if reset_loc:
        if len(ev_data['location']) == 0:
            click.echo("No location results to remove.")
            reset_loc = False
        else:
            click.echo('\n' + "=" * 20 + '\n' + "Localization Summary" + '\n' + "=" * 20)
            for loc_k, loc in enumerate(ev_data['location']):
                click.echo("#" * 14)
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

    if reset_char:
        if len(ev_data['characterization']) == 0:
            click.echo("No characterization results to remove.")
            reset_char = False
        else:
            click.echo('\n' + "=" * 24 + '\n' + "Characterization Summary" + '\n' + "=" * 24)
            for char_k, char in enumerate(ev_data['characterization']):
                click.echo("#" * 14)
                click.echo("## " + "index: " + str(char_k) + " ##")
                click.echo("#" * 14)

                click.echo("parameters" + '\n' + "-" * 10)
                for key in char['params'].keys():
                    if char['params'][key] is not None:
                        click.echo("    " + key + ": " + str(char['params'][key]))

                click.echo('\n' + "result" + '\n' + "-" * 6)
                click.echo("    maximum likelihood yield: " + str(np.round(char['result']['yld_vals'][np.argmax(char['result']['yld_pdf'])], 2)) + " tons eq. TNT")
                click.echo("    68% confidence bounds: " + str(char['result']['conf_bnds'][0]))
                click.echo("    95% confidence bounds: " + str(char['result']['conf_bnds'][1]) + '\n')

    if reset_loc or reset_char:
        click.echo('\n' + '#' * 40 + '\n' + '#' * 40 + '\n')
        while True:
            user_opt = input('WARNING!!! This action will remove existing result(s) in this event file. \nDo you want to proceed? (y/n): ').lower().strip()
            if user_opt in ['y', 'yes']:
                confirm = True
                break
            elif user_opt in ['n', 'no']:
                click.echo("")
                confirm = False
            else:
                user_opt = input('Invalid input. Proceed and remove result(s) from event file? (y/n): ').lower().strip()

        if confirm:
            if reset_loc:
                click.echo('\nRemoving localization result(s) from event file...')
                ev_data['location'] = []
            if reset_char:
                click.echo('\nRemoving characterization result(s) from event file...')
                ev_data['characterization'] = []

            with gzip.open(ev_file, 'wt', encoding='UTF-8') as zipfile:
                json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)


@click.command('ev_char_reset', short_help="Reset characterization in an event file")
@click.option("--ev-file", help="Event GZIP JSON files", default=None)
def ev_char_reset(ev_file):

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##      InfraPy Utilities      ##")
    click.echo("##       Reset Event File      ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")

    click.echo("Loading information from ev_file: " + str(ev_file))
    ev_data = data_io._load_dets_json(ev_file)[0]

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

        click.echo('\n' + '#' * 40 + '\n' + '#' * 40 + '\n')

    user_opt = input('WARNING!!! This action will remove existing characterization result(s) in this event file. \nDo you want to proceed? (y/n): ').lower().strip()
    if user_opt in ['y', 'yes']:
        confirm = True
    else:
        confirm = False

    if confirm:
        click.echo('\nRemoving characterization results from event file...')
        ev_data['characterization'] = []
        with gzip.open(ev_file, 'wt', encoding='UTF-8') as zipfile:
            json.dump(ev_data, zipfile, indent=4, cls=data_io.Infrapy_Encoder)




