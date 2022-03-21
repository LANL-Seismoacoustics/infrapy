#!/usr/bin/env python

import os
from threading import local 
import warnings 
import fnmatch
import json

import numpy as np

from obspy.clients.fdsn import Client
from obspy import read as obspy_read
from obspy import UTCDateTime, read_inventory

from ..propagation import likelihoods as lklhds

blank_sac_dict = {'delta': None, 'npts': None, 'depmin': None, 'depmax': None, 'depmen': None, 'b': 0.0, 'e': None, 'stla': None, 'stlo': None, 
                  'nzyear': None, 'nzjday': None, 'nzhour': None, 'nzmin': None, 'nzsec': None, 'nzmsec': None, 'kstnm': None, 'kcmpnm': None, 'knetwk': None}

############################
##     Data Ingestion     ##
##         Methods        ##
############################
def wvfrms_from_db(db_url, network, station, location, channel, starttime, endtime):
    # Need Jonathan or Christine to set this up...will be util.database eventually
    return None, None


def set_stream(local_opt, fdsn_opt, db_opt, network=None, station=None, location=None, channel=None, starttime=None, endtime=None, local_latlon=None):
    # check that only one option is selected and issue warning if multiple data sources are specified
    if np.sum(np.array([val is not None for val in [local_opt, fdsn_opt, db_opt]])) > 1:
        msg = '\n' + "Multiple data sources specified. Unexpected behavior is possible." + '\n' + "Priority order is [local > FDSN > DB]"
        warnings.warn(msg)

    # Check data option and populate obspy Stream
    if local_opt is not None:
        print('\n' + "Loading local data from " + local_opt)
        stream = obspy_read(local_opt)
        if local_latlon:
            latlon = np.load(local_latlon)
        else:
            latlon = None

    elif fdsn_opt is not None:
        print('\n' + "Loading data from FDSN (" + fdsn_opt + ")...")
        client = Client(fdsn_opt)
        stream = client.get_waveforms(network, station, location, channel, UTCDateTime(starttime), UTCDateTime(endtime), attach_response = True)
        stream.remove_response()  

        inventory = client.get_stations(network=network, station=station, starttime=UTCDateTime(starttime), endtime=UTCDateTime(endtime))
        latlon = []
        for network in inventory:
            for station in network:
                latlon = latlon + [[station.latitude, station.longitude]]

    elif db_opt is not None:
        print('\n' + "Loading data from database (methods not set up yet, so returning None's)...")
        stream, latlon = wvfrms_from_db(db_opt, network, station, location, channel, starttime, endtime)

    else:
        msg = "Warning: No waveform data source specified."
        warnings.warn(msg)
        stream, latlon = None, None

    return stream, latlon


def set_det_list(local_detect_label, merge=True):

    if ".dets.json" not in local_detect_label:
        local_detect_label = local_detect_label + ".dets.json"

    if "*" not in local_detect_label:
        print("Loading detections from file: " + local_detect_label)
        det_list = lklhds.json_to_detection_list(local_detect_label)
    else:
        if len(os.path.dirname(local_detect_label)) > 0:
            file_path = os.path.dirname(local_detect_label) + "/"
        else:
            file_path = ""

        file_list = []
        if "/" in local_detect_label:
            dir_files = os.listdir(os.path.dirname(local_detect_label))
        else:
            dir_files = os.listdir(".")
            
        for file in dir_files:
            if fnmatch.fnmatch(file, os.path.basename(local_detect_label)):
                file_list += [file]

        if len(file_list) == 0:
            msg = '\n' + "Detection file(s) specified not found"
            warnings.warn(msg)
            det_list = None 
        elif len(file_list) == 1:
            print("Loading detections from file: " + file_path + local_detect_label)
            det_list = [lklhds.json_to_detection_list(file_path + local_detect_label)]
        else:
            print("Loading detections from files:")
            det_list = []
            for file in file_list:
                print('\t' + file_path + file)

                if merge:
                    det_list = det_list + lklhds.json_to_detection_list(file_path + file)
                else:
                    det_list = det_list + [lklhds.json_to_detection_list(file_path + file)]

    return det_list


##########################
##     Data Writing     ##
##        Methods       ##
##########################
def write_stream(stream, latlon):
    sac_info = [blank_sac_dict] * len(stream)
    for m, tr in enumerate(stream):
        sac_info[m]['delta'] = tr.stats.delta
        sac_info[m]['npts'] = tr.stats.npts
        sac_info[m]['e'] = tr.stats.npts * tr.stats.delta

        sac_info[m]['depmin'] = min(tr.data)
        sac_info[m]['depmax'] = max(tr.data)
        sac_info[m]['depmen'] = np.mean(tr.data)

        sac_info[m]['stla'] = latlon[m][0]
        sac_info[m]['stlo'] = latlon[m][1]

        sac_info[m]['nzyear'] = tr.stats.starttime.year
        sac_info[m]['nzjday'] = tr.stats.starttime.julday
        sac_info[m]['nzhour'] = tr.stats.starttime.hour
        sac_info[m]['nzmin'] = tr.stats.starttime.minute
        sac_info[m]['nzsec'] = tr.stats.starttime.second

        sac_info[m]['knetwk'] = tr.stats.network
        sac_info[m]['kstnm'] = tr.stats.station
        sac_info[m]['kcmpnm'] = tr.stats.channel
        
        tr.stats.sac = sac_info[m]

        label = tr.stats.network + "." + tr.stats.station
        label = label + '_' + "%02d" % tr.stats.starttime.year + ".%02d" % tr.stats.starttime.month + ".%02d" % tr.stats.starttime.day
        label = label + '_' + "%02d" % tr.stats.starttime.hour + "." + "%02d" % tr.stats.starttime.minute + "." + "%02d" % tr.stats.starttime.second

        tr.write(label + ".sac", format='SAC') 


def fk_header(stream, latlon, freq_min, freq_max, back_az_min, back_az_max, back_az_step, trace_vel_min, trace_vel_max, trace_vel_step, method, signal_start, signal_end, noise_start, noise_end, window_len, sub_window_len, window_step):
    header = "InfraPy Beamforming (fk) Results" + '\n'
    header = header + '\n' + "Data summary:" + '\n'
    for tr in stream:
        header = header + "    " + tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime) + '\n'

    header = header + '\n' + "  channel_cnt: " + str(len(stream)) + '\n'

    if latlon:
        mean_lat = latlon[0][0]
        mean_lon = latlon[0][1]
    else:
        mean_lat = stream[0].stats.sac['stla']
        mean_lon = stream[0].stats.sac['stlo']

    header = header + "  t0: " + str(stream[0].stats.starttime) + '\n'
    header = header + "  latitude: " + str(mean_lat) + '\n'
    header = header + "  longitude: " + str(mean_lon) + '\n'

    header = header + '\n' + "Algorithm parameters:" + '\n'
    header = header + "  freq_min: " + str(freq_min) + '\n'
    header = header + "  freq_max: " + str(freq_max) + '\n'
    header = header + "  back_az_min: " + str(back_az_min) + '\n'
    header = header + "  back_az_max: " + str(back_az_max) + '\n'
    header = header + "  back_az_step: " + str(back_az_step) + '\n'
    header = header + "  trace_vel_min: " + str(trace_vel_min) + '\n'
    header = header + "  trace_vel_max: " + str(trace_vel_max) + '\n'
    header = header + "  trace_vel_step: " + str(trace_vel_step) + '\n'
    header = header + "  method: " + str(method) + '\n'
    header = header + "  signal_start: " + str(signal_start) + '\n'
    header = header + "  signal_end: " + str(signal_end) + '\n'
    if method == "GLS":
        header = header + "  noise_start: " + str(noise_start) + '\n'
        header = header + "  noise_end: " + str(noise_end) + '\n'
    header = header + "  window_len: " + str(window_len) + '\n'
    header = header + "  sub_window_len: " + str(sub_window_len) + '\n'
    header = header + "  window_step: " + str(window_step) + '\n'

    header = header + '\n' + "Time (rel t0) [s]      Back Az [deg]	           Tr. Velocity [m/s]       F-stat"

    return header 


def define_detection(det_info, array_loc, channel_cnt, freq_band, note=None):
    temp = lklhds.InfrasoundDetection()
    temp.latitude = float(array_loc[0])
    temp.longitude = float(array_loc[1])
    temp.array_dim = int(channel_cnt)
    temp.frequency_range = freq_band

    temp.peakF_UTCtime = det_info[0]
    temp.start = det_info[1]
    temp.end = det_info[2]
    temp.back_azimuth = np.round(det_info[3], 2)
    temp.trace_velocity = np.round(det_info[4], 2)
    temp.peakF_value = np.round(det_info[5], 4)
    temp.note = note

    return temp


def write_events(events, event_qls, det_list, local_event_label):
    for ev_n, ev in enumerate(events):
        temp = []
        for det_id in ev:
            temp = temp + [det_list[det_id]]
        lklhds.detection_list_to_json(local_event_label + "-ev" + str(ev_n) + ".dets.json", temp)


def write_locs(bisl_results, local_loc_label):
    if ".loc.json" in local_loc_label:
        with open(local_loc_label, 'w') as of:
            json.dump(bisl_results, of, indent=4, cls=lklhds.Infrapy_Encoder)
    else:
        with open(local_loc_label + ".loc.json", 'w') as of:
            json.dump(bisl_results, of, indent=4, cls=lklhds.Infrapy_Encoder)


def read_locs(local_loc_label):
    if ".loc.json" in local_loc_label:
        return json.load(open(local_loc_label))
    else:
        return json.load(open(local_loc_label + ".loc.json"))
