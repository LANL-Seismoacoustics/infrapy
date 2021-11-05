#!/usr/bin/env python

import os 
import warnings 
import fnmatch

import numpy as np

from obspy.clients.fdsn import Client
from obspy.core import read as obspy_read
from obspy import UTCDateTime

from ..propagation import likelihoods as lklhds


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
        msg = "Multiple data sources specified. Unexpected behavior is possible." + '\n' + "Priority order is [local > FDSN > DB]"
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
        stream = client.get_waveforms(network, station, location, channel, UTCDateTime(starttime), UTCDateTime(endtime))
        inventory = client.get_stations(network=network, station=station, starttime=UTCDateTime(starttime), endtime=UTCDateTime(endtime))

        latlon = []
        for network in inventory:
            for station in network:
                latlon = latlon + [[station.latitude, station.longitude]]

        return stream, latlon

    elif db_opt is not None:
        print('\n' + "Loading data from database (methods not set up yet, so returning 'None's)...")
        stream, latlon = wvfrms_from_db(db_opt, network, station, location, channel, starttime, endtime)

    else:
        print("No data source specified.")
        stream, latlon = None, None

    return stream, latlon

def set_det_list(local_det_info, merge=True):

    if "*" not in local_det_info:
        det_list = lklhds.json_to_detection_list(local_det_info)
    else:
        if len(os.path.dirname(local_det_info)) > 0:
            file_path = os.path.dirname(local_det_info) + "/"
        else:
            file_path = ""

        file_list = []
        dir_files = os.listdir(os.path.dirname(local_det_info))
        for file in dir_files:
            if fnmatch.fnmatch(file, os.path.basename(local_det_info)):
                file_list += [file]

        if len(file_list) == 0:
            msg = "Detection file(s) specified not found"
            warnings.warn(msg)
            det_list = None 
        elif len(file_list) == 1:
            det_list = [lklhds.json_to_detection_list(file_path + local_det_info)]
        else:
            det_list = []
            for file in file_list:
                if merge:
                    det_list = det_list + lklhds.json_to_detection_list(file_path + file)
                else:
                    det_list = det_list + [lklhds.json_to_detection_list(file_path + file)]

    return det_list

##########################
##     Data Writing     ##
##        Methods       ##
##########################
def write_fk_meta(stream, latlon, local_fk_out, freq_min, freq_max, back_az_min, back_az_max, back_az_step, trace_vel_min, trace_vel_max, trace_vel_step, method, 
    signal_start, signal_end, noise_start, noise_end, window_len, sub_window_len, window_step):
        file_out = open(local_fk_out + ".fk_meta.txt", 'w')

        print("InfraPy Beamforming (fk) Analysis", file=file_out)
        print("---------------------------------", file=file_out)

        print('\n' + "Data summary:", file=file_out)
        for tr in stream:
            print("  " + tr.stats.network + "." + tr.stats.station + "." + tr.stats.location + "." + tr.stats.channel + '\t' + str(tr.stats.starttime) + " - " + str(tr.stats.endtime), file=file_out)

        print('\n' + "  channel_cnt: " + str(len(stream)), file=file_out)
        if latlon:
            mean_lat = latlon[0][0]
            mean_lon = latlon[0][1]
        else:
            mean_lat = stream[0].stats.sac['stla']
            mean_lon = stream[0].stats.sac['stlo']

        print("  latitude: " + str(mean_lat), file=file_out)        
        print("  longitude: " + str(mean_lon), file=file_out)        

        print('\n' + "Algorithm parameters:", file=file_out)
        print("  freq_min: " + str(freq_min), file=file_out)
        print("  freq_max: " + str(freq_max), file=file_out)
        print("  back_az_min: " + str(back_az_min), file=file_out)
        print("  back_az_max: " + str(back_az_max), file=file_out)
        print("  back_az_step: " + str(back_az_step), file=file_out)
        print("  trace_vel_min: " + str(trace_vel_min), file=file_out)
        print("  trace_vel_max: " + str(trace_vel_max), file=file_out)
        print("  trace_vel_step: " + str(trace_vel_step), file=file_out)
        print("  method: " + str(method), file=file_out)
        print("  signal_start: " + str(signal_start), file=file_out)
        print("  signal_end: " + str(signal_end), file=file_out)
        if method == "GLS":
            print("  noise_start: " + str(noise_start), file=file_out)
            print("  noise_end: " + str(noise_end), file=file_out)
        print("  window_len: " + str(window_len), file=file_out)
        print("  sub_window_len: " + str(sub_window_len), file=file_out)
        print("  window_step: " + str(window_step), file=file_out)
        file_out.close()


def _define_deteection(det_info, array_loc, channel_cnt, freq_band, note=None):
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

def write_events(events, event_qls, det_list, local_events_out):
    print("Writing identified events into " + local_events_out)
    for ev_n, ev in enumerate(events):
        temp = []
        for det_id in ev:
            temp = temp + [det_list[det_id]]
        lklhds.detection_list_to_json(local_events_out + "-ev" + str(ev_n) + ".json", temp)




