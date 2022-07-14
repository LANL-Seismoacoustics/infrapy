#!/usr/bin/env python

import os
from threading import local 
import warnings 
import fnmatch
import json
import csv

import numpy as np

from obspy.clients.fdsn import Client
from obspy import read as obspy_read
from obspy import UTCDateTime, read_inventory

from ..propagation import likelihoods as lklhds
from . import database

blank_sac_dict = {'delta': None, 'npts': None, 'depmin': None, 'depmax': None, 'depmen': None, 'b': 0.0, 'e': None, 'stla': None, 'stlo': None, 
                  'nzyear': None, 'nzjday': None, 'nzhour': None, 'nzmin': None, 'nzsec': None, 'nzmsec': None, 'kstnm': None, 'kcmpnm': None, 'knetwk': None}

############################
##     Data Ingestion     ##
##         Methods        ##
############################

def wvfrms_from_fdsn(fdsn_opt, network, station, location, channel, starttime, endtime):
    """
    connect to an FDSN server to pull data

    Parameters
    ----------
    fsdn_opt: str
        FDSN option (e.g., IRIS); None if using another source
    network: str
        Network for FDSN and database options
    station: str
        Station for the FDSN and database options
    location: str
        Location for the FDSN and database options
    channel: str
        Channel for the FDSN and database options
    starttime: str
        Start time for the FDSN and database options; formatted to be compatible with obspy.UTCDateTime
    endtime: str
        End time for the FDSN and database options; formatted to be compatible with obspy.UTCDateTime

    Returns
    -------
    stream : obspy.core.stream.Stream
        Obspy stream containing specified waveform data
    latlon: 2darray
        Iterable with latitude and longitude info for each trace of the returned stream

    """

    client = Client(fdsn_opt)
    stream = client.get_waveforms(network, station, location, channel, UTCDateTime(starttime), UTCDateTime(endtime), attach_response = True)
    stream.remove_response()

    inventory = client.get_stations(network=network, station=station, starttime=UTCDateTime(starttime), endtime=UTCDateTime(endtime))
    latlon = []
    for network in inventory:
        for station in network:
            latlon = latlon + [[station.latitude, station.longitude]]

    return stream, latlon



def set_stream(local_opt, fdsn_opt, db_info, network=None, station=None, location=None, channel=None, starttime=None, endtime=None, local_latlon=None):
    """
    Define an ObsPy stream from a specified local, FDSN, or database source.
    1) if specifying local data, use obspy.read to set up the stream
    2) if pulling from an FDSN, use obspy.clients.fdsn.Client to pull waveforms and station info
    3) if pulling from a database...this is still in development

    Parameters
    ----------
    local_opt: str
        Local waveform files (must be readable by obspy.read); None is using another source
    fsdn_opt: str
        FDSN option (e.g., IRIS); None if using another source
    db_info: str
        Database info to pull data; None if using another source
    network: str
        Network for FDSN and database options
    station: str
        Station for the FDSN and database options
    location: str
        Location for the FDSN and database options
    channel: str
        Channel for the FDSN and database options
    starttime: str
        Start time for the FDSN and database options; formatted to be compatible with obspy.UTCDateTime
    endtime: str
        End time for the FDSN and database options; formatted to be compatible with obspy.UTCDateTime
    local_latlon: str
        File containing latlon info for local waveform data (need to add an option/method for a site file)


    Returns
    -------
    stream : obspy.core.stream.Stream
        Obspy stream containing specified waveform data
    latlon: 2darray
        Iterable with latitude and longitude info for each trace of the returned stream

    """

    # check that only one option is selected and issue warning if multiple data sources are specified
    if np.sum(np.array([val is not None for val in [local_opt, fdsn_opt, db_info]])) > 1:
        msg = '\n' + "Multiple data sources specified. Unexpected behavior is possible." + '\n' + "Priority order is [local > FDSN > DB]"
        warnings.warn(msg)

    # Check data option and populate obspy Stream
    if local_opt is not None:
        print('\n' + "Loading local data from " + local_opt)
        stream = obspy_read(local_opt)
        if local_latlon:
            latlon = np.load(local_latlon)
        else:
            latlon = [[tr.stats.sac['stla'], tr.stats.sac['stlo']] for tr in stream]

    elif fdsn_opt is not None:
        print('\n' + "Loading data from FDSN (" + fdsn_opt + ")...")
        stream, latlon = wvfrms_from_fdsn(fdsn_opt, network, station, location, channel, starttime, endtime)

    elif db_info is not None:
        print('\n' + "Loading data from database (" + db_info['url'] + ")...")
        stream, latlon = database.wvfrms_from_db(db_info, station, channel, starttime, endtime)

    else:
        msg = "Warning: No waveform data source specified."
        warnings.warn(msg)
        stream, latlon = None, None

    return stream, latlon


def set_det_list(local_detect_label, merge=True):
    """
    Read detections from a file (or files) using the [...].dets.json format used to output detections

    Parameters
    ----------
    local_detect_label: str
        String denoting detection file(s) to be loaded for analysis
    merge: bool
        Control for merging files into a single list (for event ID) or creating nested lists (for multiple localization analyses)

    Returns
    -------
    det_list : list
        List containing infrapy.propagation.likelihoods.InfrasoundDetection instances for analysis; if merge=False, returns list of lists of detections

    """

    if ".dets.json" not in local_detect_label:
        local_detect_label = local_detect_label + ".dets.json"

    if "*" not in local_detect_label:
        print("Loading detections from file: " + local_detect_label)
        det_list = json_to_detection_list(local_detect_label)
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
            det_list = [json_to_detection_list(file_path + local_detect_label)]
        else:
            print("Loading detections from files:")
            det_list = []
            for file in file_list:
                print('\t' + file_path + file)

                if merge:
                    det_list = det_list + json_to_detection_list(file_path + file)
                else:
                    det_list = det_list + [json_to_detection_list(file_path + file)]

    return det_list


##########################
##     Data Writing     ##
##        Methods       ##
##########################
def write_stream_to_sac(stream, latlon):
    """
    Write info from an obspy.core.stream.Stream instance into local sac files with populated header info.  Defines the output label from the network, station, and start/end times of the stream

    Parameters
    ----------
    stream: obspy.core.stream.Stream
        Stream of waveform data to be output
    latlon: 2darray
        Iterable containing latitude and longitude info for each trace of the stream
    """

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
    """
    Write fk (beamforming) analysis parameter info into a header for output of results

    Parameters
    ----------
    stream: obspy.core.stream.Stream
        Stream of waveform data used in analysis
    latlon: 2darray
        Iterable containing latitude and longitude info for each trace of the stream
    freq_min: float
        Minimum frequency used in fk analysis
    freq_max: float
        Maximum frequency used in fk analysis
    back_az_min: float
        Minimum back azimuth used in fk analysis
    back_az_max: float
        Maximum back azimuth used in fk analysis
    trace_vel_min: float
        Minimum trace velocity used in fk analysis
    trace_vel_max: float
        Maximum_trace_velocity used in fk analysis
    method: str
        Method (e.g., Bartlett, Capon) used in fk analysis
    signal_start: float
        Signal start [s] if not analyzing the entire stream
    signal_end: float
        Signal end [s] if not analyzing the entire stream
    noise_start: float
        Noise start [s] if using adaptive beamforming (GLS)
    noise_end: float
        Noise end [s] if using adaptive beamforming (GLS)
    window_len: float
        Analysis window length [s]
    sub_window_len: float
        Sub-window length if using a covariance matrix method (e.g., Bartlett_Covar, MUSIC)
    window_step: float
        Step between analysis windows [s] adjustable for overlapping windows


    Returns
    -------
    header : str
        Header for numpy.savetxt output of fk results

    """
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
    """
    Write detection info from fd analysis into a infrapy.propagation.likelihoods.InfrasoundDetection instance for output into a [...].dets.json file

    # I expanded the InfrasoundDetection constructor to include everything here, so maybe this is now redundant?


    Parameters
    ----------
    det_info: ndarray
        Detection info containing peak F-stat time, relative start/end, and direction of arrival info
    array_loc: iterable
        Latitude and longitude of the detecting array
    channel_cnt: int
        Number of channels in the detecting array
    freq_band: iterable
        Minimum and maximum frequencies used in analysis
    note: str
        Any note about the detection (e.g., 'InfraPy CLI Detection')

    Returns
    -------
    detection : infrapy.propagation.likelihoods.InfrasoundDetection
        Detection info in expected format

    """

    return lklhds.InfrasoundDetection(lat_loc=float(array_loc[0]), 
                                      lon_loc=float(array_loc[1]), 
                                      time=det_info[0], 
                                      azimuth=np.round(det_info[3], 2), 
                                      f_stat=np.round(det_info[5], 4), 
                                      array_d=int(channel_cnt),
                                      f_range=freq_band,
                                      start_end=(det_info[1], det_info[2]),
                                      note=note,
                                      traceV=np.round(det_info[4],2)
                                      )


def write_events(events, event_qls, det_list, local_event_label):
    """
    Write detections from event ID analysis into individual output files

    # TODO: event_qls isn't used in this function. Are we saving it for later?


    Parameters
    ----------
    events: iterable
        List of event labels (detection indices)
    event_qls: iterable
        Event cluster qualities (not currently used, not sure how to include in output .dets.json files)
    det_list: list
        List of infrapy.propagation.likelihoods.InfrasoundDetection instances for the full analysis
    local_event_label: str
        Path for output file(s)

    """
    for ev_n, ev in enumerate(events):
        temp = []
        for det_id in ev:
            temp = temp + [det_list[det_id]]
        detection_list_to_json(local_event_label + "-ev" + str(ev_n) + ".dets.json", temp)


def write_locs(bisl_results, local_loc_label):
    """
    Write localization results to file

    Parameters
    ----------
    bisl_results: dict
        Dictionary of Bayesian Infrasonic Source Localization (BISL) results
    local_loc_label: str
        Path for output file

    """

    if ".loc.json" in local_loc_label:
        with open(local_loc_label, 'w') as of:
            json.dump(bisl_results, of, indent=4, cls=Infrapy_Encoder)
    else:
        with open(local_loc_label + ".loc.json", 'w') as of:
            json.dump(bisl_results, of, indent=4, cls=Infrapy_Encoder)


def read_locs(local_loc_label):
    """
    Ingest a localization result (likely for visualization)

    Parameters
    ----------
    local_loc_label: str
        Path for file

    """

    if ".loc.json" in local_loc_label:
        return json.load(open(local_loc_label))
    else:
        return json.load(open(local_loc_label + ".loc.json"))


def export_beam_results_to_csv(filename, t, f_stats, back_az, trace_v):
    """
    Export the results of the beamforming operation to a csv file for external analysis/plotting
    
    # t, f_stats, back_az, and trace_v are all lists, and they must be the same length

    Parameters
    ----------
    filename: str
        Path for file
    t: iterable
        Analysis times
    f_stats: iterable
        Fisher statistic values
    back_az: iterable
        Back azimuth values
    trace_v: iterable
        Trace velocity values

    """

    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(["Datetime", "Fstat", "TraceV", "BackAz"])
        for idx, t in enumerate(t):
            writer.writerow(t[idx], f_stats[idx], trace_v[idx], back_az[idx])


def export_waveform_to_csv(filename, time, waveform_data):
    """
    Export the timeseries data to a csv file for external analysis/plotting

    # t and data are lists, and they must be the same length

    Parameters
    ----------
    filename: str
        Path for file
    time: iterable
        Waveform times
    waveform_data: iterable
        Waveform values (e.g., overpressure)

    """

    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow("DateTime", "Waveform")
        for t, data in zip(time, waveform_data):
            writer.writerow(t, data)


# ####################################### #
#        Load InfrasoundDetections        #
#           From File                     #
# ####################################### #
def file2dets(file_name):
    """
    Load detection info from a flat (ascii) file

    Parameters
    ----------
    filename: str
        Path for file
    """

    det_list = []
    input = np.genfromtxt(file_name, dtype=None)
    for line in input:
        det_list += [lklhds.InfrasoundDetection(line[0], line[1], np.datetime64(line[2].astype(str)), line[3], line[4], line[5])]

    return det_list


# ############################# #
#   Save detections to a json   #
#   file                        #
# ############################# #
class Infrapy_Encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int64):
            return int(obj)
        elif isinstance(obj, np.float64):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, str):
            return str(obj)
        else:
            return str(obj)


def detection_list_to_json(filename, detections, stream_info=None):
    """
    Write detection info into a .dets.json file

    Parameters
    ----------
    filename: str
        Path for file
    detections: list
        List of infrapy.propagation.likelihoods.InfrasoundDetection instances
    stream_info: list
        Network, station, and channel info
    """

    output = []
    for entry in detections:
        output.append(entry.generateDict())
        if stream_info:
            output[-1]['Network'] = stream_info[0]
            output[-1]['Station'] = stream_info[1]
            output[-1]['Channel'] = stream_info[2]

    with open(filename, 'w') as of:
        json.dump(output, of, indent=4, cls=Infrapy_Encoder)


# ############################# #
#   Load detections from a json   #
#   file                        #
# ############################# #


def json_to_detection_list(filename):
    """
    Read detection info from a .dets.json file

    Parameters
    ----------
    filename: str
        Path for file
    """

    detection_list = []
    with open(filename, 'r') as infile:
        newdata = json.load(infile)
        for entry in newdata:
            detection = lklhds.InfrasoundDetection()
            detection.fillFromDict(entry)
            detection_list.append(detection)
    return detection_list


# ############################# #
#   Load detections from    #
#   database processing         #
# ############################# #

def db2dets(file_name):
    """
    Read detection info from the database (not working yet...not sure we need it really)

    Parameters
    ----------
    filename: str
        Path for file
    """
    det_list = []
    for line in file_name:
        det_list += [lklhds.InfrasoundDetection(line[0], line[1], np.datetime64(UTCDateTime(line[2])), line[3], line[4], line[5])]

    return det_list
