#!/usr/bin/env python

import os
import warnings
import fnmatch
import json
import gzip
import csv

import numpy as np

from obspy.clients.fdsn import Client
from obspy import read as obspy_read
from obspy import UTCDateTime

from ..propagation import likelihoods as lklhds
from . import database


#########################
##  Meta Data Methods  ##
#########################

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


blank_sac_dict = {'delta': None,
                  'npts': None,
                  'depmin': None,
                  'depmax': None,
                  'depmen': None,
                  'b': 0.0,
                  'e': None,
                  'stla': None,
                  'stlo': None,
                  'nzyear': None,
                  'nzjday': None,
                  'nzhour': None,
                  'nzmin': None,
                  'nzsec': None,
                  'nzmsec': None,
                  'kstnm': None,
                  'kcmpnm': None,
                  'knetwk': None}


def stream_label(st):
    label = os.path.commonprefix([tr.stats.network for tr in st])
    label = label + "." + os.path.commonprefix([tr.stats.station for tr in st])
    label = label + '_' + st[0].stats.starttime.strftime('%Y.%m.%dT%H.%M.%S')

    return label


def wvfrm_info(st, latlon):
    info = []
    for n, tr in enumerate(st):
        info = info + [{}]
        info[-1]['trace id'] = tr.id
        info[-1]['starttime'] = str(tr.stats.starttime)
        info[-1]['endtime'] = str(tr.stats.endtime)
        info[-1]['latitude'] = float(latlon[n][0])
        info[-1]['longitude'] = float(latlon[n][1])

    return info



##############################
##  Data Ingestion Methods  ##
##############################
def wvfrms_from_fdsn(fdsn_opt, network, station, location, channel, starttime, endtime):
    """
    Connect to an FDSN server and pull waveform data into an ObsPy Stream

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
    t1 = UTCDateTime(starttime)
    t2 = UTCDateTime(endtime)

    stream = client.get_waveforms(network, station, location, channel, t1, t2)
    stream.merge(fill_value=0)

    inventory = client.get_stations(network=network, station=station, location=location, channel=channel, starttime=t1, endtime=t2, level="response")
    stream.remove_response(inventory=inventory, output="DEF")

    latlon = []
    for tr in stream:
        coords = inventory.get_coordinates(tr.get_id(), UTCDateTime(tr.stats.starttime))
        latlon = latlon + [[coords['latitude'], coords['longitude']]]

    return stream, latlon


def set_stream(local_opt, fdsn_opt, db_info, network=None, station=None, location=None, channel=None, starttime=None, endtime=None, local_latlon=None):
    """
    Define an ObsPy stream from a specified local, FDSN, or database source.
    1) if specifying local data, use obspy.read to set up the stream
    2) if pulling from an FDSN, use obspy.clients.fdsn.Client to pull waveforms and station info
    3) if pulling from a database...this needs to be updated

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

    # if local data is specified, load using ObsPy's read function
    if local_opt is not None:
        print('\n' + "Loading local data from " + local_opt)
        stream = obspy_read(local_opt)
        if local_latlon:
            latlon = np.load(local_latlon)
        else:
            latlon = [[tr.stats.sac['stla'], tr.stats.sac['stlo']] for tr in stream]

    # if FDSN, pass to above function
    elif fdsn_opt is not None:
        print('\n' + "Loading data from FDSN (" + fdsn_opt + ")...")
        stream, latlon = wvfrms_from_fdsn(fdsn_opt, network, station, location, channel, starttime, endtime)

    ## if database extraction, pass to database methods
    elif db_info is not None:
        print('\n' + "Loading data from database...")
        session, db_tables = database.prep_session(db_info)
        stream, latlon = database.wvfrms_from_db(session, db_tables, station, channel, UTCDateTime(starttime), UTCDateTime(endtime))

    # return error if no waveform data was specified
    else:
        msg = "Warning: No waveform data source specified."
        warnings.warn(msg)
        stream, latlon = None, None

    return stream, latlon


def _load_dets_json(dets_files):
    """
    Read in multiple [...].dets.json(.gz) files specified by either a comma
    separated list or a wild card glob

    Parameters
    ----------
    dets_files: str
        Detections file(s), [...].dets.json.gz, to be ingested

    Returns
    -------
    det_list : list of dictionary instances
        Iterable list of dictionaries containing detection information

    """

    def temp_open_json(dets_file):
        if os.path.splitext(dets_file)[-1] == ".gz":
            return json.load(gzip.open(dets_file, 'rt'))
        else:
            return json.load(open(dets_file))

    det_list = []
    if "*" not in dets_files:
        # define file list from comma or space separated list
        for file in dets_files.replace(" ","").split(","):
            det_list = det_list + [temp_open_json(file)]
    else:
        # define file list from wild card glob
        if "/" in dets_files:
            file_path = os.path.dirname(dets_files) + "/"
            dir_files = os.listdir(os.path.dirname(dets_files))
        else:
            file_path = ""
            dir_files = os.listdir(".")
        dir_files = np.sort(dir_files)

        file_list = []
        for file in dir_files:
            if fnmatch.fnmatch(file, os.path.basename(dets_files)):
                file_list += [file]

        if len(file_list) == 0:
            msg = '\n' + "Detection file(s) specified not found"
            warnings.warn(msg)
            det_list = None
        else:
            for file in file_list:
                det_list = det_list + [temp_open_json(file_path + file)]

    return det_list


def _det_dict_to_likelihood(det_dict):
    """
    Initialize a infrapy.propagation.likelihoods.InfrasoundDetection instance
    and load from an extracted dictionary into it for use in BISL, TRIBL, or SpYE

    Parameters
    ----------
    det_dict: dictionary
        Dictionary containing InfraPy detection information


    Returns
    -------
    detection: infrapy.propagation.liklihoods.InfrasoundDetection instance
        InfrasoundDetection instance containing the information from the dictionary

    """

    detection = lklhds.InfrasoundDetection()
    detection.fillFromDict2(det_dict)

    return detection


############################
##  Data Writing Methods  ##
############################

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

    labels = [tr.id for tr in stream]
    if len(np.unique(labels)) < len(stream):
        print("Warning!  Non-unique labels.  Adding indexing...")
        labels = [label + "-" + str(n) for n, label in enumerate(labels)]

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

        label = labels[m] + tr.stats.starttime.strftime('_%Y.%m.%d_%H.%M.%S')

        tr.write(label + ".sac", format='SAC')


#######################
##  InfraView Write. ##
##   to CSV Methods  ##
#######################
def export_beam_results_to_csv(filename, time, f_stats, back_az, trace_v):
    """
    Export the results of the beamforming operation to a csv file for external analysis/plotting

    # t, f_stats, back_az, and trace_v are all lists, and they must be the same length

    Parameters
    ----------
    filename: str
        Path for file
    time: iterable
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
        for t, fs, tv, ba in zip(time, f_stats, back_az, trace_v):
            writer.writerow([t, fs, tv, ba])


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
        writer.writerow(["DateTime", "Waveform"])
        for t, data in zip(time, waveform_data):
            writer.writerow([t, data])







