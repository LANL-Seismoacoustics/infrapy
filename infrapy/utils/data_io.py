#!/usr/bin/env python

import warnings 

import numpy as np

from obspy.clients.fdsn import Client
from obspy.core import read as obspy_read
from obspy import UTCDateTime


def wvfms_from_fdsn(fdsn_opt, network, station, location, channel, starttime, endtime):
    client = Client(fdsn_opt)
    stream = client.get_waveforms(network, station, location, channel, UTCDateTime(starttime), UTCDateTime(endtime))
    inventory = client.get_stations(network=network, station=station, starttime=UTCDateTime(starttime), endtime=UTCDateTime(endtime))
    latlon = []
    for sta in inventory[0]:
        latlon = latlon + [[sta.latitude, sta.longitude]]

    return stream, latlon


def wvfrms_from_db(db_url, network, station, location, channel, starttime, endtime):
    # Need Jonathan or Christine to set this up...
    return None, None


def set_stream(local_opt, fdsn_opt, db_opt, network=None, station=None, location=None, channel=None, starttime=None, endtime=None, local_latlon=None):
    # check that only one option is selected and issue warning if multiple data sources are specified
    if np.sum(np.array([val is not None for val in [local_opt, fdsn_opt, db_opt]])) > 1:
        msg = "Multiple data sources specified. Unexpected behavior is possible." + '\n' + "Preferential order is local > FDSN > DB"
        warnings.warn(msg)

    # Check data option and populate obspy Stream
    if local_opt is not None:
        print('\n' + "Loading local data...")
        stream = obspy_read(local_opt)
        if local_latlon:
            latlon = np.load(local_latlon)
        else:
            latlon = None

    elif fdsn_opt is not None:
        print('\n' + "Loading data from FDSN (" + fdsn_opt + ")...")
        stream, latlon = wvfms_from_fdsn(fdsn_opt, network, station, location, channel, starttime, endtime)

    elif db_opt is not None:
        print('\n' + "Loading data from database (methods not set up yet, so returning 'None's)...")
        stream, latlon = wvfrms_from_db(db_opt, network, station, location, channel, starttime, endtime)

    else:
        print("No data source specified.")
        stream, latlon = None, None

    return stream, latlon