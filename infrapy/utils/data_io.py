#!/usr/bin/env python

import numpy as np

from obspy.clients.fdsn import Client
from obspy.core import read as obspy_read




def wvfms_from_fdsn(fdsn_opt, network, station, location, channel, starttime, endtime):
    # Define FDSN client
    client = Client(fdsn_opt)

    # Get waveforms
    stream = client.get_waveforms(network, station, location, channel, starttime, endtime)

    # Get array geometry
    inventory = client.get_stations(network=network, station=station)
    latlon = None

    return stream, latlon


def wvfrms_from_db(db_url, network, station, location, channel, starttime, endtime):
    return 0


def set_stream(local_opt, fdsn_opt, db_opt, network=None, station=None, location=None, channel=None, starttime=None, endtime=None):
    # check that only one option is selected and issue warning if multiple data sources are specified



    # Check data option and populate obspy Stream
    if local_opt is not None:
        print('\n' + "Loading local data.")
        stream = obspy_read(local_opt)
        latlon = None

    elif fdsn_opt is not None:
        print('\n' + "FDSN methods not set up yet...")
        # stream, latlon = wvfms_from_fdsn(fdsn_opt, network, station, location, channel, starttime, endtime)
        stream, latlon = None, None

    elif db_opt is not None:
        print('\n' + "Database methods not set up yet...")
        # stream, latlon = wvfrms_from_db(db_opt, network, station, location, channel, starttime, endtime)
        stream, latlon = None, None

    else:
        print("No data source specified.")
        stream, latlon = None, None

    return stream, latlon