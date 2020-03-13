#!/usr/bin/env python
from obspy.core import UTCDateTime

def short_time(timeSTR):
     return UTCDateTime(timeSTR).strftime('%y-%m-%d %H:%M:%S')

