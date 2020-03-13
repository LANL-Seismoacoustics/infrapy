import sys, pdb
import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
from obspy.core import trace
from obspy.core import Stream
from obspy.core.util import AttribDict
from datetime import datetime

import numpy as np
from numpy import append
import pisces as ps
from pisces import request
from IPython import embed

import time



def get_channel(session, Site, Wfdisc, array):
    ssSite = session.query(Site).filter(Site.refsta==array).all()
    wf_T
    for ssi in ssSite:
        wf = session.query(Wfdisc).filter(Wfdisc.sta==ssi).distinct()
        wf_T.append(wf)
    return np.unique(wf_T)

def get_arraywaveforms(session, Site, Wfdisc, array, t0=None, te=None, channel=None):
    ssSite = session.query(Site).filter(Site.refsta==array).all()
    wf = session.query(Wfdisc)
    aa = Stream()

    for ssi in ssSite:
        if t0 == None:
            wfT = wf.filter(Wfdisc.sta==ssi.sta).all()
            timeI = []
            timeF = []
            for wfTi in wfT:
                timeI.append(wfTi.time)
                timeF.append(wfTi.endtime)
                timeI = np.asanyarray(timeI)
                timeF = np.asanyarray(timeF)
                t0 = max(timeI)
                te = min(timeF)

        for t1 in range(5):
            try:
                aaT = request.get_waveforms(session, Wfdisc, station=ssi.sta, starttime=t0,
                                            endtime=te,channel=channel)
                break
            except:
                print('try get data:',t1)
                print('go to sleep for 5 seconds and try again')
                time.sleep(5)

            if t1 == 4:
                print('There is problem connecting to the data waveforms')
                exit()

        if len(aaT) == 0:
            print('maybe this is a ref name, sta:',ssi.sta)
            continue

        aaT.merge(fill_value=0)
        #if not len(aaT) == 1:
        #    print('there is a problem with data retrieving; there is more than one trace for this station')
        #    sys.exit(0)

        aaT[0].stats.coordinates = AttribDict({'latitude': ssi.lat,'elevation': ssi.elev,'longitude': ssi.lon})
        aa = aa + aaT

    aa.merge(fill_value=0)

    return aa
