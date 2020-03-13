#!/usr/bin/env python
import sys, pdb

from sqlalchemy.orm import Session
from sqlalchemy.ext.declarative import declarative_base
#from pisces.io.trace import read_waveform
from obspy.core import UTCDateTime
from datetime import datetime

import pisces as ps


import infrapy.database.schema as schema
import argparse
from IPython import embed

class Fk_params(schema.fk_params):
    __tablename__ = 'FK_PARAMS'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Read combination of parameters used for FK processing")
    parser.add_argument('-d', dest='sq',required=True,help="name of the database connection, e.g.: -d sqlite:///mydb.sqlite ")
    args = parser.parse_args()
    if args.sq:
        sq=args.sq

    session=ps.db_connect(sq)

    try:
        fk_par=session.query(Fk_params).all()
    except Exception as ex1:
        class Fk_params(schema.fk_params):
            __tablename__ = 'FK_PARAMS'

        fk_par=session.query(Fk_params).all()

    for fk_i in fk_par:
        print('pfkid:',str(fk_i.pfkid), 'domain:', str(fk_i.domain), 'freqmin:',str(fk_i.freqmin), 'freqmax:', str(fk_i.freqmax), 'window length:', str(fk_i.beamwinlen), 'window step:', str(fk_i.beamwinstep))

    sys.exit(0)
