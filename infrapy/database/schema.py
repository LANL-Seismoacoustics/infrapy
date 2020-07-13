

from datetime import datetime
from sqlalchemy import Date, DateTime, Float, Numeric, String, Integer
from sqlalchemy import Column, Table, func
from sqlalchemy import PrimaryKeyConstraint, UniqueConstraint, Index
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.declarative import declared_attr

from obspy.core import UTCDateTime
from pisces.schema.util import PiscesMeta
from pisces.schema.util import parse_int, parse_float, parse_str

from pisces.io.trace import wfdisc2trace


import sqlalchemy as sa
import pisces.schema.css3 as schem
#import pisces.schema.kbcore as scheme


algorithm = Column(String(15),info={'default': '-', 'parse': parse_str, 'dtype': 'a15', 'width': 15, 'format': '15.15s'})
azimuth = Column(Float(53), info={'default': -1.0, 'parse': parse_float, 'dtype': 'float', 'width': 7, 'format': '7.2f'})
id = sa.Column(Integer,info={'default': -1, 'parse': parse_int, 'dtype': 'int', 'width': 9, 'format': '9d'})
name = Column(Float(53),info={'default': -9999999999.999, 'parse': parse_float, 'dtype': 'float', 'width': 17, 'format': '17.5f'})
overlapwlen = sa.Column(Integer,info={'default': -1, 'parse': parse_int, 'dtype': 'int', 'width': 9, 'format': '9d'})
time = Column(Float(53),info={'default': -9999999999.999, 'parse': parse_float, 'dtype': 'float', 'width': 17, 'format': '17.5f'})
wlen = sa.Column(Float(53), info={'default': -1, 'parse': parse_int, 'dtype': 'int', 'width': 8, 'format': '8d'})
gen_float =sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
dist = Column(Float(53), info={'default': -999.0, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
geoloc = Column(Float(53), info={'default': -999.0, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
net = Column(String(8),info={'default': '-', 'parse': parse_str, 'dtype': 'a8', 'width': 8, 'format': '8.8s'})
chan = Column(String(8),info={'default': '-', 'parse': parse_str, 'dtype': 'a8', 'width': 8, 'format': '8.8s'})


class noise_params(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('pnoiseid'),sa.UniqueConstraint('overlapwlen','filterl', 'filterh', 'wlen','name','stepslowness','algorithm','fkfreql','fkfreqh'))

    pnoiseid = id.copy()
    name=name.copy()

    filtertype = sa.Column(String(15), info={'default': '-', 'parse': parse_str, 'dtype': 'a32', 'width': 15, 'format': '15.15s'})
    filterl = sa.Column(Float(24), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    filterh = sa.Column(Float(24), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    wlen = wlen.copy()
    overlapwlen = overlapwlen.copy()
    n_octave= id.copy()
    algorithm=algorithm.copy()

class EDetectParams(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('pedetectid'),sa.UniqueConstraint('overlapwlen','filterl', 'filterh', 'wlen','algorithm'))

    pedetectid = id.copy()
    name=name.copy()
    filtertype  = sa.Column(String(15), info={'default': '-', 'parse': parse_str, 'dtype': 'a32', 'width': 15, 'format': '15.15s'})
    filterl     = sa.Column(Float, info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    filterh     = sa.Column(Float, info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    wlenshort   = sa.Column(Float, info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    wlenlong    = sa.Column(Float, info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    overlapwlen = sa.Column(Float, info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    algorithm   = algorithm.copy()

class fk_params(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('pfkid'),sa.UniqueConstraint('beamwinlen','freqmin', 'freqmax', 'beamwinstep','stepslowness','algorithm','minslowness','maxslowness','numsources','domain'))

    pfkid = id.copy()
    name=algorithm.copy()

    domain = algorithm.copy()

    filter = sa.Column(String(15), info={'default': '-', 'parse': parse_str, 'dtype': 'a32', 'width': 15, 'format': '15.15s'})
    freqmin = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    freqmax = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    beamwinlen = wlen.copy()
    beamwinstep = overlapwlen.copy()
    backazmin = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    backazmax = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    backazstep = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    trvelmin = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    trvelmax = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    trvelstep = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    numsources = sa.Column(Integer,info={'default': 1, 'parse': parse_int, 'dtype': 'int', 'width': 9, 'format': '9d'})
    minslowness = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    maxslowness = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    stepslowness = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    minslowness = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    maxslowness = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    stepslowness = sa.Column(Float(53), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    numsources = sa.Column(Integer,info={'default': 1, 'parse': parse_int, 'dtype': 'int', 'width': 9, 'format': '9d'})
    additional1=sa.Column(Float(53),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    additional2=sa.Column(Float(53),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    algorithm=algorithm.copy()


class fk_results(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('fkid'),)

    pfkid = id.copy()
    fkid = id.copy()
    chan=chan.copy()
    nchan = id.copy()
    sta=schem.sta.copy()
    timeini=time.copy()
    timeend=time.copy()
    sx=sa.Column(Float(24), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    sy=sa.Column(Float(24), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    esx=sa.Column(Float(24), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    esy=sa.Column(Float(24), info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 5, 'format': '5.3f'})
    bz=azimuth.copy()
    slofk=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    trvel=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    ebz=azimuth.copy()
    eslofk=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    etrvel=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    fval=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    xcorrvalmax=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    xcorrvalmin=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    xcorrvalmean=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    coher=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    rmsval=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    sourcenum = id.copy()
    additional1=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    additional2=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})


class fd_params(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('pfdid'),sa.UniqueConstraint('pthreshold', 'cthr', 'detwinlen','minlen','algorithm','numsources'))

    pfdid = id.copy()
    name=name.copy()
    pthreshold = sa.Column(Float(53),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    cthr = sa.Column(Float(53),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    numsources = sa.Column(Integer,info={'default': 1, 'parse': parse_int, 'dtype': 'int', 'width': 9, 'format': '9d'})
    adapwlen = wlen.copy()
    minlen = wlen.copy()
    #detthresh = sa.Column(Float(53),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    dsegmin = sa.Column(Float(53),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    numsources = sa.Column(Integer,info={'default': 1, 'parse': parse_int, 'dtype': 'int', 'width': 9, 'format': '9d'})
    tb_prod = sa.Column(Integer,info={'default': 1, 'parse': parse_int, 'dtype': 'int', 'width': 9, 'format': '9d'})
    backazlim = sa.Column(Integer,info={'default': 1, 'parse': parse_int, 'dtype': 'int', 'width': 9, 'format': '9d'})
    detwinlen = wlen.copy()
    algorithm=algorithm.copy()
    additional1=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    additional2=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})

class fd_results(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('fdid'),sa.UniqueConstraint('fdid','pfdid', 'pfkid'))

    fdid = id.copy()
    pfdid = id.copy()
    pfkid = id.copy()
    sta=schem.sta.copy()
    sourcenum = id.copy()
    timeini=time.copy()
    timeend=time.copy()
    tonset=time.copy()
    tend=time.copy()
    #fval=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    #bz=azimuth.copy()
    trvel=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    ebz=azimuth.copy()
    etrvel=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    #fval=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    c=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    maxfc=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    etimeini=time.copy()
    etimeend=time.copy()
    ec=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    maxfc=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    #maxpfc=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    maxfo=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    #c=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    fktablename=Column(String(15),info={'default': '-', 'parse': parse_str, 'dtype': 'a15', 'width': 15, 'format': '15.15s'})


class ASSOC_params(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('passocid'),sa.UniqueConstraint('algorithm','beamwidth','rangemax','clusterthresh','trimthresh','mindetpop','minarraypop','duration'))
    passocid = id.copy()
    name=name.copy()
    algorithm=algorithm.copy()
    beamwidth=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    rangemax=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    clusterthresh=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    trimthresh = Column(String(15),info={'default': '-', 'parse': parse_str, 'dtype': 'a15', 'width': 15, 'format': '15.15s'})
    trimthreshscalar =sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    mindetpop=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    minarraypop=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})
    duration=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width':9, 'format': '9.4f'})


class ASSOC_results(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('associd'),sa.UniqueConstraint('net', 'fdid','associd','passocid','timeini','timeend'))
    associd = id.copy()
    passocid = id.copy()
    eventid = id.copy()
    fdid = id.copy()
    fdtable = Column(String(15),info={'default': '-', 'parse': parse_str, 'dtype': 'a15', 'width': 15, 'format': '15.15s'})
    sta=schem.sta.copy()
    qdetcluster=gen_float.copy()
    qassoc=gen_float.copy()
    net=net.copy()
    timeini=time.copy()
    timeend=time.copy()

class loc_params(schem.Base):
    __abstract__ = True

    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('plocid'),sa.UniqueConstraint('maxlat','minlat','maxlon','minlon','priors','algorithm'))
    plocid = id.copy()
    name=name.copy()
    maxlat=geoloc.copy()
    minlat=geoloc.copy()
    maxlon=geoloc.copy()
    minlon=geoloc.copy()
    priors=id.copy()
    algorithm=algorithm.copy()

class loc_results(schem.Base):
    __abstract__ = True
    @declared_attr
    def __table_args__(cls):
        return (sa.PrimaryKeyConstraint('locid'),sa.UniqueConstraint('plocid','eventid','timeini','timeend'))

    locid = id.copy()
    plocid = id.copy()
    eventid = id.copy()
    net = net.copy()
    numstations=id.copy()
    timeorigmap=time.copy()
    timeorigmean=time.copy()
    timeorigvar=time.copy()
    timeorigmin=time.copy()
    timeorigmax=time.copy()

    latorigmean=geoloc.copy()
    latorigvar=geoloc.copy()

    lonorigmean=geoloc.copy()
    lonorigvar=geoloc.copy()

    latlonorigcovar=geoloc.copy()

    lonorigmap=geoloc.copy()
    latorigmap=geoloc.copy()

    par_a=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    par_b=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    par_theta=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})

    timeini=time.copy()
    timeend=time.copy()

    additional1=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    additional2=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    additional3=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
    additional4=sa.Column(Float(24),info={'default': -1, 'parse': parse_float, 'dtype': 'float', 'width': 9, 'format': '9.4f'})
