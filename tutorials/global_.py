from sqlalchemy import MetaData
from sqlalchemy.ext.declarative import declarative_base
import pisces.schema.kbcore as kb

# make a base that targets global schema
GlobalBase = declarative_base(metadata=MetaData(schema='global'))

class Affiliation(GlobalBase, kb.Affiliation):
    __tablename__ = 'AFFILIATION'

class Amplitude(GlobalBase, kb.Amplitude):
    __tablename__ = 'AMPLITUDE'

class Arrival(GlobalBase, kb.Arrival):
    __tablename__ = 'ARRIVAL'

class Assoc(GlobalBase, kb.Assoc):
    __tablename__ = 'ASSOC'

class Event(GlobalBase, kb.Event):
    __tablename__ = 'EVENT'

class Gregion(GlobalBase, kb.Gregion):
    __tablename__ = 'GREGION'

class Instrument(GlobalBase, kb.Instrument):
    __tablename__ = 'INSTRUMENT'

class Lastid(GlobalBase, kb.Lastid):
    __tablename__ = 'LASTID'

class Netmag(GlobalBase, kb.Netmag):
    __tablename__ = 'NETMAG'

class Network(GlobalBase, kb.Network):
    __tablename__ = 'NETWORK'

class Origerr(GlobalBase, kb.Origerr):
    __tablename__ = 'ORIGERR'

class Origin(GlobalBase, kb.Origin):
    __tablename__ = 'ORIGIN'

class Remark(GlobalBase, kb.Remark):
    __tablename__ = 'REMARK'

class Sensor(GlobalBase, kb.Sensor):
    __tablename__ = 'SENSOR'

class Site(GlobalBase, kb.Site):
    __tablename__ = 'SITE'

class Sitechan(GlobalBase, kb.Sitechan):
    __tablename__ = 'SITECHAN'

class Sregion(GlobalBase, kb.Sregion):
    __tablename__ = 'SREGION'

class Stamag(GlobalBase, kb.Stamag):
    __tablename__ = 'STAMAG'

class Wfdisc(GlobalBase, kb.Wfdisc):
    """
    Waveform description and pointer table.

    Generally shorter, segmented waveforms.

    """
    __tablename__ = 'WFDISC'

class Wfdisc_raw(GlobalBase, kb.Wfdisc):
    """
    Contains raw, generally longer waveforms.

    """
    __tablename__ = 'WFDISC_RAW'

class Wftag(GlobalBase, kb.Wftag):
    __tablename__ = 'WFTAG'
