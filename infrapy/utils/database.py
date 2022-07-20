"""
  put database related functions here
  generic db functions ONLY.  No LANL specific code at all
"""

import numpy as np
import fnmatch 

import sqlalchemy as sa
from sqlalchemy.ext.declarative import declarative_base

import pisces as ps
import pisces.schema.kbcore as kb
import pisces.tables.css3 as css_tables
import pisces.tables.kbcore as kb_tables
import pandas as pd

from obspy import Stream, UTCDateTime

DIALECT_LIST = ['oracle', 'mysql', 'mssql', 'sqllite', 'postgresql']

def db_connect_url(url):
    """
    connect to a database to do database things...

    Parameters
    ----------
    url: str
        Properly formed string containing the connection url for the database

    Returns
    -------
    session : bound SQLAlchemy session instance
    """
    session = ps.db_connect(url)
    return session


def db_connect(dialect, hostname, db_name, port=None, username="", password="", driver=""):
    '''
        Connect to a database to do database things...

        Parameters
        ----------
        dialect: str
            Type of database.   (examples: oracle, mysql)
        hostname : str
            The url of the database. (example: mydb.home.org)local branch

        Returns
        -------
        session : bound SQLAlchemy session instance

    '''
    return ps.db_connect(assemble_db_url(dialect, hostname, db_name, port, username, password, driver))


def assemble_db_url(dialect, hostname, db_name, port=None, username="", password="", driver=""):
    '''
        Assemble a database connection url given the supplied parts.
        All inputs are strings
    '''
    if driver:
        driver = '+' + driver
    return dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + port + "/" + db_name
 

def db_connect2(db_info):
    dialect = db_info['DATABASE']['dialect']
    hostname = db_info['DATABASE']['hostname']
    db_name = db_info['DATABASE']['database_name']
    port = db_info['DATABASE']['port']

    try:
        username = db_info['DATABASE']['username']
    except:
        username = ''

    try:
        password = db_info['DATABASE']['password']
    except:
        password = ''

    try:
        driver = db_info['DATABASE']['driver']
    except:
        driver = ''

    return ps.db_connect(assemble_db_url(dialect, hostname, db_name, port, username, password, driver))

def check_connection(session):
    """
    Simple function to check that there is a valid connection by 
    calling engine.connect() to see if it returns true
    """
    engine = session.get_bind()
    try:
        engine.connect()
        return True
    except Exception as e:
        print(e)
        return False

'''
    class Site(kb.Site):
        __tablename__ = 'global.site'


    class Wfdisc(kb.Wfdisc):
        __tablename__ = 'global.wfdisc_raw'
        

    class Affiliation(kb.Affiliation):
        __tablename__ = 'global.affiliation'

    class Origin(kb.Origin):
        __tablename__ = 'global.origin'

    class Event(kb.Event):
        __tablename__ = 'global.event'
'''

def query_db(session, start_time, end_time, sta="%", cha="%", return_type='dataframe'):
    """
    function to query a database using an existing session.  

    return_type values can be 'dataframe' for a pandas dataframe, or 'wfdisc_rows' for wfdisc rows
    """
    db_tables = make_tables_from_dict(tables={'Wfdisc':'wfdisc_raw', 'Site': 'site'}, schema='kbcore', owner='global')

    # db_tables['Site'].__table__

    if session is None:
        return None
    
    print("return type = {}".format(return_type))
    if return_type == 'dataframe':
        my_query = session.query(db_tables['Wfdisc']).filter(db_tables['Wfdisc'].sta == sta)\
                                        .filter(db_tables['Wfdisc'].time < end_time.timestamp)\
                                        .filter(db_tables['Wfdisc'].endtime > start_time.timestamp)\
                                        .filter(db_tables['Wfdisc'].chan.like(cha))

        return pd.read_sql(my_query.statement, session.bind)
    elif return_type == 'wfdisc_rows':
        return ps.request.get_wfdisc_rows(session, db_tables['Wfdisc'], sta=sta, t1=start_time, t2=end_time)
    else:
        return None

def wvfrms_from_db(db_info, stations, channel, starttime, endtime):
    # Set up db connection and table info
    # session = ps.db_connect( db_info['url'])
    session = db_connect2(db_info)
    
    db_tables = make_tables_from_dict(tables=db_info['DBTABLES'], schema=db_info['DATABASE']['schema'], owner=db_info['DATABASE']['owner'])

    Site = db_tables['site']
    Wfdisc = db_tables['wfdisc']

    # convert station wildcards to SQL and check that channel is not None
    if type(stations) is str:
        stations = stations.replace('*','%')

    if channel is None:
        channel = "*"

    # get station info
    if "%" in stations:
        # Load data specified with a while card (e.g., 'I26H*') via a Site table query
        sta_list = session.query(Site).filter(Site.sta.contains(stations))
    elif ',' in stations:
        # Load data specified by a string list of stations (e.g., 'I26H1, I26H2, I26H3, I26H4') with get_stations
        sta_list = ps.request.get_stations(session, Site, stations=stations.strip(' ()[]').split(','))
    else:
        # Load data specified by a Python list of strings (e.g., ['I26H1', 'I26H2', 'I26H3', 'I26H4']) with get_stations
        sta_list = ps.request.get_stations(session, Site, stations=stations)

    # pull data into the stream and merge to combine time segments
    st = Stream()
    for sta_n in sta_list:
        temp_st = ps.request.get_waveforms(session, Wfdisc, station=sta_n.sta, starttime=UTCDateTime(starttime).timestamp, endtime=UTCDateTime(endtime).timestamp)
        for tr in temp_st:
            if  fnmatch.fnmatch(tr.stats.channel, channel.replace("%","*")):
                tr.stats.sac = {'stla': sta_n.lat, 'stlo': sta_n.lon}
                if len(tr.stats.network) == 0:
                    tr.stats.network = "__"
                st.append(tr)
    st.merge()
    st.split()
    
    # Set the latlon info
    latlon = [[tr.stats.sac['stla'], tr.stats.sac['stlo']] for tr in st]

    return st, latlon


def make_tables_from_dict(tables=None, schema=None, owner=None):
    # first handle the bailout conditions
    if tables is None and schema is None:
        msg = "Not enough information to generate tables"
        raise ValueError(msg)
        return
    if schema.lower() not in ['kbcore', 'css3']:
        msg = "Unsupported schema: {}".format(schema)
        raise ValueError(msg)
        return

    if schema.lower() == 'kbcore':
        core_tables = kb_tables.CORETABLES
    elif schema.lower() == 'css3' or schema.lower() == 'css':
        core_tables = css_tables.CORETABLES
    
    if tables is None:
        return core_tables
    else:
        if owner:
            for key, value in tables.items():
                tables[key] = owner + "." + value

        dict_of_classes = {}
        for table, tablename in tables.items():
            prototype = core_tables[table.lower()].prototype
            dict_of_classes[table] = type(table, (prototype,), {'__tablename__': tablename})
    
    return dict_of_classes


def eventID_query(session, eventID):
    print("Querying for event id: {}".format(eventID))
    events = ps.request.get_events(session, Origin, Event, eventID)
    if events:
        print(events)
    else:
        print("no event found")