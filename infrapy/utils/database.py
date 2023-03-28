"""
    put database related functions here
    generic db functions ONLY.  No LANL specific code at all
"""


import numpy as np
import fnmatch 
import os

import sqlalchemy as sa
from sqlalchemy.orm import Session

import pisces as ps
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


def db_connect(dialect="", hostname="", db_name="", port="", username="", password="", driver=""):
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
    # my_dialect = dialect + "://"
    # url = sa.engine.url.make_url(my_dialect)
    # print(url)
    # url.username = username
    # # url.drivername = driver
    # url.password = password
    # url.host = hostname
    # url.port = port
    # url.database = db_name

    # print(url)
    # engine = sa.create_engine(url)
    # return Session(bind=engine)

    return ps.db_connect(assemble_db_url(dialect, hostname, db_name=db_name, port=port, username=username, password=password, driver=driver))


def assemble_db_url(dialect, hostname, db_name, port=None, username="", password="", driver=""):
    '''
        Assemble a database connection url given the supplied parts.
        All inputs are strings
    '''
    if driver:
        driver = '+' + driver
    return dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + port + "/" + db_name
 
def set_db_env_variables(env_vars):
    # env_vars is a dictionary containing the environment variable as the key, and what to set it to as the value
    # Note that the environment variables will only last for the duration of the session.
    for key, value in env_vars.items():
            os.environ[key] = value

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
    try:
        session.get_bind().connect()
        return True
    except Exception as e:
        return False


def query_db(session, tables, start_time, end_time, sta="%", cha="%", return_type='dataframe', asquery=False):

    db_tables = make_tables_from_dict(tables=tables, schema='kbcore')

    if session is None:
        return None
    
    if return_type == 'dataframe':
        my_query = session.query(db_tables['wfdisc']).filter(db_tables['wfdisc'].sta.like(sta)\
                                        .filter(db_tables['wfdisc'].time < end_time.timestamp)\
                                        .filter(db_tables['wfdisc'].endtime > start_time.timestamp)\
                                        .filter(db_tables['wfdisc'].chan.like(cha)))
        if asquery:
            return my_query.statement
        else:
            return pd.read_sql(my_query.statement, session.bind)

    elif return_type == 'wfdisc_rows':
        my_query =  ps.request.get_wfdisc_rows(session, db_tables['wfdisc'], chan=cha, t1=start_time, t2=end_time, asquery=True)
        my_query = my_query.filter(db_tables['wfdisc'].sta.like(sta))
        if asquery:
            return my_query
        else:
            return my_query.all()
    else:
        return None

def prep_session(db_info, check_connection=False):
    if 'url' in db_info['DATABASE'].keys():
        session = ps.db_connect( db_info['DATABASE']['url'])
    else:
        session = db_connect2(db_info)

    if check_connection:
        try:
            session.get_bind().connect()
            print("Database connection check passed")
        except Exception as e:
            print("Database connection check failed")
        
    db_tables = make_tables_from_dict(tables=db_info['DBTABLES'], schema=db_info['DATABASE']['schema'])

    return session, db_tables


def wvfrms_from_db(session, db_tables, stations, channel, starttime, endtime):
    ''' 
        function to pull obspy streams from the database.  
        stations: str
        channel: str
        starttime: UTCDateTime 
    '''

    Site = db_tables['site']
    wfdisc = db_tables['wfdisc']

    # convert station wildcards to SQL and check that channel is not None
    if type(stations) is str:
        stations = stations.replace('*','%')

    if channel is None:
        channel = "*"

    # the pisces.request.get_stations function requires the start/end days to be an integer in jdate format (YYYYDDD)
    # so we have to massage our starttime and endtimes to that form
    julian_start = starttime.year * 1000 + starttime.julday
    julian_end = endtime.year * 1000 + endtime.julday
    wtime = (julian_start, julian_end)

    # get station info
    # if "%" in stations:
    #     # Load data specified with a while card (e.g., 'I26H*') via a Site table query
    #     sta_list = session.query(Site).filter(Site.sta.contains(stations))
    # elif ',' in stations:
    #     # Load data specified by a string list of stations (e.g., 'I26H1, I26H2, I26H3, I26H4') with get_stations
    #     sta_list = ps.request.get_stations(session, Site, stations=stations.strip('()[]').replace(" ", "").split(','))
    # else:
    #     # Load data specified by a Python list of strings (e.g., ['I26H1', 'I26H2', 'I26H3', 'I26H4']) with get_stations
    #     sta_list = ps.request.get_stations(session, Site, stations=stations)

    sta_list = ps.request.get_stations(session, Site, stations=stations, time_span=wtime)

    # pull data into the stream and merge to combine time segments
    st = Stream()
    for sta_n in sta_list:
        temp_st = ps.request.get_waveforms(session, wfdisc, station=sta_n.sta, starttime=UTCDateTime(starttime).timestamp, endtime=UTCDateTime(endtime).timestamp)
        for tr in temp_st:
            tr.data = tr.data - np.mean(tr.data)
            tr.stats['_format'] = 'SAC'
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

def make_tables_from_dict(tables=None, schema=None):
    # first handle the bailout conditions
    if tables is None and schema is None:
        msg = "Not enough information to generate tables"
        raise ValueError(msg)

    if schema.lower() not in ['kbcore', 'css3', 'css']:
        msg = "Unsupported schema: {}".format(schema)
        raise ValueError(msg)

    if schema.lower() == 'kbcore':
        core_tables = kb_tables.CORETABLES
    elif schema.lower() == 'css3' or schema.lower() == 'css':
        core_tables = css_tables.CORETABLES
    
    if tables is None:
        return core_tables
    else:
        dict_of_classes = {}
        for table, tablename in tables.items():
            prototype = core_tables[table.lower()].prototype
            dict_of_classes[table] = type(table, (prototype,), {'__tablename__': tablename})
    
    return dict_of_classes

def eventID_query(session, eventID, db_tables, asquery):
    # session is a current active session
    # eventID is the event id to search for
    # db_tables is a dictionary of available mapped tables (needs to contain Event and Origin tables)

    evIDs = [int(eventID)]  # for now only query one evid at a time
    if asquery:
        return ps.request.get_events(session, db_tables['origin'], event=db_tables['event'], evids=evIDs, asquery=True)
    else:
        events = ps.request.get_events(session, db_tables['origin'], event=db_tables['event'], evids=evIDs)

    if events:
        return events
    else:
        print("no event found")
        return None

def event_query_area(session, center_lat, center_lon, minr, maxr, db_tables):

    if 'Origin' not in db_tables or 'Event' not in db_tables:
        raise KeyError
    
    events = ps.request.get_events(session, db_tables['origin'], db_tables['events'], km=(center_lat, center_lon, minr, maxr), etime=(startt, endt))
