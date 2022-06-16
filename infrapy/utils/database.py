"""
  put database related functions here
  generic db functions ONLY.  No LANL specific code at all
"""

import numpy as np
import fnmatch 

import pisces as ps
import pisces.schema.kbcore as kb
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
    """
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

    """

    if driver:
        driver = '+' + driver

    connect_str = dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + port + "/" + db_name

    session = ps.db_connect(connect_str)

    return session


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


class Site(kb.Site):
    __tablename__ = 'global.site'


class Wfdisc(kb.Wfdisc):
    __tablename__ = 'global.wfdisc_raw'
    

class Affiliation(kb.Affiliation):
    __tablename__ = 'global.affiliation'


def query_db(session, start_time, end_time, evid="", sta="%", cha="%", return_type='dataframe'):
    """
    function to query a database using an existing session.  

    return_type values can be 'dataframe' for a pandas dataframe, or 'wfdisc_rows' for wfdisc rows
    """

    if session is None:
        return None

    my_query = session.query(Wfdisc).filter(Wfdisc.sta == sta)\
                                    .filter(Wfdisc.time < end_time.timestamp)\
                                    .filter(Wfdisc.endtime > start_time.timestamp)\
                                    .filter(Wfdisc.chan.like(cha))

    if return_type == 'dataframe':
        return pd.read_sql(my_query.statement, session.bind)
    elif return_type == 'wfdisc_rows':
        return ps.request.get_wfdisc_rows(session, Wfdisc, sta=sta, t1=start_time, t2=end_time)
    else:
        return None


def wvfrms_from_db(db_info, stations, channel, starttime, endtime):
    # Set up db connection and table info
    session = ps.db_connect( db_info['url'])

    class Site(kb.Site):
        __tablename__ = db_info['site']

    class Wfdisc(kb.Wfdisc):
        __tablename__ = db_info['wfdisc']

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
