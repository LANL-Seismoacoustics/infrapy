"""
    put database related functions here
    generic db functions ONLY.  No LANL specific code at all
"""


import numpy as np
import fnmatch
import os

import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy import create_engine

import pisces as ps
import pisces.tables.css3 as css_tables
import pisces.tables.kbcore as kb_tables
import pandas as pd

from obspy import Stream, UTCDateTime

DIALECT_LIST = ['oracle', 'mysql', 'mssql', 'sqllite', 'postgresql']
AUTHENTICATION_LIST = ['password', 'os', 'wallet']

def db_connect(db_info):
    """
        Function to set up a database connection session using a dictionary of connection information.
        The current configuration file format is below ([DATABASE] block passed to this function):

            [DATABASE]
            hostname = [...]
            username = [...]
            database_name = [...]
            port = [...]
            driver = [...]
            dialect = [...] (required)
            alias = [...]
            authentication = [...] (required)
            url = [...]

        Though a dictionary of key/value pairs matching the above keys can also be passed into the function.
        
        Dialect and authentication fields are always required and must be in the DIALECT_LIST and 
        AUTHENTICATION_LIST as follows:

            DIALECT_LIST = ['oracle', 'mysql', 'mssql', 'sqllite', 'postgresql']
            AUTHENTICATION_LIST = ['password', 'os', 'wallet']

        If a url is not provided in the dictionary block, an attempt will be made to construct one based on 
        the other variables, in which case, hostname, database_name, and port are required unless 
        an alias is provided.  If the authentication method is 'password', a username is also required.

        Additional information for Oracle databases:
            1) Driver is required for Oracle databases and must be one of cx_oracle or oracledb if a url is 
               not provided.  Note that cx_oracle is deprecated in it's entirety in favor of oracledb.
            2) Authentication methods 'os' and 'wallet' are specific to oracle databases.  If your 
               authentication method is 'wallet' and that wallet is pointed to an alias, that alias must be 
               provided or used in the provided url.
            3) For 'os' and 'wallet' authentications, username and password are automatically set to ''.

        Parameters
        ----------
        db_dict: dict
            Dictionary of database connection information following the above constraints and requirements.

        Returns
        -------
        session : SQLAlchemy session
            A bound SQLAlchemy session instance that can be passed to pisces for queries

    """

    # assign dialect, authenticaion, and driver at start since they are needed for all steps of
    # decision tree
    dialect = db_info['dialect']
    authentication = db_info['authentication']
    
    if 'driver' not in db_info.keys():
        driver = ''
        driverStr= driver
    else:
        driver = db_info['driver']
        driverStr = '+{}'.format(driver)

    # is dialect oracle? Does thick mode need to be enabled?
    # this is true for os authentication and wallets
    if dialect == 'oracle' and authentication in ['os', 'wallet'] and driver == 'oracledb':
        ex_auth = True
    else:
        ex_auth = False

    # find or assemble url
    if 'url' not in db_info.keys() or db_info['url'] == '':
        if authentication in ['os', 'wallet']:
            username = ''
            password = ''
        else:
            username = db_info['username']
        
        hostname = db_info['hostname']
        db_name = db_info['database_name']
        port = db_info['port']

        if 'password' not in db_info.keys():
            password = ''
        else:
            password = db_info['password']

           

        if 'alias' not in db_info.keys():
            alias = None
        else:
            alias = db_info['alias']

        if alias is None or alias == '':
            dburl = '{}{}://{}:{}@{}:{}/{}'.format(dialect, driverStr, username, password, hostname, port, db_name)

        else:
            dburl = '{}{}://{}:{}@{}'.format(dialect, driverStr, username, password, alias)
        
    else:
        dburl = db_info['url']

    print("  Connecting to database through url: " + dburl)

    session = db_connect_url(dburl, ex_auth = ex_auth)

    return session

def db_connect_url(url, ex_auth = False):
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
    if ex_auth == True:
        engine = create_engine(url, thick_mode = True, connect_args = {"externalauth": True} )
        session = Session(bind=engine)
    else:
        engine = create_engine(url)
        session = Session(bind=engine)
    
    return session

def set_db_env_variables(env_vars):
    # env_vars is a dictionary containing the environment variable as the key, and what to set it to as the value
    # Note that the environment variables will only last for the duration of the session.
    for key, value in env_vars.items():
            os.environ[key] = value

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
        origins = ps.request.get_events(session, db_tables['origin'], event=db_tables['event'], evids=evIDs)
        prefor = ps.request.get_events(session, db_tables['origin'], event=db_tables['event'], evids=evIDs, prefor=True)

    if origins:
        return prefor, origins
    else:
        print("no event found")
        return None


def ev_query_area(session, center_lat, center_lon, minr, maxr, db_tables):

    if 'Origin' not in db_tables or 'Event' not in db_tables:
        raise KeyError

    events = ps.request.get_events(session, db_tables['origin'], db_tables['events'], km=(center_lat, center_lon, minr, maxr), etime=(startt, endt))



#####################################
## New stuff with Christine's help ##
#####################################
def set_session(db_info):
    """
        Function to set up a database connection session using a dictionary of connection information.
        The current configuration file format is below ([DATABASE] block passed to this function):

            [DATABASE]
            hostname = [...]
            username = [...]
            database_name = [...]
            port = [...]
            driver = [...]
            dialect = [...]
            url = [...]

            [DBTABLES]
            schema = KBCore
            wfdisc = global.wfdisc_raw
            site = global.site
            origin = global.origin
            event = global.event

        Parameters
        ----------
        db_info: dict
            Dictionary of database connection information

        Returns
        -------
        session : SQLAlchemy session
            A bound SQLAlchemy session instance that can be passed to pisces for queries

    """

    print("Setting up database configuration...")
    for key in db_info.keys():
        print('  ' + key + '\t' + db_info[key])

    session = db_connect(db_info)

    # check the session works
    try:
        session.get_bind().connect()
    except Exception as e:
        print("Database connection failed")
        return

    return session


# Methods to write arrivals into JSON files for use in InfraPy
def _query2det_dict(query_line, array_dim):

    wvfrm_info = [{"trace id" : "__." + query_line.arrival.sta + ".." + query_line.arrival.chan,
                   "latitude" : query_line.site.lat,
                   "longitude" : query_line.site.lon}]

    if "I" in query_line.arrival.iphase:
        # infrasound detection entry

        # if direciton-of-arrival is defined, convert slowness s/deg to s/km and invert to get trace velocity
        if query_line.arrival.azimuth > 0.0:
            back_az = query_line.arrival.azimuth
            tr_vel = 111.32 / query_line.arrival.slow
        else:
            back_az = None
            tr_vel = None

        # check if a valid snr value is defined and convert to f-stat
        if query_line.arrival.snr > 0.0:
            f_stat = query_line.arrival.snr * array_dim + 1.0
        else:
            f_stat = None

        det = {'wvfrm_info' : wvfrm_info,
               'peak f-stat time': UTCDateTime(query_line.arrival.time),
               'f-stat': f_stat,
               'tr vel': tr_vel,
               'back az': back_az,
               'phase id' : query_line.arrival.iphase}
    else:
        # seismic pick entry
        det = {'wvfrm_info' : wvfrm_info,
               'pick time' : UTCDateTime(query_line.arrival.time),
               'deltim' : query_line.arrival.deltim,
               'slow' : query_line.arrival.slow,
               'azimuth' : query_line.arrival.azimuth,
               'snr' : query_line.arrival.snr,
               'phase id' : query_line.arrival.iphase,
               }

    return det


def db2dets(db_session, db_tbls, lat_bnds, lon_bnds, starttime, endtime, phase_list="I"):
    """
        Function to pull all detections within lat/lon bounds and between start and end times into a dictionary form usable in InfraPy (writable into event JSON)
        Uses the same config file structure as above, and expects the [DBTABLES] block.

            [DATABASE]
            hostname = [...]
            username = [...]
            database_name = [...]
            port = [...]
            driver = [...]
            dialect = [...]
            url = [...]

            [DBTABLES]
            schema = KBCore
            wfdisc = global.wfdisc_raw
            site = global.site
            origin = global.origin
            event = global.event

        Parameters
        ----------
        session : SQLAlchemy session
            A bound SQLAlchemy session instance that can be passed to pisces for queries
        db_tbls: dict
            Dictionary of database connection information
        lat_bnds: iterable float
            Minimum and maximum latitude bounds for search
        lon_bnds: iterable float
            Minimum and maximum longitude bounds for search
        starttime: obspy.UTCDatetime
            Minimum time for search
        endtime: obspy.UTCDatetime
            Maximum time for search
        phase_list: string
           Comma separated string of phases to include (e.g., 'I,P,S')

        Returns
        -------
        session : SQLAlchemy session
            A bound SQLAlchemy session instance that can be passed to pisces for queries

    """

    array_dim = 5

    # define tables needed for query
    if 'css3' in db_tbls['schema'].lower():
        db_schema = ps.schema.css3
    elif 'kbcore' in db_tbls['schema'].lower():
        db_schema = ps.schema.kbcore
    else:
        print("Invalid database schema.")
        return None

    # need a catch if these aren't in the user config file...
    class arrival(db_schema.Arrival):
        __tablename__ = db_tbls['arrival']

    class sitechan(db_schema.Sitechan):
        __tablename__ = db_tbls['sitechan']

    class site(db_schema.Site):
        __tablename__ = db_tbls['site']

    class sensor(db_schema.Sensor):
        __tablename__ = db_tbls['sensor']

    arrQuery = db_session.query(arrival, site, sitechan, sensor) \
                .filter(arrival.time.between(starttime, endtime)) \
                .filter(site.lat.between(lat_bnds[0], lat_bnds[1])) \
                .filter(site.lon.between(lon_bnds[0], lon_bnds[1])) \
                .filter(arrival.chanid == sensor.chanid) \
                .filter(arrival.time.between(sensor.time, sensor.endtime)) \
                .filter(sensor.chanid == sitechan.chanid) \
                .filter(site.sta == sitechan.sta) \
                .filter(sitechan.ondate.between(site.ondate, site.offdate))

    # Extract pick information
    det_dicts = []
    for line in arrQuery:
        if line.arrival.iphase in phase_list.replace(" ","").split(","):
            det_dicts = det_dicts + [_query2det_dict(line, array_dim)]

    return det_dicts


def db2ev(db_session, db_tbls, evid, phase_list="I"):
    """
        Function to pull all detections associated with an event and return in a dictionary form usable in InfraPy (writable into event JSON)
        Uses the same config file structure as above, and expects the [DBTABLES] block.

            [DATABASE]
            hostname = [...]
            username = [...]
            database_name = [...]
            port = [...]
            driver = [...]
            dialect = [...]
            url = [...]

            [DBTABLES]
            schema = KBCore
            wfdisc = global.wfdisc_raw
            site = global.site
            origin = global.origin
            event = global.event

        Parameters
        ----------
        session : SQLAlchemy session
            A bound SQLAlchemy session instance that can be passed to pisces for queries
        db_tbls: dict
            Dictionary of database connection information
        evid: int
            Event identifier in the database for event of interest
        phase_list: string
           Comma separated string of phases to include (e.g., 'I,P,S')

        Returns
        -------
        session : SQLAlchemy session
            A bound SQLAlchemy session instance that can be passed to pisces for queries

    """

    array_dim = 5

    # define tables needed for query
    if 'css3' in db_tbls['schema'].lower():
        db_schema = ps.schema.css3
    elif 'kbcore' in db_tbls['schema'].lower():
        db_schema = ps.schema.kbcore
    else:
        print("Invalid database schema.")
        return None

    # need a catch if these aren't in the user config file...
    class event(db_schema.Event):
        __tablename__ = db_tbls['event']

    class origin(db_schema.Origin):
        __tablename__ = db_tbls['origin']

    class assoc(db_schema.Assoc):
        __tablename__ = db_tbls['assoc']

    class arrival(db_schema.Arrival):
        __tablename__ = db_tbls['arrival']

    class sitechan(db_schema.Sitechan):
        __tablename__ = db_tbls['sitechan']

    class site(db_schema.Site):
        __tablename__ = db_tbls['site']

    class sensor(db_schema.Sensor):
        __tablename__ = db_tbls['sensor']

    # build query
    arrQuery = db_session.query(event, origin, assoc, arrival, site, sitechan, sensor) \
                .filter(event.evid == evid) \
                .filter(event.prefor == origin.orid) \
                .filter(origin.orid == assoc.orid) \
                .filter(arrival.arid == assoc.arid) \
                .filter(arrival.chanid == sensor.chanid) \
                .filter(arrival.time.between(sensor.time, sensor.endtime)) \
                .filter(sensor.chanid == sitechan.chanid) \
                .filter(site.sta == sitechan.sta) \
                .filter(sitechan.ondate.between(site.ondate, site.offdate))

    ev_output = {'ground truth' : {},
                 'det_info' : [],
                 'location' : [],
                 'characterization' : []}

    # Extract GT info from prefor origin info
    ev_output['ground truth']['latitude'] = arrQuery[0].origin.lat
    ev_output['ground truth']['longitude'] = arrQuery[0].origin.lon
    ev_output['ground truth']['origin time'] = UTCDateTime(arrQuery[0].origin.time)
    ev_output['ground truth']['name'] = arrQuery[0].event.evname

    # Extract pick information
    for line in arrQuery:
        if line.arrival.iphase in phase_list.replace(" ","").split(","):
            ev_output['det_info'] = ev_output['det_info'] + [_query2det_dict(line, array_dim)]

    return ev_output


