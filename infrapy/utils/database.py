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

##################################
## NEW STUFF PHIL IS WORKING ON ##
##################################
def set_session(db_info, db_schema='css3'):

    print("Setting up database configuration...")
    # set up the session and check connection
    if 'url' in db_info['DATABASE'].keys():
        print("  Connecting to database through url: " + db_info['DATABASE']['url'])
        db_session = ps.db_connect( db_info['DATABASE']['url'])
    else:
        # clean up the above to simplify this or just require a url?
        db_session = db_connect2(db_info)

    # check the session works
    try:
        db_session.get_bind().connect()
    except Exception as e:
        print("Database connection failed")
        return 

    db_tables = {}

    if 'css3' in db_schema:
        core_tbls = ps.tables.css3.CORETABLES
        print(core_tbls)

        # core_tbls["sensor"] = ps.tables.css3.Sensor

    elif 'kbcore' in db_schema:
        core_tbls = ps.tables.kbcore.CORETABLES
        print(core_tbls)

        # core_tbls["sensor"] = ps.tables.kbcore.Sensor
    else:
        print("Unrecognized database schema! Options are 'css3' and 'kbcore'")

    for tbl_name in core_tbls:
        print(tbl_name)
        if tbl_name in db_info['DBTABLES']:
            print("  Setting " + tbl_name + " table from: " + db_info['DBTABLES'][tbl_name])
            class tbl(core_tbls[tbl_name].prototype):
                __tablename__ = db_info['DBTABLES'][tbl_name]
            db_tables[tbl_name] = tbl
    
    if 'sensor' in db_info['DBTABLES']:
        print("  Setting sensor table from: " + db_info['DBTABLES']['sensor'])
        if 'css3' in db_schema:
            class tbl(ps.schema.css3.Sensor):
                __tablename__ = db_info['DBTABLES']['sensor']
        else: 
            class tbl(ps.schema.kbcore.Sensor):
                __tablename__ = db_info['DBTABLES']['sensor']
        db_tables['sensor'] = tbl

    return db_session, db_tables



# Methods to write arrivals into JSON files for use in InfraPy
def query2det_dict(query_line, array_dim):

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
    

def db2dets_json(db_session, db_tbls, lat_bnds, lon_bnds, starttime, endtime, phase_list="I", db_schema="css3"):

    array_dim = 5 
    
    # define tables needed for query
    if 'css3' in db_schema.lower():
        db_schema = ps.schema.css3
    elif 'kbcore' in db_schema.lower(): 
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
            det_dicts = det_dicts + [query2det_dict(line, array_dim)]

    return det_dicts


def db2ev_json(db_session, db_tbls, evid, phase_list="I", db_schema="css3"):

    array_dim = 5

    # define tables needed for query
    if 'css3' in db_schema.lower():
        db_schema = ps.schema.css3
    elif 'kbcore' in db_schema.lower(): 
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
            ev_output['det_info'] = ev_output['det_info'] + [query2det_dict(line, array_dim)]

    return ev_output


