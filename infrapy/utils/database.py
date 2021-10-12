# put database related functions here
# generic db functions ONLY.  No LANL specific code at all

import sqlalchemy as sa
import pisces as ps
import pisces.schema.kbcore as kb

import urllib.parse

from obspy.core.utcdatetime import UTCDateTime


valid_dialects = ["postgresql", "mysql", "oracle", "mssql", "sqlite"]


def db_connect(dialect, hostname, db_name, port=None, username="", password="", driver=""):
    """
    Connect to a database to do database things...

    Parameters
    ----------
    dialect: str
        Type of database.   (examples: oracle, mysql)
    hostname : str
        The url of the database. (example: mydb.home.org)
    db_name: str
        Name of your database
    port: int or string
        Port number of your database
    username: str
        login name for your database connection.  For self authenticated sites, leave this as an empty string
    password: str
        used in conjunction with username.  For self authenticated sites, leave this as an empty string
    driver: str
        driver api used for connection/  (example: pymysql) 

    Returns
    -------
    session : bound SQLAlchemy session instance

    """
    # before we do anything, lets clean up some of the input...

    # in case someone passed the port as an integer, first convert it to a string
    str_port = str(port)

    # deal with special characters in the username and password
    if password:
        password = urllib.parse.quote_plus(password)

    #make sure things that need to be lower case are lower case
    dialect = dialect.lower()

    if driver:
        driver = driver.lower()
        driver = '+' + driver
    
    connect_str = dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + str_port + "/" + db_name

    return ps.db_connect(connect_str)



if __name__ == '__main__':
    dialect = 'oracle'
    hostname = 'gnemdb.lanl.gov'
    db_name = 'gnem19c'
    port = 1523

    session = db_connect(dialect, hostname, db_name, port=port)
    print(session)
    print(session.get_bind())

    class Site(kb.Site):
        __tablename__ = 'global.site'

    class Wfdisc(kb.Wfdisc):
        __tablename__ = 'global.wfdisc_raw'

    staname = 'ELK'
    arrname = 'NVAR'
    netname = 'USGS'

    origintime = UTCDateTime(2019,7,6,3,19,53)
    starttime = origintime - (5*60)
    endtime = origintime + (10*60)

    q_station=session.query(Wfdisc).filter(Wfdisc.sta == staname).\
        filter(Wfdisc.time<endtime.timestamp).\
        filter(Wfdisc.endtime>starttime.timestamp).\
        filter(Wfdisc.chan.like('%H%'))