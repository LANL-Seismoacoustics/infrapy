"""
  put database related functions here
  generic db functions ONLY.  No LANL specific code at all
"""

import pisces as ps
import pisces.schema.kbcore as kb
import pandas as pd

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
    except:
        return False


class Site(kb.Site):
    __tablename__ = 'global.site'

class Wfdisc(kb.Wfdisc):
    __tablename__ = 'global.wfdisc_raw'

def query_db(session, start_time, end_time, sta="*", loc="*", cha="*"):

    if session is None:
        return None

    my_query = session.query(Wfdisc).filter(Wfdisc.sta == sta)\
                                    .filter(Wfdisc.time < end_time.timestamp)\
                                    .filter(Wfdisc.endtime > start_time.timestamp)\
                                    .filter(Wfdisc.chan.like(cha))

    df = pd.read_sql(my_query.statement, session.bind)
    
    return df
