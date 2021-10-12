# put database related functions here
# generic db functions ONLY.  No LANL specific code at all

import sqlalchemy as sa
import pisces as ps




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
    
    # in case someone passed the port as an integer, first convert it to a string
    str_port = str(port)

    if driver:
        driver = '+' + driver
    
    connect_str = dialect + driver + "://" + username + ":" + password + "@" + hostname + ":" + str_port + "/" + db_name

    session = ps.db_connect(connect_str)

    print(connect_str)
    return session