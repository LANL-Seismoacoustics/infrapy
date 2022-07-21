.. _pisces:

====================
Database Interfacing
====================


Infrapy leverages pisces to connect with and process data in databases. More information about pisces can be found at https://jkmacc-lanl.github.io/pisces/.  The methods described here are still in active development; therefore, if something doesn't work please contact both the infrapy POCs as well as the pisces POCs for assistance and clarification.

---------------------------------
Defining Database Connection Info
---------------------------------

- In order to connect and pull data from a database, the pisces methods need to establish a database session and perform a number of table queries.  The information needed to connect and perform these queries is defined within a separate configuation file specifically for the database that can be defined either on the command line using :code:`--db-config` or in a general infrapy configuration file as below:

    .. code-block:: none

        [WAVEFORM IO]
        db_config = /path/to/database_info.config

        station = I53*
        channel = *DF
        starttime = 2018-12-19T01:00:00
        endtime = 2018-12-19T03:00:00

        ...

- The database configuration file includes information about the database itself (e.g., schema, dialect, hostname) as well as a list of database table names that will be queried and used in analysis.  An example database configuration file contents is below:

    .. code-block:: none

        [DATABASE]
        schema = KBCore
        dialect = oracle
        hostname = my_host.com
        database_name = my_db
        port = 1234
        driver = 
        username = 
        password =  
        owner = my_account

        [DBTABLES]
        wfdisc = wfdisc
        site = site

    - Note that the driver, username, and password fields are optional and if accessing your database doesn't require them they can be left out of the config file or simply left blank as in the above

    - The owner information is used in combination with the database table list to define the actual table names (e.g., the above would define the site table reference as :code:`my_account.site`).

    - Add notes about other fields?

-----------------------------------------
Accessing Waveform Data from the Database
-----------------------------------------

- The utility function, :code:`infrapy utils check-db-wvfrms`, can be used to check what waveform data will be extracted to ensure that your database connection is working and that data is being pulled correctly.  An example output of this function is shown below:

    .. code-block:: none


        #################################
        ##                             ##
        ##      InfraPy Utilities      ##
        ##       check_db_wvfrms       ##
        ##                             ##
        #################################


        Loading configuration info from: my_configuration.config

        Data parameters:
          db_config: /path/to/database_info.config
          network: None
          station: I53*
          location: None
          channel: *DF
          starttime: 2018-12-19T01:00:00
          endtime: 2018-12-19T03:00:00

        Loading data from database...

        Data summary:
        __.I53H1..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        __.I53H2..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        __.I53H3..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        __.I53H4..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        __.I53H5..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        __.I53H6..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        __.I53H7..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        __.I53H8..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z

        Location info:
        64.87500	-147.86114
        64.87236	-147.83828
        64.86208	-147.84281
        64.85914	-147.86615
        64.86805	-147.87858
        64.86722	-147.86120
        64.86617	-147.85670
        64.86541	-147.85990
        
    - Note that all 8 channels of the I53 station are found in the database and the displayed start and end times indicate that the requested 2 hours of data are able to be pulled.  The latitude and longitude of the various station elements are summarized for completeness.

    - Also, pulling waveform data from a database requires only the :code:`site` table to identify station information and :code:`wfdisc` to access waveform data.  Additional tables will be required for attaching response information (currently data is returned in counts and not physical units) and other advances in accessing waveform data.

- Running and visualizing beamforming analysis with this data can now be done using the same :code:`my_configuration.config` file with parameter modifications and other details as discussed in the :ref:`quickstart`.  Essentially, when pulling waveform data for analysis, our aim is to be able to point the infrapy methods to a database via pisces and a :code:`database_info.config` file as easily as pointing it at an FDSN.

------------------------------------------
Writing Analysis Results into the Database
------------------------------------------

- Currently, analysis results using waveform data pulled from a database are written into [...].fk_results.dat and [...].dets.json files as discussed in the :ref:`quickstart`.  In future developments, methods to write those results back into database tables will be implemented.  
  
- The current plan is to implement such functionality as part of the :ref:`utilities` methods and include methods to both write detection results into the appropriate tables as well as pull detections from a database for event ID and similar analysis.  


