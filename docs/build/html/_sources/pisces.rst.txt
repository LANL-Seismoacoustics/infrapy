.. _pisces:

=====================================
Interfacing with Pisces
=====================================

Infrapy leverages pisces to connect with and process data in either local sqlite databases or oracle databases. More information about pisces can be found at https://jkmacc-lanl.github.io/pisces/.

-------------------------------------
Converting Data into Sqlite Databases
-------------------------------------
Data in miniseed or sac formats can be loaded into a sqlite database for pipeline processing using commands from pisces.

1. mseed to database (ms2db.py)

.. code-block:: python

    >> ms2db.py sqlite:///example.sqlite mslist.txt

2. sac to database (sac2db)


.. code-block:: python

    >> pisces sac2db sqlite:///example.sqlite *.sac

As infrapy is an array processing tool, after your sqlite database is created, you will need to update the REFSTA for each array using update_refsta.py

.. code-block:: python

    >> update_refsta.py sqlite:///example.sqlite <array name>

You can update the calibration for each array using update_calib.py

.. code-block:: python

    >> update_calib.py sqlite:///example.sqlite <array name> <calibration>

-------------------------------------
Connecting to a SQL Database
-------------------------------------

Infrapy employs two main methods for connecting to either Oracle or sqlite databases.  Example files to facilitate these connections are found in tutorials/.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Defining Schema Specific Tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pipeline processing in infrapy utilizes information from CSS3.0 Site and Wfdisc tables.  If your database schema differs from the CSS3.0 schema in any way, you can define the differences using a _global.py file.  An example _global.py file is found in tutorial/ .

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Connection within pipeline processing configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first three lines of your configuration file define the database you will connect to:

**Example Configuration File for Sqlite Processing (Sqlite_Config.txt)**

.. code-block:: none

    [database] # required
    # url to database where you have the pointers to data and metadata
    url = sqlite:///example.sqlite
    # schema specific tables for your site and wfdisc files.  If you are processing in a sqlite database, these variables will refer to schema specified in pisces. If you are processing in an oracle database, these variables will refer to schema specified in your global_.py file
    site = pisces.tables.css3:Site
    wfdisc = pisces.tables.css3:Wfdisc


**Example Configuration File for Oracle DB Processing (Oracle_Config.txt)**

.. code-block:: none

    [database] # required
    url = oracle://<database name>:<port>
    site = global_:Site
    wfdisc = global_:Wfdisc_raw

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Connection with a db.cfg file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some modules in infrapy (db2sac) require a .cfg file to establish connection with a database.  Examples are found in tutorial/ . More information can be found in the pisces documentation.

**Example Configuration File for Oracle DB Processing (oracle_connection.cfg)**

.. code-block:: none

    [database] # required
    url = oracle://<db name>:<db port>
    site = global_:Site
    wfdisc = global_:Wfdisc_raw
    origin = global_:Origin

**Example Configuration File for Sqlite Processing (sqlite_connection.cfg)**

.. code-block:: none

    [database] # required
    url = sqlite:///example.sqlite
    site = pisces.tables.css3:Site
    wfdisc = pisces.tables.css3:Wfdisc
    origin = pisces.tables.css3:Origin
