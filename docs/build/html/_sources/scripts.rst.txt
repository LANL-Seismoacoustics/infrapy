.. _scripts:

=====================================
Scripts
=====================================

Scripts to manipulate Infrapy FK results
________________________________________

Use "any_command".py -h to get specific information to run the script

1. read_pfk.py to see the different set of configuration parameters used previously for FK analysis

| -h    --help     show this help message and exit
| -d    SQ      name of the database connection, e.g.: -d sqlite:///UT_tutorial.sqlite


.. code-block:: python

    >> read_pfk.py [-h] -d SQ


2. print_rfk.py to see fk results for a specific array and FK parameter ID

| -h --help     show this help message and exit
| -d SQ         name of the database connection, e.g.: -d sqlite:///UT_tutorial.sqlite
| -a ARRAY      array name, e.g.: -HWU4
| -t FKRESULTS  specific table with results, e.g.: -t fk_HWU
| -s TS        starttime plot, e.g.: -s /'2014-03-02T00:00:00/'
| -e TE        endtime plot, e.g.: -s /'2014-03-03T00:00:00/'
| -F FVAL      limit Fvalue, e.g.: -F 0
| -o OUTTEXT    print fk data, e.g.: -o res_FILE

.. code-block:: python

    >> print_rfk.py [-h] -d SQ -a ARRAY -f PFK_ID -t FKRESULTS [-s TS] [-s TE] [-F FVAL] -o res_FILE

3. plot1_rfk.py to plot FK results

| -h --help    show this help message and exit
| -d SQ        name of the database connection, e.g.: -d sqlite:///UT_tutorial.sqlite
| -a ARRAY     array name, e.g.: -a HWU4
| -f PFK_ID    FK parameter ID to be plot, e.g.: -f 3
| -t FKRESULTS specific table with results, e.g.: -t fk_res_HWU
| -w WAVEPLOT  plot waveforms, e.g.: -w 0
| -s TS        starttime plot, e.g.: -s /'2014-03-02T00:00:00/'
| -e TE        endtime plot, e.g.: -s /'2014-03-03T00:00:00/'
| -F FVAL      limit Fvalue, e.g.: -F 0
| -slo SLOFK   limit slofk, e.g.: -slo 0
| -bzmin BZMIN limit min bz, e.g.: -bzmin 0
| -bzmax BZMAX limit max bz, e.g.: -bzmax 360

.. code-block:: python

    >>plot1_rfk.py [-h] -d SQ -a ARRAY -f PFK_ID [-t FKRESULTS] [-w WAVEPLOT] [-s TS] [-e TE] [-F FVAL] [-slo SLOFK] [-bzmin BZMIN] [-bzmax BZMAX]


.. image:: _static/_images/FKres_Fval.png

.. image:: _static/_images/FKres_trcvel.png

.. image:: _static/_images/FKres_bz.png



Scripts to manipulate Infrapy FD results
________________________________________


1. read_pfd.py to see the different set of configuration parameters used previously for detection analysis

| -h, --help  show this help message and exit
| -d SQ       name of the database connection, e.g.: -d sqlite:///mydb.sqlite

.. code-block:: python

    >> read_pfd.py [-h] -d SQ


2. read_rfd.py to see the available detection results

| -h, --help            show this help message and exit
| -d SQ                 name of the database connection, e.g.: -db sqlite:///UT_tutorial.sqlite
| -a ARRAY              array name, e.g.: -a HWU4
| -f PFK_ID, --pfkid PFK_ID FK parameter ID to be plot, e.g.: -f 0
| -j PFDID, --pfdid PFDID fd parameter id, e.g.: -j 0
| -t FDRESULTS          specific table with results, e.g.: -t fd_res_HWU

.. code-block:: python

    >> read_rfd.py [-h] -d SQ -a ARRAY -f PFK_ID -j PFDID [-t FDRESULTS]
    fdid: 1  pfdid: 0 pfkid: 0   timeini: 12-08-14 00:24:30   timeend: 12-08-14 00:26:00   maxf0: 5.54282460217 0.60800443287
    fdid: 2  pfdid: 0 pfkid: 0   timeini: 12-08-14 00:47:00   timeend: 12-08-14 00:53:30   maxf0: 3.67253815208 0.46716764663
    fdid: 3  pfdid: 0 pfkid: 0   timeini: 12-08-14 00:54:30   timeend: 12-08-14 00:56:30   maxf0: 2.61098286001 0.372113714177
    fdid: 4  pfdid: 0 pfkid: 0   timeini: 12-08-14 01:47:30   timeend: 12-08-14 01:49:00   maxf0: 9.91099311406 0.749628285504


3. read_rfd_fast.py to write the available detection results in text file

| -h, --help            show this help message and exit
| -d SQ                 name of the database connection, e.g.: -d sqlite:///UT_tutorial.sqlite
| -a ARRAY              array name, e.g.: -a I37NO
| -f PFKID, --pfkid PFKID fk parameter id, e.g.: -f 0
| -j PFDID, --pfdid PFDID fd parameter id, e.g.: -j 0
| -t FDRESULTS          specific table with results, e.g.: -t fd_I37
| -o OUTTEXT            fd parameter id, e.g.: -o res_FILE

.. code-block:: python

    >> read_rfd_fast.py [-h] -d SQ -a ARRAY -f PFKID -j PFDID [-t FDRESULTS] [-o OUTTEXT]
