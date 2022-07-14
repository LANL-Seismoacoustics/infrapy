.. _utilities:

=====================================
Utility Functions
=====================================

InfraPy includes an interactive graphical user interface (GUI) called InfraView.  This interface includes waveform analysis methods, beamforming and detection tools, and event identification and localization capabilities.  The InfraView GUI can be accessed from the command line by simply running,

    .. code-block:: bash

        >> infrapy utils --help

    .. code-block:: none

        Usage: infrapy utils [OPTIONS] COMMAND [ARGS]...

          infrapy utils - various utility functions for infrapy usage

        Options:
          -h, --help  Show this message and exit.

        Commands:
          arrival-time     Estimate the arrival time for a source-receiver pair
          arrivals2json    Convert infraGA/GeoAc arrivals to detection file
          best-beam        Compute the best beam via shift/stack
          calc-celerity    Compute the celerity for an arrival from a known source
          check-db-wvfrms  Check waveform pull from database
          write-wvfrms     Save waveforms from FDSN or database

- Estimating Arrival Time Bounds

    - Stuff...

        .. code-block:: bash

            infrapy utils arrival-time --help

        .. code-block:: none

            Usage: infrapy utils arrival-time [OPTIONS]

              Compute the range of possible arrivals times for a source-receiver pair
              given a range of celerity values. Can use a receiver latitude/longitude or
              reference from a list (currently only IMS stations)

              Example usage (requires InfraGA/GeoAc arrival output):
                  infrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr-lat 40.0 --rcvr-lon -110.0
                  infrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr I57US

            Options:
              --src-lat TEXT        Source latitude
              --src-lon TEXT        Source longitude
              --src-time TEXT       Source time
              --rcvr-lat TEXT       Receiver latitude
              --rcvr-lon TEXT       Receiver longitude
              --rcvr TEXT           Reference IMS station (e.g., 'I53')
              --celerity-min FLOAT  Minimum celerity
              --celerity-max FLOAT  Maximum celerity
              -h, --help            Show this message and exit.

    - Running the first example...

        .. code-block:: bash

            infrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr-lat 40.0 --rcvr-lon -110.0

        .. code-block:: none

            #################################
            ##                             ##
            ##      InfraPy Utilities      ##
            ##         arrival-time        ##
            ##                             ##
            #################################

              Source Time: 2020-12-25T00:00:00
              Source Location: (30.0, -110.0)
              Receiver Location: (40.0, -110.0)

              Celerity Range: (0.24, 0.35)

              Propagation range: 1111.95 km
              Propagation azimuth: 0.0 degrees
  
              Estimated arrival back azimuth: 180.0 degrees
              Estimated arrival time range:
                2020-12-25T00:52:57
                2020-12-25T01:17:13
                

    - Running with a reference IMS station...

        .. code-block:: bash

            infrapy utils arrival-time --src-lat 30.0 --src-lon -110.0 --src-time "2020-12-25T00:00:00" --rcvr I57

        .. code-block:: none

            #################################
            ##                             ##
            ##      InfraPy Utilities      ##
            ##         arrival-time        ##
            ##                             ##
            #################################

              Source Time: 2020-12-25T00:00:00
              Source Location: (30.0, -110.0)

              User specified reference receiver: I57
              Reference IMS station match: I57US
              Receiver Location: (33.6064, -116.455)

              Celerity Range: (0.24, 0.35)

              Propagation range: 729.75 km
              Propagation azimuth: -55.01 degrees
  
              Estimated arrival back azimuth: 121.59 degrees
              Estimated arrival time range:
                2020-12-25T00:34:45
                2020-12-25T00:50:41


- Computing the best beam waveform from fk analysis results...
  
    - Stuff...

        .. code-block:: bash

            infrapy utils best-beam --help

    - More stuff...

        .. code-block:: bash

            infrapy utils best-beam --config-file config/detection_local.config


    - More stuff...

        .. code-block:: bash

            infrapy utils best-beam --config-file config/detection_local.config --back-az -39.0 --trace-vel 358.0


    - More stuff...

        .. code-block:: bash

            infrapy utils best-beam --config-file config/detection_local.config --signal-start '2012-04-09T18:13:00' --signal-end '2012-04-09T18:15:00'


- Calculating the arrival celerity from an arrival with a known source
  
    - Stuff...

        .. code-block:: bash

            infrapy utils calc-celerity --help

    - More stuff...


- Write waveforms from an FDSN or database source
    
    - Stuff...

        .. code-block:: bash

            infrapy utils write-wvfrms --help

    - More stuff...

        .. code-block:: bash

            infrapy utils write-wvfrms --config-file config/detection_fdsn.config

    - More stuff...
