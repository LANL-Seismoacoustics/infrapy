.. _quickstart:

=====================================
Quickstart
=====================================

****************************
Command Line Interface (CLI) 
****************************

Most of InfraPy's analysis methods are accessible through a command line interface (CLI) with configuation options either via command line parameter definitions or a configuation file.  Waveform data can be ingested from local files (eg., SAC or similar format that can be ingested via :code:`obspy.core.read`) or downloaded from FDSN clients via :code:`obspy.clients.fdsn`.  Array- and network-level analyses can be performed from the command line enabling a full pipeline of analysis from beamforming/detection to event identification and localization.  Visualization methods are also included to quickly interrogate analysis results.  The Quickstart summarized here steps through these various CLI methods and demonstrates the usage of InfraPy from the command line.

--------------------
Array-Level Analyses
--------------------

- The beamforming methods in InfraPy can be run via the :code:`run_fk` CLI option.  For a local data source such as the included SAC files in the data directory, this is simply,

    .. code-block:: bash

        infrapy run_fk --local-wvfrms 'data/YJ.BRP*.SAC'

    Note that the data path must be in quotes in order to properly parsed.  As the methods are run, data and algorithm parameters are summarized and a progress bar shows how much of the data has been analyzed:

    .. code-block:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##    Beamforming (fk) Analysis    ##
        ##                                 ##
        #####################################

        Data parameters:
          local_wvfrms: data/YJ.BRP*.SAC
          local_latlon: None
          local_fk_label: None

        Algorithm parameters:
          freq_min: 0.5
          freq_max: 5.0
          back_az_min: -180.0
          back_az_max: 180.0
          back_az_step: 2.0
          trace_vel_min: 300.0
          trace_vel_max: 600.0
          trace_vel_step: 2.5
          method: bartlett
          signal_start: None
          signal_end: None
          window_len: 10.0
          sub_window_len: None
          window_step: 5.0
          multithread: False

        Loading local data from data/YJ.BRP*.SAC

        Data summary:
        YJ.BRP1..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP2..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP3..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP4..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z

        Running fk analysis...
	        Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]

        Writing results to specified output...

- Once completed, this analysis produces three output files: 

+---------------------------------------------------+-----------------------------------------------+
| YJ.BRP4_2012.04.09_18.00.00-18.19.59.fk_meta.txt  | Meta-data summarizing the fk calculations     |
+---------------------------------------------------+-----------------------------------------------+
| YJ.BRP4_2012.04.09_18.00.00-18.19.59.fk_times.npy | Analysis window times                         |
+---------------------------------------------------+-----------------------------------------------+
| YJ.BRP4_2012.04.09_18.00.00-18.19.59.fk_peaks.npy | Peak info (Fisher stat, direction of arrival) |
+---------------------------------------------------+-----------------------------------------------+

- The beamforming results from the :code:`infrapy run_fk` analysis can be visualized using:

    .. code-block:: bash

        infrapy plot_fk --local-wvfrms 'data/YJ.BRP*.SAC'

    The resulting plot of the included example data set is shown below for comparison:

    .. image:: _static/_images/plot_fk.png
        :width: 1200px
        :align: center

- The default beamforming parameters in :code:`run_fk` are useful, but in many cases the frequency band for a signal of interst or the window length appropriate for a given frequency band needs to be modified.  From the command line, this can be done by specifying a number of options in the algorthm as summarized in the :code:`--help` information.  For example, the analysis of data from BRP can be completed using a modified frequency band via:

    .. code-block:: bash

        infrapy run_fk --local-wvfrms 'data/YJ.BRP*.SAC' --freq-min 1.0 --freq-max 8.0

    The fk output files are automatically named from the data file (network and station codes plus start and end times), but a label can be specified as :code:`--local_fk-label example`.

- In the case that multiple analysis parameters are changed from their default values, a configuration file is useful to simplify running analysis and keep a record of what was used for future review of anlaysis.  Create a text file called :code:`BRP_analysis.config` and enter the following:

    .. code-block:: none

        [WAVEFORM IO]
        local_wvfrms = data/YJ.BRP*.SAC

        [DETECTION IO]
        local_fk_label = BRP_analysis

        [FK]
        freq_min = 1.0
        freq_max = 8.0
        window_len = 5.0
        window_step = 2.5
        cpu_cnt = 8

    Adjust the CPU count value to whatever number of available threads you have on your machine.  The analysis can now be completed by simply running:

    .. code-block:: bash

        infrapy run_fk --config-file BRP_analysis.config

    When using a config file for analysis, any additional parameters set on the command line will overwrite the values from the config file.  For example, to run the analysis with a maximum frequency of 10 Hz instead of 8 Hz, one can simply run:

    .. code-block:: bash

        infrapy run_fk --config-file BRP_analysis.config --freq-max 10.0

    If a parameter is not included in a config file or via the command line, a default value is used and can be found in the ouput at the time of the analysis or in the meta-data file.

- From the beamforming results, detection analysis can be conducted via the :code:`run_fd` method.  This anlaysis requires the fk output label and can use a custom detection label or automatically use the fk label if none is specified.

    .. code-block:: bash

        infrapy run_fd --local-fk-label data/BRP_analysis


    Similarly to the :code:`run_fk` methods, parameter summaries are provided; however, because this anlaysis is relatively quick there is no progress bar:

    .. code-block:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##     Detection (fd) Analysis     ##
        ##                                 ##
        #####################################

        Data parameters:
          local_fk_label: data/BRP_analysis
          local_detect_label: data/BRP_analysis

        Algorithm parameters:
          window_len: 3600.0
          p_value: 0.99
          min_duration: 10.0
          back_az_width: 15.0
          fixed_thresh: None
          return_thresh: False

        Running fd...
        Writing detections to data/BRP_analysis.dets.json

    As noted in the output, a new file named :code:`BRP_analysis.dets.json` is created containing all of the detections identified in the fk results.  This file contains the information summarizing each detection in a format that can be ingested for further CLI analysis and can also be loaded into the :ref:`infraview` GUI.  The first detection from this analysis of the included BRP data is shown below:

    .. code-block:: none

        [
            {
            "Name": "",
                "Time (UTC)": "2012-04-09T18:10:17.008300",
                "F Stat.": 13.4034,
                "Trace Vel. (m/s)": 335.08,
                "Back Azimuth": -111.3,
                "Latitude": 39.4727,
                "Longitude": -110.741,
                "Elevation (m)": null,
                "Start": 0.0,
                "End": 15.0,
                "Freq Range": [
                    1.0,
                    10.0
                ],
                "Array Dim.": 4,
                "Method": "",
                "Event": "",
                "Note": "InfraPy CLI detection"
            }, ...


- In some cases, the parameters in the detection analysis are modified without changing the beamforming configuration and the :code:`run_fd` is useful in such scenarios.  However, most of the time, the beamforming and detection analysis are run together.  This can be accomplished in the InfraPy CLI via the :code:`run_fkd` option.  

    .. code-block:: bash
    
        infrapy run_fkd --config-file BRP_analysis.config

    This option essentially combines the :code:`run_fk` and :code:`run_fd` options into a single analysis run.

- In addition to analysis of local data, InfraPy's use of :code:`obspy.clients.fdsn` methods enables analysis of data available on IRIS and similar FDSNs.  Instead of specifying local waveform files, this requires defining the FDSN (e.g., IRIS, USGS) as well as the network, station, channel, and location information of the array.  Lastly, the start and end time are also needed to identify the segment of data to download for analysis.  This information can be entered on the command line, but it's easier to simply write up a config file in most cases (recall that individual parameters can be overwritten on the command line, so the station or start/end times can be modified as needed).  An example analysis from the IMS I53US array can be specified as:

    .. code-block:: none

        [WAVEFORM IO]
        fdsn = IRIS
        network = IM
        station = I53*
        location = *
        channel = *DF
        starttime = 2018-12-19T01:00:00
        endtime = 2018-12-19T03:00:00

        [DETECTION IO]
        local_fk_label = I53US_analysis
        local_detect_label = I53US_analysis

    Although not yet included in the CLI methods, an FDSN station browser is available in the :ref:`infraview` GUI to search for available data given a reference location, radius, and time bounds.

- Analysis of data from a local database is also available through the InfraPy CLI, and is covered in a separate tutorial on :ref:`pisces`.

--------------------
Network-Level Analyses
--------------------

- Once fk and fd analysis are run and detections are identified across a network of infrasound arrays, event identification and localization can be completed.  The detection set used in the Blom et al. (2020) evalution of a pair-based, joint-likelihood association algorithm are included as an example to demonstrate these analysis steps.  Detection files are in the examples/data/Blom_etal_2020/ directory and contain detections on each of 4 regional array in the western US (see the manuscript for a full discussion of the generation of this synthetic data set).  Analysis of these detections and identification of events can be completed by running:

    .. code-block:: bash
    
        infrapy run_assoc --local-detect-label 'data/Blom_etal2020_GJI/*' --local-event-label example

    Note that once again quotes are needed to define multiple files for ingestion.  This analysis can be on the slow side, so it's recommended to add on a :code:`--cpu-cnt` option and multithread the computation of the joint-likelihood values.  The analysis results will be summarized to the screen,

    .. code-block:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##       Association Analysis      ##
        ##                                 ##
        #####################################

        Data summary:
          local_detect_label: data/Blom_etal2020_GJI/*
          local_event_label: example
          starttime: None
          endtime: None

        Parameter summary:
          back_az_width: 10.0
          range_max: 2000.0
          resolution: 180
          distance_matrix_max: 8.0
          cluster_linkage: weighted
          cluster_threshold: 5.0
          trimming_threshold: 3.8

        Loading detections from files:
        	data/Blom_etal2020_GJI/NVIAR.dets.json
        	data/Blom_etal2020_GJI/I57US.dets.json
        	data/Blom_etal2020_GJI/DLIAR.dets.json
        	data/Blom_etal2020_GJI/PDIAR.dets.json

        Running event identification for: 2010-01-01T09:35:59.773Z - 2010-01-01T13:23:14.773Z
        	Computing joint-likelihoods...
		        Progress: 	[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
        	Clustering detections into events...
        	Trimming poor linkages and repeating clustering analysis...

        Running event identification for: 2010-01-01T10:51:44.773Z - 2010-01-01T14:38:59.773Z
        	Computing joint-likelihoods...
        		Progress: 	[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
	        Clustering detections into events...
        	Trimming poor linkages and repeating clustering analysis...

        Running event identification for: 2010-01-01T12:07:29.773Z - 2010-01-01T15:54:44.773Z
        	Computing joint-likelihoods...
        		Progress: 	[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
        	Clustering detections into events...
        	Trimming poor linkages and repeating clustering analysis...

        Running event identification for: 2010-01-01T13:23:14.773Z - 2010-01-01T17:10:29.773Z
        	Computing joint-likelihoods...
		        Progress: 	[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
        	Clustering detections into events...
	        Trimming poor linkages and repeating clustering analysis...

        Cleaning up and merging clusters...

    The analysis breaks the detection list into segments defined by the maximum propagation distance allows in order to avoid including detections in one analysis that will not be associated with others due to differences in detection times and typical infrasonic propagation velocities.  For each event identified in the analysis, a new .dets.json file is written that includes the subset of the original detections that have been identified as originating from a common event.  The naming convention of these files is :code:`local_event_label_ev-#.dets.json` and the example analysis here should have identified 3 events.

- Detection sets can be visualized on a map using the :code:`plot_dets` option.  This is useful in determining a useful maximum range for event identification and localization analysis.  For the above analysis of the Blom et al. (2020) synthetic data set, the full data set can be visualized with,

    .. code-block:: bash
    
        infrapy plot_dets --local-detect-label 'data/Blom_etal2020_GJI/*'

    .. image:: _static/_images/plot_dets1.png
        :width: 1200px
        :align: center

    This result is rather busy, but plotting each individual event's detections shows that the association algorithm correctly identified the events,

    .. code-block:: bash

        infrapy plot_dets --local-detect-label 'example1-ev0.dets.json'  --range-max 1000

    .. image:: _static/_images/plot_dets2.png
        :width: 1200px
        :align: center


- Once an event has been identified, the detections can be analyzed using the Bayesian Infrasonic Source Localization (BISL) methods as discussed in Blom et al (2015).  This requires specifying the detection list file as well as an output location file label,

    .. code-block:: bash

        infrapy plot_dets --local-detect-label example1-ev0  --local-loc-label example1-ev0

    The analysis steps are udpated as localization is performed and the resulting location and origin time information is printed to screen as well as written into an output file (the output file for InfraPy's localization is also a .json format file, but it's naming convention uses ".loc.json" to distinguish it from a detection result file)

    .. code-block:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##      Localization Analysis      ##
        ##                                 ##
        #####################################

        Data summary:
          local_event_label: example1-ev0
          local_loc_label: example1-ev0

        Parameter summary:
          back_az_width: 10.0
          range_max: 2000.0
          resolution: 180
          src_est: None
          pgm_file: None

        Loading detections from file: example1-ev0.dets.json

        Running Bayesian Infrasonic Source Localization (BISL) Analysis...
        	Identifying integration region...
        	Computing marginalized spatial PDF...
        	Computing confidence ellipse parameters...
        	Computing marginalized origin time PDF...

        BISL Summary:
        Maximum a posteriori analysis: 
        	Source location: 41.092, -113.144 
        	Source time: 2004-06-02T17:20:14.150 
        Source location analysis:
	        Latitude (mean and standard deviation): 40.977 +/- 29.236 km. 
	        Longitude (mean and standard deviation): -113.256 +/- 30.177 km.
	        Covariance: -0.139.
	        Area of 95 confidence ellipse: 16606.375 square kilometers
        Source time analysis:
        	Mean and standard deviation: 2004-06-02T17:19:14.339 +/- 99.051 second
        	Exact 90% confidence bounds: [2004-06-02T17:16:26.109, 2004-06-02T17:21:50.533]

        Writing localization result into example1-ev0.loc.json

- The localization result can be visualized in a number of ways.  Firstly, the detecting arrays and location estimate can be plotted on map using,

    .. code-block:: bash

        infrapy plot_loc --local-detect-label example1-ev0 --local-loc-label example1-ev0 --range-max 1000.0

    .. image:: _static/_images/plot_loc1.png
        :width: 1200px
        :align: center

    For Visualization of the source region in more detail, the :code:`--zoom` option can be set to true and the map zooms in to show only the estimated source region.

    .. code-block:: bash

        infrapy plot_loc --local-detect-label example1-ev0 --local-loc-label example1-ev0 --range-max 1000.0 --zoom true

    .. image:: _static/_images/plot_loc2.png
        :width: 900px
        :align: center

    Lastly, the origin time is estimated as part of the BISL analysis and can be visualized as,

    .. code-block:: bash

        infrapy plot_origin_time --local-loc-label example1-ev0 


    .. image:: _static/_images/plot_origin_time.png
        :width: 1200px
        :align: center

*************************************
Scripting and Notebook-Based Analysis 
*************************************

- A series of scripts illustrating how to use infrapy subroutines as stand-alone modules are found in the /examples folder.
The jupyter notebook documenting these steps is found in /tutorials/Quick Start.ipynb.  The notebook can be run by installing jupyter notebook via conda.

.. code-block:: python

    >> conda install jupyter notebook

*Example Analysis Scripts:*

+--------------------------------+------------------------------------------------------+
| examples/test_beamforming.py   | Run beamforming (fk) analysis on an Obspy stream     |
+--------------------------------+------------------------------------------------------+
| examples/test_slowness-grid.py | Run beamforming (fk) analysis on an Obspy stream     |
+--------------------------------+------------------------------------------------------+
| examples/test_detection.py     | Run beamforming (fk) analysis on an Obspy stream     |
+--------------------------------+------------------------------------------------------+
| examples/test_assoc.py         | Run beamforming (fk) analysis on an Obspy stream     |
+--------------------------------+------------------------------------------------------+
| examples/test_bisl.py          | Run beamforming (fk) analysis on an Obspy stream     |
+--------------------------------+------------------------------------------------------+
| examples/test_yield.py         | Run beamforming (fk) analysis on an Obspy stream     |
+--------------------------------+------------------------------------------------------+



1. Run Bartlett, Capon or Generalized Least Squares beamforming processes on an hour-long dataset from the BRP array in Utah

.. code-block:: python

    >> python test_beamforming.py

2. Visualize beamforming results in the sx/sy space

.. code-block:: python

    >> python test_slowness-grid.py

Detection:

1. Run detection on the series of beamforming results produced in the above step

.. code-block:: python

    >> python test_detection.py

Association

1. Associate a number of detections contained in a .dat file (/data/detection_set1.dat or /data/detection_set2.dat)

.. code-block:: python

    >> python test_assoc.py

Location

1. Test the Bayesian Infrasonic Source Localization (BISL) methodology using a set of provided detections (/data/detection_set1.dat or /data/detection_set2.dat).  Location will be run twice, once assuming uniform atmospheric propagation and a second time applying provided atmospheric propagation priors for the Western US (see Blom et al., 2015 for further explanation)

.. code-block:: python

    >> python test_bisl.py
