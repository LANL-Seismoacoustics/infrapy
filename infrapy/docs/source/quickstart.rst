.. _quickstart:

==========
Quickstart
==========

Most of InfraPy's analysis methods are accessible through a command line interface (CLI) with parameters specified either via command line flags or a configuration file.  Waveform data can be ingested from local files (eg., SAC or similar format that can be ingested via :code:`obspy.core.read`) or downloaded from FDSN clients via :code:`obspy.clients.fdsn`.  Detection and event analyses can be performed from the command line enabling a full pipeline of analysis from beamforming/detection to event identification and localization.  Visualization methods are also included to quickly interrogate analysis results.  The Quickstart summarized here steps through these various CLI methods and demonstrates the usage of InfraPy from the command line.

------------------
Detection Analyses
------------------
- Infrasonic signatures can be detected using the InfraPy algorithm on single channels or using coherence analysis across a spatially distributed set of sensors.  Detection methods are access usedin :code:`infrapy detect` and output into detection JSON files (e.g., '[ANALYSIS_ID].det.json.gz').

- Beamforming methods in InfraPy can be run via the :code:`infrapy detect beam` CLI option.  For a local data source such as the included SAC files in the data directory, this is simply,

    .. code-block:: bash

        infrapy detect beam --local-wvfrms 'data/YJ.BRP*.SAC' --cpu-cnt 4

    Note that the data path must be in quotes in order to be properly parsed and that this Quickstart assumes you are in the infrapy/examples directory (if you are getting an error that the waveform data isn't found, make sure you're in the correct directory).  As the methods are run, data and algorithm parameters are summarized and a progress bar shows how much of the data has been analyzed:

    .. code-block:: none

        ######################################
        ##                                  ##
        ##              InfraPy             ##
        ##  Beamforming Detection Analyses  ##
        ##                                  ##
        ######################################


        Data parameters:
          local_wvfrms: data/YJ.BRP*.SAC
          local_latlon: None
          det_label: None

        fk (beam) parameters:
          freq_min: 0.5
          freq_max: 5.0
          back_az_min: -180.0
          back_az_max: 180.0
          back_az_step: 2.0
          trace_vel_min: 300.0
          trace_vel_max: 600.0
          trace_vel_step: 2.5
          method: bartlett
          window_len: 10.0
          window_step: 5.0
          cpu_cnt: 4

        detection parameters:
          window_len: 3600.0
          p_value: 0.01
          min_duration: 10.0
          back_az_width: 15.0
          merge_dets: False

        Loading local data from data/YJ.BRP*.SAC

        Data summary:
        YJ.BRP1..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP2..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP3..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z
        YJ.BRP4..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z

        Running fk analysis...
            Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
        Running adaptive f-detector...

        Writing beamforming and detection results into data/YJ.BRP_2012.04.09T18.00.00.dets.json.gz


- Once completed, this analysis produces an output file containing a summary of the analysis and any identified detection, :code:`data/YJ.BRP_2012.04.09T18.00.00.dets.json.gz`. The naming convention of the output file uses the network, station, and time associated with the waveform data, but can be overwritten via the :code:`--det-label` parameter.

- Uncompressing and interrogating this file, one can find the waveform info and parameter information, computed fk (i.e., direction-of-arrival and Fisher statistic values) for the entire analysis duration, as well as information for individual detections.

    .. code-block:: none

       {
       "wvfrm_info": [
            [{
                "trace id": "YJ.BRP1..EDF",
                "starttime": "2012-04-09T18:00:00.008300Z",
                "endtime": "2012-04-09T18:19:59.998300Z",
                "latitude": 39.47269821166992,
                "longitude": -110.74089813232422
            }, ...]
        ],
        "fk_params": [
            {"freq_min": 0.5, "freq_max": 5.0, ...}
        ],
        "det_params": [
            {"window_len": 3600.0, "p_value": 0.01, ...}
        ],
        "fk": {
            "time": [5.0, 10.0, ..., 1185.0, 1190.0],
            "back az": [-138.72961493387828, -102.31457933170114, ...,-66.15065005843925, 1.032026780110567],
            "tr vel": [299.31252980459686, 500.4425296042015, ..., 356.0196938960705, 298.776690619962],
            "f-stat": [1.7871106151747223, 1.432899229285171, ..., 1.8014862125127102, 2.2871847705410375],
            "thresh": [3.877682287563356, 3.877682287563356, ..., 3.877682287563356, 3.877682287563356]
        },
        "det_info": [
            {
                "peak f-stat time": "2012-04-09T18:10:20.008300",
                "start/end": [
                    [-5.0, 15.0]
                ],
                "f-stat": 7.598484850361515,
                "back az": -108.51593202520215,
                "tr vel": 334.9964881969233,
                "fk": [
                    {
                        "time": [-20.0, -15.0, ..., 25.0, 30.0],
                        "back az": [-111.43647403399356, -111.71492846584691, ..., -110.81575881997433, -102.20248576801993],
                        "tr vel": [336.17337928063006, 347.1485162377091, ..., 376.06290789009705, 341.4288289753299],
                        "f-stat": [5.053869620047696, 2.8133754341174946, ..., 3.0936454658111896, 2.6011926111268546]
                    }
                ],
                "beam": [
                    {
                        "time": [-20.0, -19.99, ..., 29.98, 29.99],
                        "signal": [-118.61456758973551, -154.31485492857715, ..., -77.2684996439051, -67.25749776831691],
                        "resid": [381.6423602015375, 235.87965220970239, ..., 446.07608621837704, 506.4357230653115]
                    }
                ],
                "spec": [
                    {
                        "freq": [0.0, 0.048828125, ..., 49.951171875, 50.0],
                        "signal": [5.474957275388414, 436.0711990694136, ..., 2.216215682088617, 2.213116253142397],
                        "resid": [10.732711181641648, 267.84490402778556, ..., 1.532004515697567, 1.4349766456268886]
                    }
                ]
            },...

        ]


- The beamforming results from the :code:`infrapy detect beam` analysis can be visualized using :code:`infrapy plot detect` to show the full beam results or summaries for individual detections.  Providing just the detection visualizes the beam summary:

    .. code-block:: bash

        infrapy plot detect --dets-file data/YJ.BRP_2012.04.09T18.00.00.dets.json.gz


    .. image:: _static/_images/plot_beam1.png
        :width: 1500px
        :align: center

    The default behavior of the plotting methods in InfraPy are to generate a :code:`matplotlib` window and print the image to screen.  This can be overwritten by specifying an output file and turning the print to screen off:

    .. code-block:: bash

        infrapy plot detect --dets-file data/YJ.BRP_2012.04.09T18.00.00.dets.json.gz --figure-out "BRP_beam.png" --show-figure false

    Summaries of specific detections identified in the analysis can be visualized by specifying a detection index number (the colored segments in the above result),

    .. code-block:: bash

        infrapy plot detect --dets-file data/YJ.BRP_2012.04.09T18.00.00.dets.json.gz --det-index 1

    .. image:: _static/_images/plot_beam2.png
        :width: 1200px
        :align: center

    This plot summarizes the waveform data used as well as the detection time, f-stat, direction of arrival, duration, and frequency band used in analysis.

- In many cases the frequency band for a signal of interest or the window length appropriate for a given frequency band needs to be modified.  From the command line, this can be done by specifying a number of options in the algorithm as summarized in the :code:`--help` (:code:`-h`) information.  For example, the analysis of data from BRP can be completed using a higher frequency band and shortened analysis windows via:

    .. code-block:: bash

        infrapy detect beam --local-wvfrms 'data/YJ.BRP*.SAC' --freq-min 1.0 --freq-max 8.0 --fk-window-len 4.0 --fk-window-step 2.0 --cpu-cnt 4

    If you've run the above after already running the previous example, you'll get a response that looks like this:

    .. code-block:: bash

        ######################################
        ##                                  ##
        ##              InfraPy             ##
        ##  Beamforming Detection Analyses  ##
        ##                                  ##
        ######################################


        Data parameters:
          local_wvfrms: data/YJ.BRP*.SAC
          local_latlon: None
          det_label: None

        ...

        WARNING!!! Detection results file (data/YJ.BRP_2012.04.09T18.00.00.dets.json.gz) already exists and will be overwritten.
        Do you want to proceed? (y/n):

    The intended output :code:`[...].dets.json.gz`` file for this analysis already exists from a previous run, and so a catch is included to warn and possibly prevent this overwriting.  Though not generally recommended, a CLI flag is available to automatically overwrite any existing results file (:code:`--auto-overwrite True`).

- In the case that multiple analysis parameters are changed from their default values, a configuration file is useful to simplify running analysis and keep a record of what was used for future review of analysis.  Within the :code:`examples/config` directory are several example configuration files.  The :code:`detection_local.config` file has a configuration to run detection (fk and fd) analysis on local waveform data:

    .. code-block:: none

        [DATA IO]
        local_wvfrms = data/YJ.BRP*.SAC
        det_label = BRP_hf

        [FK]
        freq_min = 2.0
        freq_max = 8.0
        window_len = 4.0
        window_step = 2.0

        [FD]
        p_value = 0.02
        min_duration = 20.0
        merge_dets = True

- Note that the parameter specifications use underscores in the config file and hyphens in the command line flags and also that headers in the config file are used to distinguish beamforming (FK) and f-stat detection (FD) parameters.  For example, the frequency minimum is specified using :code:`--freq-min` on the command line and :code:`freq_min` in a configuration file.  The beamforming window length is set by :code:`--fk-window-len` on the command line and as :code:`window_len` under :code:`[FK]` in the configuration file while the detection window length is set by :code:`--fd-window-len` on the command line and as :code:`window_len` under :code:`[FD]` in the configuration file.  There is a :code:`default.config` file in the InfraPy resources directory with a full set of all tunable parameters and which header they fall under.  Data input or output information is expected to be listed under :code:`[DATA IO]`.

- Using the configuration file, the analysis can be completed by simply running:

    .. code-block:: bash

        infrapy detect beam --cnfg-file config/detection_local.config

    .. code-block:: bash

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##    Beamforming (fk) Analysis    ##
        ##                                 ##
        #####################################

        Loading configuration info from: config/detection_local.config

        Data parameters:
          local_wvfrms: data/YJ.BRP*.SAC
          local_latlon: None
          det_label: BRP_hf

        ...

        Running beamforming analysis...
	        Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]

        Writing beamforming and detection results into BRP_hf.dets.json.gz

- In addition to ingesting local :code:`.SAC` or similar files (anything that can read by ObsPy), the FDSN Client in ObsPy can be called to pull data.  This is best done via a configuration file given the FDSN, network, station, location, channel, and start/end times of the data of interest must all be specified for the data pull.  The :code:`[DATA IO]` block to run a beam using data from the I53 infrasound station is,


    .. code-block:: bash

        [DATA IO]
        fdsn = IRIS
        network = IM
        station = I53*
        location = *
        channel = *DF
        starttime = 2018-12-19T01:00:00
        endtime = 2018-12-19T03:00:00

- An example analysis using this time period (during which a bolide signal was recorded on the station) can be completed using the included FDSN configuration file,

    .. code-block:: bash

        infrapy detect beam --cnfg-file config/detection_fdsn.config --cpu-cnt 4


    .. code-block:: bash

        ######################################
        ##                                  ##
        ##              InfraPy             ##
        ##  Beamforming Detection Analyses  ##
        ##                                  ##
        ######################################

        Loading configuration info from: config/detection_fdsn.config

        Data parameters:
          fdsn: IRIS
          network: IM
          station: I53*
          location: *
          channel: *DF
          starttime: 2018-12-19T01:00:00
          endtime: 2018-12-19T03:00:00
          det_label: None

        fk (beam) parameters:
          freq_min: 0.08
          freq_max: 2.0
          back_az_min: -180.0
          back_az_max: 180.0
          back_az_step: 2.0
          trace_vel_min: 300.0
          trace_vel_max: 600.0
          trace_vel_step: 2.5
          method: bartlett
          window_len: 20.0
          window_step: 10.0
          cpu_cnt: 4

        detection parameters:
          window_len: 3600.0
          p_value: 0.05
          min_duration: 40.0
          back_az_width: 15.0
          merge_dets: True

        Loading data from FDSN (IRIS)...

        Data summary:
        IM.I53H1..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        IM.I53H2..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        IM.I53H3..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        IM.I53H4..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        IM.I53H5..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        IM.I53H6..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        IM.I53H7..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z
        IM.I53H8..BDF	2018-12-19T01:00:00.000000Z - 2018-12-19T03:00:00.000000Z

        Running beamforming analysis...
            Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
        Running adaptive f-detector...

        Writing beamforming and detection results into IM.I53H_2018.12.19T01.00.00.dets.json.gz

- The above beamforming detection analysis requires a spatially distributed set of microbarometers and uses the coherence to identify signals of interest.  When only a single data stream is available, a spectrogram-based detection method can be applied to identify possible signals of interest. The spectrogram-based detection methods can be accessed using the :code:`infrapy detect spectral` CLI option.  For a local data source such as the included SAC files in the data directory, this is simply (note: if multiple files are selected, only the first file will be used given the single-channel nature of this analysis),

    .. code-block:: bash

        infrapy detect spectral --local-wvfrms 'data/YJ.BRP1..EDF.SAC' --cpu-cnt 4

    .. code-block:: bash

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##   Spectral Detection Analyses   ##
        ##                                 ##
        #####################################


        Data parameters:
          local_wvfrms: data/YJ.BRP1..EDF.SAC
          local_latlon: None
          det_label: None
          cpu_cnt: 4

        sd (spectral detector) parameters:
          spectral_option: spectrogram
          freq_min: 1.0
          freq_max: 20.0
          window_len: 900.0
          window_step: 450.0
          p_value: 0.01
          freq_tm_factor: 35.0
          cluster_eps: 10.0
          cluster_min_samples: 40
          cluster_window_len: 600.0
          cpu_cnt: 4

        Loading local data from data/YJ.BRP1..EDF.SAC

        Data summary:
        YJ.BRP1..EDF	2012-04-09T18:00:00.008300Z - 2012-04-09T18:19:59.998300Z

        Running spectral detection (sd) analysis...
            Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
        Clustering into detections...
        Identified 5 detections.

- The resulting spectral analysis can be visualized similarly to a beam-based detection result and shows the original spectrogram, the normalized spectrogram with the background removed, and the clustered time-frequency points that were grouped into detections.  Logic is written into the detection plotting method to identify whether beam- or spectral-based parameters are saved in the file, so it will automatically generate the appropriate visualization.

    .. code-block:: bash

            infrapy plot detect --dets-file data/YJ.BRP1_2012.04.09T18.00.00.dets.json.gz

    .. image:: _static/_images/plot_spec1.png
        :width: 1200px
        :align: center

- Once again, the information related to an individual detection can be visualized by specifying a detection index:

    .. code-block:: bash

            infrapy plot detect --dets-file data/YJ.BRP1_2012.04.09T18.00.00.dets.json.gz --det-index 2

    .. image:: _static/_images/plot_spec2.png
        :width: 1200px
        :align: center


- Detections obtained from this spectral analysis do not include direction-of-arrival information (e.g., back azimuth and trace velocity), but the detection time, waveform, spectrogram, averaged signal spectra curve, and background spectra curve are all written into the resulting JSON detection file.

    .. code-block:: none

       {
       "wvfrm_info": [
            {
                "trace id": "YJ.BRP1..EDF",
                "starttime": "2012-04-09T18:00:00.008300Z",
                "endtime": "2012-04-09T18:19:59.998300Z",
                "latitude": 39.47269821166992,
                "longitude": -110.74089813232422
            }
        ],
        "sd_params": [
            {"spectral_option": "spectrogram", ...}
        ],
        "spectrogram": [...],
        "history:", [...],
        "det_info": [
            {
                "peak f-stat time": "2012-04-09T18:10:20.008300",
                "waveform": [...],
                "spec pnts": [...],
                "spectrogram": [...],
                "spec": [...],
                "bg spec": [...]
            }, ...
        ]

- Depending on the width of the frequency band of the signal, it might be more insightful to visualize the detection results using line or logarithmic scaled frequency.  The default behavior of the visualization uses linear scaling, but logarithic scaling can be applied by including the option :code:`--log-scale-freq True``.

    .. image:: _static/_images/plot_spec3.png
        :width: 1200px
        :align: center

    .. image:: _static/_images/plot_spec4.png
        :width: 1200px
        :align: center

--------------
Event Analyses
--------------

- Once detections are identified across a network of infrasound arrays, event identification, localization, and characterization can be completed.  The various event-level analysis methods are accessible through :code:`infrapy event`.  The detection set used in the Blom et al. (2020) evaluation of a pair-based, joint-likelihood association algorithm are included as an example to demonstrate these analysis steps.  Detection files are in the examples/data/Blom_etal_2020/ directory and contain detections on each of 4 regional array in the western US (see the manuscript for a full discussion of the generation of this synthetic data set).  Analysis of these detections and identification of events can be completed by running:

    .. code-block:: bash

        infrapy event build --dets-files 'data/Blom_etal2020_GJI/SY*dets.json.gz' --ev-label Blom_etal2020_GJI --range-max 1500.0 --cpu-cnt 4

    Note that once again quotes are needed to define multiple files for ingestion.  This analysis can be on the slow side, so it's recommended to use the :code:`--cpu-cnt` option and multithread the computation of the joint-likelihood values.  For this analysis, multi-threading distributes the individual joint-likelihood calculations between pairs of detections to available threads.  The analysis results will be summarized to the screen,

    .. code-block:: none

        ####################################
        ##                                ##
        ##             InfraPy            ##
        ##         Event Building         ##
        ##                                ##
        ####################################


        Data summary:
          dets_files: data/Blom_etal2020_GJI/SY*dets.json.gz
          ev_file: Blom_etal2020_GJI

        association parameters:
          celerity_model: regional_lf
          back_az_width: 10.0
          range_max: 1500.0
          resolution: 180
          distance_matrix_max: 8.0
          cluster_linkage: weighted
          cluster_threshold: 5.0
          trimming_threshold: 3.6
          ev_population_min: 3
          ev_station_min: 2
          cpu_cnt: 4


        Running event identification for: 2010-01-01T10:13:51.773000Z - 2010-01-01T13:04:18.773000Z
            Computing joint-likelihoods...
                Progress: 	[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
            Clustering detections into events...
            Trimming poor linkages and repeating clustering analysis...

        Running event identification for: 2010-01-01T11:10:40.773000Z - 2010-01-01T14:01:07.773000Z
            Computing joint-likelihoods...
                Progress: 	[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
            Clustering detections into events...
            Trimming poor linkages and repeating clustering analysis...

        Running event identification for: 2010-01-01T12:07:29.773000Z - 2010-01-01T14:57:56.773000Z
            Computing joint-likelihoods...
                Progress: 	[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
            Clustering detections into events...
            Trimming poor linkages and repeating clustering analysis...

        Running event identification for: 2010-01-01T13:04:18.773000Z - 2010-01-01T15:54:45.773000Z
            Computing joint-likelihoods...
                Progress: 	[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
            Clustering detections into events...
            Trimming poor linkages and repeating clustering analysis...

        Cleaning up and merging clusters...
        Identified 3 event(s).

    The analysis breaks the detection list into segments defined by the maximum propagation distance allows in order to avoid including detections in one analysis that will not be associated with others due to differences in detection times and typical infrasonic propagation velocities.

- For each cluster of detections identified in the analysis, an event JSON file output is written that includes the subset of the original detections which have been identified as originating from a common event.  The naming convention of these files is :code:`ev_label_ev-#.ev.json.gz` and the example analysis here should have identified 3 events.  In addition to the detection list for the event, the association analysis parameters and distance matrix of the event (see Blom et al., 2020 for information about the distance matrix) are saved for re-producibility and to document the event cluster quality.  Empty entries for ground truth information, location results, and characterization (yield estimation) results are also defined when the file is created.  Also note that since the source waveform data and detection parameters can be distinct for different stations and signals, those entries are copied into the individual detection entries when writing information into the event file.

    .. code-block:: none

        {
            "ground truth": {},
            "det_info": [
                {
                    "wvfrm_info": [...],
                    "fk_params": [...],
                    "fk_params": [...],
                    "det_params": [...],
                    "peak f-stat time": "2010-01-01T12:58:17",
                    "start/end":[...],
                    "f-stat": 20.0,
                    "back az": -77.1,
                    "tr vel": 330.0
                }, ...
            ],
            "assoc_params": {
                "celerity_model": 'regional_lf',
                "back_az_width": 10.0,
                ...
            },
            "dist_matrix": [
                [...],
                [...],
                ...
            ],
            "location": []
            "characterization": []
        }

- Detection sets can be visualized on a map using the :code:`plot map_dets` option.  This is insightful in determining a maximum range for event identification and localization analysis.  The synthetic detection set used in the Blom et al. (2020) evaluation of the event building algorithm can be visualized as by defining the :code:`--dets-files`

    .. code-block:: bash

        infrapy plot map_dets --dets-files 'data/Blom_etal2020_GJI/SY.DLIAR_2010.01.01T12.00.00-14.00.00.dets.json.gz' --range-max 1500

    .. image:: _static/_images/map_dets-DLIAR.png
        :width: 1200px
        :align: center

    The full set of detections used in the analysis can be specified using wild cards to ingest multiple detection files,

    .. code-block:: bash

        infrapy plot map_dets --dets-files 'data/Blom_etal2020_GJI/SY.*.dets.json.gz' --range-max 1500

    .. image:: _static/_images/map_dets-all.png
        :width: 1200px
        :align: center

    Finally, the detections included in a given event files can be visualized by specifying :code:`--ev-file` instead of a detections file.  Note that this visualization can be made before any localization analysis has been applied to inform such analyses or as quality control before attempting localization (to ensure a spurious detection hasn't been included).

    .. code-block:: bash

        infrapy plot map_dets --ev-file data/Blom_etal2020_GJI/Blom_etal2020_GJI-0.ev.json.gz

    .. image:: _static/_images/map_dets-ev0.png
        :width: 1200px
        :align: center


- Once an event has been identified, the detections can be analyzed using the Bayesian Infrasonic Source Localization (BISL) methods as discussed in Blom et al. (2015) or using the Time-Reversed Infrasonic Bayesian Localization (TRIBL) algorithm more recently developed in Blom et al. (2025).  Localization analysis is run using, :code:`infrapy event locate`,

    .. code-block:: bash

        infrapy event locate --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz

    The default set of parameters uses the maximum range from the event building stage (1500 km for this example) and a 10.0 degree back azimuth width to identify the spatial region for analysis.  The default celerity model used in analysis is tuned for regional low-frequency signals (including contributions for tropospheric, stratospheric, and thermospheric waveguide travel times).

    .. code-block:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##      Localization Analysis      ##
        ##                                 ##
        #####################################


        Data summary:
          ev_file: data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz

        localization parameters:
          back_az_width: 10.0
          range_max: 1500.0
          grid_resol: 180
          celerity_model: regional_lf


        Running Bayesian Infrasonic Source Localization (BISL) Analysis...
            Identifying integration region...
            Evaluating localization probability on grid...
                Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
            Analyzing localization pdf...
                Normalizing and marginalizing...
                Analyzing spatial PDF...
                Analyzing temporal PDF...

        Localization Summary:
        Maximum a posteriori analysis:
            Source location: 41.489, -112.015
            Source time: 2010-01-01T12:10:11.515500
        Source location analysis:
            Latitude (mean and standard deviation): 41.341 +/- 25.451 km.
            Longitude (mean and standard deviation): -112.119 +/- 25.387 km.
            Covariance: 0.428.
            Area of 90% confidence ellipse: 9347.869 square kilometers
        Source time analysis:
            Mean and standard deviation: 2010-01-01T12:08:04.896500 +/- 125.167 second
            Exact 90% confidence bounds: [2010-01-01T12:04:58.593500, 2010-01-01T12:11:27.064500]


    The localization result is written into the existing event file and can be visualized by specifying the location solution index.  If no location index is provided, only the back azimuth projections are visualized reproducing the above visualization from :code:`maps_dets` and a summary of available event analysis results is printed to screen.


        .. code-block:: bash

            infrapy plot event --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz

        .. code-block:: none


            ###############################
            ##                           ##
            ##          InfraPy          ##
            ##    Event Visualization    ##
            ##                           ##
            ###############################


            Data summary:
              ev_file: data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz
              loc_index: None (1 location result(s) in file
              char_index: None (0 characterization result(s) in file


            Drawing map with detection back azimuth projections...


    Specifying the first localization result in the visualization command produces a figure summarizing the localization result and a summary printed to screen matching the above informaiton when the analysis was initially run.

    .. code-block:: bash

        infrapy plot event --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz  --loc-index 0

    .. image:: _static/_images/plot_ev1.png
        :width: 1200px
        :align: center


    Several different celerity (horizonal group velocity) models are built into InraPy and can be used by specifying them via the :code:`--celerity-model` option.  The Infrasound Global Empirical Model (infGEM) that is tuned for use in larger propagation distance scenarios (see Nippress et al, 2023) can be used and the result visualized via,

    .. code-block:: bash

        infrapy event locate --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz --celerity-model infgem


    .. code-block:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##      Localization Analysis      ##
        ##                                 ##
        #####################################


        Data summary:
        ev_file: data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz

        localization parameters:
          back_az_width: 10.0
          range_max: 1500.0
          grid_resol: 180
          celerity_model: infgem


        Running Bayesian Infrasonic Source Localization (BISL) Analysis...
            Identifying integration region...
            Evaluating localization probability on grid...
                Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
            Analyzing localization pdf...
                Normalizing and marginalizing...
                Analyzing spatial PDF...
                Analyzing temporal PDF...

        Localization Summary:
        Maximum a posteriori analysis:
            Source location: 41.357, -112.037
            Source time: 2010-01-01T12:10:11.515500
        Source location analysis:
            Latitude (mean and standard deviation): 41.349 +/- 14.105 km.
            Longitude (mean and standard deviation): -111.952 +/- 14.286 km.
            Covariance: -0.273.
            Area of 90% confidence ellipse: 2915.259 square kilometers
        Source time analysis:
            Mean and standard deviation: 2010-01-01T12:10:35.531500 +/- 43.012 second
            Exact 90% confidence bounds: [2010-01-01T12:09:30.144500, 2010-01-01T12:11:47.106500]


    Once again, visualization without a specified location index just plots the back projections and summarizes the available informaiton (2 location results now).  Specifying this new location index,

    .. code-block:: bash

        infrapy plot event --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz  --loc-index 1

    .. image:: _static/_images/plot_ev2.png
        :width: 1200px
        :align: center

- The Time-Reversed Infrasonic Bayesian Localization (TRIBL) methods require ray tracing back projection paths through an atmosphere specification and a more limited grid to evaluate results on.  Using the general region identified by BISL and accounting for some cross wind deviations, one can define the lower-left (ll) and upper-right (ur) corners as well as a range of origin times.  These various quantities could be defined item by item on the command line, but it's once again easier to save them into a configuration file and point the method at that.  An example configuration file is included here with the following informaiton:

    .. code-block:: bash

        [LOC]
        atmo_data = data/Blom_etal2024_GJI/g2stxt_2010010100_41.1310_-112.8960.dat

        range_max = 1250
        ll_corner = 40.5, -113.5
        ur_corner = 41.75, -112.0
        alt_bounds = 0, 0

        latlon_resol = 0.04
        alt_resol = 1.0

        tm_min = 2010-01-01T12:03:00
        tm_max = 2010-01-01T12:09:00
        tm_resol = 10.0

        grnd_snd_spd_stdev = 5.0
        det_tm_stdev = 5.0

        local_temp_dir = data/Blom_etal2024_GJI/temp


    TRIBL can be run using the :code:`infrapy event locate` method and pointing at the configuration file.  When the :code:`--atmo-data` parameter is specified, the implementation recognizes that atmospheric data is being used and applies the TRIBL algorithm instead of BISL.

    .. code-block:: bash

        infrapy event locate --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz --cnfg-file config/tribl_example.config

    .. code-block:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##      Localization Analysis      ##
        ##                                 ##
        #####################################


        Loading configuration info from: config/tribl_example.config

        Data summary:
        ev_file: data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz

        localization parameters:
          ll_corner: [  40.5 -113.5]
          ur_corner: [  41.75 -112.  ]
          latlon_resol: 0.04
          tm_min: 2010-01-01T12:03:00
          tm_max: 2010-01-01T12:09:00
          tm_resol: 10.0
          atmo_data: data/Blom_etal2024_GJI/g2stxt_2010010100_41.1310_-112.8960.dat
          alt_resol: 1.0
          c0_stdev: 5.0
          det_tm_stdev: 5.0
          az_limit: 2.0
          local_temp_dir: data/Blom_etal2024_GJI/temp


        Running Time-Reversed Infrasonic Bayesian Localization (TRIBL) Analysis...
            Identifying integration region and building grid...
            Computing back projections for detection list...
            Evaluating localization probability on grid...
                Progress: [>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>]
            Analyzing localization pdf...
                Normalizing and marginalizing...
                Analyzing spatial PDF...
                Analyzing temporal PDF...

        Localization Summary:
        Maximum a posteriori analysis:
            Source location: 41.18, -112.86
            Source time: 2010-01-01T12:06:20.000
        Source location analysis:
            Latitude (mean and standard deviation): 41.185 +/- 9.521 km.
            Longitude (mean and standard deviation): -112.854 +/- 8.363 km.
            Covariance: 0.094.
            Area of 90% confidence ellipse: 1151.876 square kilometers
        Source time analysis:
            Mean and standard deviation: 2010-01-01T12:06:21.542 +/- 28.952 second
            Exact 90% confidence bounds: [2010-01-01T12:05:33.634, 2010-01-01T12:07:08.872]


    Visualization is once again done through :code:`plot event` with the new location index specified.

    .. code-block:: bash

        infrapy plot event --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz  --loc-index 2

    .. image:: _static/_images/plot_ev3.png
        :width: 1200px
        :align: center

- Once these various location results completed, it's useful to interrogate the event file and identify what has been done.  A utility function is available to summarize the contents of an event file.  Running this on the event file we've generated localization results for,

    .. code-block:: bash

        infrapy utils ev_summary --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz

    .. code-block:: none

            #################################
            ##                             ##
            ##      InfraPy Utilities      ##
            ##     Summarize Event File    ##
            ##                             ##
            #################################

            Loading information from ev_file: data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz

            =================
            Detection Summary
            =================
            SY.PDIAR..BDF	  loc: 42.7670, -109.5940	  time: 2010-01-01T12:22:23	  back az [deg]:  -122.30	  tr vel [m/s]: 341.1	  f-stat:  25.0
            SY.PDIAR..BDF	  loc: 42.7670, -109.5940	  time: 2010-01-01T12:25:26	  back az [deg]:  -127.24	  tr vel [m/s]: 388.3	  f-stat:  25.0
            SY.PDIAR..BDF	  loc: 42.7670, -109.5940	  time: 2010-01-01T12:26:37	  back az [deg]:  -126.32	  tr vel [m/s]: 471.5	  f-stat:  25.0
            SY.DLIAR..BDF	  loc: 35.8570, -106.3150	  time: 2010-01-01T12:50:00	  back az [deg]:   -37.52	  tr vel [m/s]: 363.9	  f-stat:  25.0
            SY.DLIAR..BDF	  loc: 35.8570, -106.3150	  time: 2010-01-01T12:55:44	  back az [deg]:   -34.75	  tr vel [m/s]: 391.0	  f-stat:  25.0
            SY.NVIAR..BDF	  loc: 38.4300, -118.3040	  time: 2010-01-01T12:50:10	  back az [deg]:    58.03	  tr vel [m/s]: 365.5	  f-stat:  25.0


            ====================
            Localization Summary
            ====================

            ##############
            ## index: 0 ##
            ##############
            parameters
            ----------
                back_az_width: 10.0
                range_max: 1500.0
                grid_resol: 180
                celerity_model: regional_lf

            result
            ------
                latitude: 41.341 deg +/- 25.45 km.
                longitude: -112.119 deg +/- 25.39 km.
                origin time: 2010-01-01T12:08:04.896500 +/- 125.2 s.

            ##############
            ## index: 1 ##
            ##############
            parameters
            ----------
                back_az_width: 10.0
                range_max: 1500.0
                grid_resol: 180
                celerity_model: infgem

            result
            ------
                latitude: 41.349 deg +/- 14.1 km.
                longitude: -111.952 deg +/- 14.29 km.
                origin time: 2010-01-01T12:10:35.531500 +/- 43.0 s.

            ##############
            ## index: 2 ##
            ##############
            parameters
            ----------
                ll_corner: [40.5, -113.5]
                ur_corner: [41.75, -112.0]
                latlon_resol: 0.04
                tm_min: 2010-01-01T12:03:00
                tm_max: 2010-01-01T12:09:00
                tm_resol: 10.0
                atmo_data: data/Blom_etal2024_GJI/g2stxt_2010010100_41.1310_-112.8960.dat
                alt_resol: 1.0
                c0_stdev: 5.0
                det_tm_stdev: 5.0
                az_limit: 2.0
                local_temp_dir: data/Blom_etal2024_GJI/temp

            result
            ------
                latitude: 41.185 deg +/- 9.52 km.
                longitude: -112.854 deg +/- 8.36 km.
                origin time: 2010-01-01T12:06:21.542 +/- 29.0 s.


- When a location analysis is attempted, but results already exist in the file, the result is simply printed to screen.  For running examples, debugging, or related work, an event file reset utility is avilable.  This will reset the JSON file localization and/or characterization fields to empty lists, :code:`[]`.  A warning is given to confirm that this removal of existing analysis is desired,

    .. code-block:: bash

        infrapy utils ev_reset --ev-file data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz --reset-loc True


    .. code-block:: none

        #################################
        ##                             ##
        ##      InfraPy Utilities      ##
        ##       Reset Event File      ##
        ##                             ##
        #################################

        Loading information from ev_file: data/Blom_etal2024_GJI/SY.UTTR_2010.01.01T12.00.00.ev.json.gz
          3 location result(s) in file
          0 characterization result(s) in file

        ====================
        Localization Summary
        ====================

        ##############
        ## index: 0 ##
        ##############
        parameters
        ----------
            back_az_width: 10.0
            range_max: 1500.0
            grid_resol: 180
            celerity_model: regional_lf

        result
        ------
            latitude: 41.341 deg +/- 25.45 km.
            longitude: -112.119 deg +/- 25.39 km.
            origin time: 2010-01-01T12:08:04.896500 +/- 125.2 s.

        ##############
        ## index: 1 ##
        ##############
        parameters
        ----------
            back_az_width: 10.0
            range_max: 1500.0
            grid_resol: 180
            celerity_model: infgem

        result
        ------
            latitude: 41.349 deg +/- 14.1 km.
            longitude: -111.952 deg +/- 14.29 km.
            origin time: 2010-01-01T12:10:35.531500 +/- 43.0 s.

        ##############
        ## index: 2 ##
        ##############
        parameters
        ----------
            ll_corner: [40.5, -113.5]
            ur_corner: [41.75, -112.0]
            latlon_resol: 0.04
            tm_min: 2010-01-01T12:03:00
            tm_max: 2010-01-01T12:09:00
            tm_resol: 10.0
            atmo_data: data/Blom_etal2024_GJI/g2stxt_2010010100_41.1310_-112.8960.dat
            alt_resol: 1.0
            c0_stdev: 5.0
            det_tm_stdev: 5.0
            az_limit: 2.0
            local_temp_dir: data/Blom_etal2024_GJI/temp

        result
        ------
            latitude: 41.185 deg +/- 9.52 km.
            longitude: -112.854 deg +/- 8.36 km.
            origin time: 2010-01-01T12:06:21.542 +/- 29.0 s.

        ########################################
        ########################################

        WARNING!!! This action will remove existing result(s) in this event file.
        Do you want to proceed? (y/n): y

        Removing localization analysis result(s) from event file...

- Infrasonic signals produced by above-ground explosive sources can be used to estimate the explosive yield via source models such as the Kinney & Graham blastwave scaling laws. InfraPy's Spectral Yield Estimate (SpYE) methods can be applied to relate regional infrasonic signal spectral amplitude to a near-source estimate and then the blastwave model to estimate yield.  Usage of these methods requires an event file for which the detections have spectral amplitude data as well as transmission loss models (TLMs) generated using the *stochprop* software.  An example set of TLMs are included in the data/TLMs/ directory for the western US during the month of August.  These models were developed during testing of the SpYE methods and applied to estimate the yield of the Humming Roadrunner experiments (see Blom et al., 2018).  An event file for HRR-5 is included in the example data directory to demonstrate yield estimation.

    .. code:: bash

        infrapy event locate --ev-file data/HRR-5.ev.json.gz

    Visualizating the localization result shows that this event was detected by 6 stations off to the west and three to the north.

    .. image:: _static/_images/HRR5_loc.png
        :width: 1200px
        :align: center


    Printing the summary of this event to screen, we can identify which of the stations are to the north (not in the stratospheric waveguide) and which don't have good SNR,


    .. code:: none

        #################################
        ##                             ##
        ##      InfraPy Utilities      ##
        ##     Summarize Event File    ##
        ##                             ##
        #################################

        Loading information from ev_file: data/HRR-5.ev.json.gz

        =================
        Detection Summary
        =================
        2D.N30CN..CDF	  loc: 36.2494, -105.8898	  time: 2012-08-27T23:19:28	  back az [deg]:  -169.43	  tr vel [m/s]: 359.9	  f-stat:   7.4	  high snr band: 0.40 - 181.49 Hz
        2D.N34C..CDF	  loc: 36.5727, -105.7564	  time: 2012-08-27T23:23:51	  back az [deg]:  -168.06	  tr vel [m/s]: 357.2	  f-stat:  41.5	  high snr band: 0.10 - 1.05 Hz
        2D.W22CW..CDF	  loc: 34.3109, -108.6895	  time: 2012-08-27T23:14:56	  back az [deg]:   109.19	  tr vel [m/s]: 358.2	  f-stat: 380.4	  high snr band: 0.17 - 7.00 Hz
        2D.W24CE..CDF	  loc: 34.2289, -108.8677	  time: 2012-08-27T23:15:32	  back az [deg]:   108.97	  tr vel [m/s]: 362.4	  f-stat:  48.8	  high snr band: 0.17 - 5.40 Hz
        2D.W34CE..CDF	  loc: 34.2877, -110.0220	  time: 2012-08-27T23:22:23	  back az [deg]:   100.72	  tr vel [m/s]: 364.5	  f-stat:  12.5	  high snr band: 0.60 - 1.01 Hz
        2D.W42CE..CDF	  loc: 33.1100, -110.9228	  time: 2012-08-27T23:24:29	  back az [deg]:    73.53	  tr vel [m/s]: 368.1	  f-stat:  18.8	  high snr band: 0.42 - 0.79 Hz
        2D.W46CE..CDF	  loc: 33.0742, -111.3568	  time: 2012-08-27T23:28:04	  back az [deg]:    81.67	  tr vel [m/s]: 362.5	  f-stat:  60.0	  high snr band: 0.29 - 3.34 Hz
        2D.W49CE..CDF	  loc: 32.9786, -111.7163	  time: 2012-08-27T23:29:36	  back az [deg]:    82.26	  tr vel [m/s]: 357.8	  f-stat:  66.8	  high snr band: 0.18 - 3.30 Hz

        ====================
        Localization Summary
        ====================

        ##############
        ## index: 0 ##
        ##############
        parameters
        ----------
            back_az_width: 10.0
            range_max: 1000.0
            grid_resol: 180
            celerity_model: regional_lf

        result
        ------
            latitude: 33.588 deg +/- 17.29 km.
            longitude: -106.499 deg +/- 23.79 km.
            origin time: 2012-08-27T23:00:55.894000 +/- 91.4 s.

    From this information, a configuration for running Spectral Yield Estimation (SpYE) analysis can be defined by masking out the stations we wnat to exclude, defining the frequency band for spectral anlaysis, and other relevant quanties for the explosion model (e.g., ambient pressure, temperature, chemical or nuclear explosion model, doubling of yield for a surface burst vs. air burst).  The blast wave models used is the Friedlander blastwave with parameterization defined by Kinney & Graham scaling laws.  The configuration then looks like this,

    .. code:: none

        [YIELD]

        det_mask = 0, 0, 1, 1,  0, 0, 1, 1

        freq_min = 0.2
        freq_max = 1.0

        yld_min = 1
        yld_max = 1e3

        amb_press = 101.325
        amb_temp = 288.15
        grnd_burst = True
        exp_type = chemical

        tlm_label = data/TLMs/2007_08-

-   Spectral yield estimation is completed using the :code:`characterization` option in the event analysis methods.  As with localization, the event file is specified and configuration file used to relay the various parameters,

    .. code:: bash

        infrapy event characterize --ev-file data/HRR-5.ev.json.gz  --cnfg-file config/SpYE_HRR.cnfg

    As with other analysis methods, parameter information will be summarized and high level results:

    .. code:: none

        #####################################
        ##                                 ##
        ##             InfraPy             ##
        ##    Characterization Analysis    ##
        ##                                 ##
        #####################################


        Loading configuration info from: config/SpYE_HRR.cnfg

        Data summary:
          ev_file: data/HRR-5.ev.json.gz

        characterization parameters:
          det_mask: [0, 0, 1, 1, 0, 0, 1, 1]
          loc_index: 0
          tlm_label: data/TLMs/2007_08-
          freq_min: 0.2
          freq_max: 1.0
          yld_min: 1.0
          yld_max: 1000.0
          resolution: 200
          ref_rng: 1.0
          amb_press: 101.325
          amb_temp: 288.15
          grnd_burst: True
          exp_type: chemical

        =================
        Detection Summary
        =================
        2D.W22CW..CDF	  loc: 34.3109, -108.6895	  time: 2012-08-27T23:14:56	  back az [deg]:   109.19	  tr vel [m/s]: 358.2	  f-stat: 380.4	  high snr band: 0.17 - 7.00 Hz
        2D.W24CE..CDF	  loc: 34.2289, -108.8677	  time: 2012-08-27T23:15:32	  back az [deg]:   108.97	  tr vel [m/s]: 362.4	  f-stat:  48.8	  high snr band: 0.17 - 5.40 Hz
        2D.W46CE..CDF	  loc: 33.0742, -111.3568	  time: 2012-08-27T23:28:04	  back az [deg]:    81.67	  tr vel [m/s]: 362.5	  f-stat:  60.0	  high snr band: 0.29 - 3.34 Hz
        2D.W49CE..CDF	  loc: 32.9786, -111.7163	  time: 2012-08-27T23:29:36	  back az [deg]:    82.26	  tr vel [m/s]: 357.8	  f-stat:  66.8	  high snr band: 0.18 - 3.30 Hz

        Loading transmission loss statistics...
          2007_08-0.025Hz.tlm.json.gz
          2007_08-0.030Hz.tlm.json.gz
          2007_08-0.037Hz.tlm.json.gz
          2007_08-0.044Hz.tlm.json.gz
          2007_08-0.054Hz.tlm.json.gz
          2007_08-0.065Hz.tlm.json.gz
          2007_08-0.079Hz.tlm.json.gz
          2007_08-0.096Hz.tlm.json.gz
          2007_08-0.116Hz.tlm.json.gz
          2007_08-0.141Hz.tlm.json.gz
          2007_08-0.170Hz.tlm.json.gz
          2007_08-0.206Hz.tlm.json.gz
          2007_08-0.250Hz.tlm.json.gz
          2007_08-0.303Hz.tlm.json.gz
          2007_08-0.367Hz.tlm.json.gz
          2007_08-0.445Hz.tlm.json.gz
          2007_08-0.539Hz.tlm.json.gz
          2007_08-0.653Hz.tlm.json.gz
          2007_08-0.791Hz.tlm.json.gz
          2007_08-0.958Hz.tlm.json.gz
          2007_08-1.160Hz.tlm.json.gz
          2007_08-1.406Hz.tlm.json.gz
          2007_08-1.703Hz.tlm.json.gz
          2007_08-2.064Hz.tlm.json.gz
          2007_08-2.500Hz.tlm.json.gz

        Estimating yield using spectral amplitudes...

        Results Summary (tons eq. TNT):
            Maximum a Posteriori Yield: 45.52935074866948
            68% Confidence Bounds: [22. 81.]
            95% Confidence Bounds: [  4. 197.]

    The location from the BISL analysis is used to estimate the source-receiver distances and azimuths for evaluating the transmission loss models.  For comparison, the ground truth explosive yield of HRR-5 was 45.3 tons equivalent TNT.  Visualization of the SpYE analysis result can be done by referencing the event file, localization index, and characterization index,

    .. code:: bash

        infrapy plot event --ev-file data/HRR-5.ev.json.gz --loc-index 0 --char-index 0

    In this case, the localization and characterization results are shown along with the various waveforms sorted by distance from the location estimate.  Currently, all included detections are visualized; though, ongoing work on InfraPy will consider methods of denoting which detections were excluded from a localization or characterization analysis.

    .. image:: _static/_images/HRR5_char.png
        :width: 1200px
        :align: center
