.. _config:

=====================================
Configuration Files
=====================================
We run the different steps of the processing by running a script along with a configuration file. This applies to both station and network analysis levels. Configuration files set the parameters that are required to perform a specific analysis. The configuration file for the array processing (station level analysis) has the following structure:

An example configuration file is provided in tutorials/.  Parameters for each field within the configuration file are outlined below.

.. code-block:: none

    # configuration file to run array (FK) processing and (FD) detection

    [database] # required
    url = sqlite:///example.sqlite
    # database for processing
    site = pisces.tables.css3:Site
    wfdisc = pisces.tables.css3:Wfdisc
    affiliation = pisces.tables.css3:Affiliation
    # schemas for tables utilized in processing
    # if processing in sqlite database, tables remain the same
    # if processing in Oracle database, tables should point to your custom global_.py file

    [GeneralParams]
    year=2012
    # year for processing
    dayofyearini=100
    # Julian Day to begin processing
    dayofyearend=102
    # Julian Day to stop processing
    station=BRP
    # REFSTA of station for processing
    channel=EDF
    # channel
    name=example
    # name of processing parameters
    cpucnt=30
    # number of cpu cores to use for processing
    domain=time
    # domain (time or frequency) to run FK and FD processing


    [FKParams]
    name=mid band fk test
    # name for fk processing
    freqmin=0.5
    # minimum frequency
    freqmax=5.0
    # maximum frequency
    beamwinlen=60
    # beam window length
    beamwinstep=30
    # beam window step
    backazmin=-180.0
    # minimum bz for processing
    backazmax=180.0
    # maximum bz for processing
    backazstep=1.5
    # bz step
    trvelmin=300.0
    # minimum trace velocity
    trvelmax=600.0
    # maximum trace velocity
    trvelstep=2.5
    # trace velocity step
    beammethod=bartlett
    # beam method
    fkresults=fk_res_brp
    # where fk processing results are stored
    numsources = 1
    func_fk = None

    [FDetectParams]
    back_az_lim=10
    # limit of bz deviation between consecutive fk results
    detwinlen=300.0
    # window length for adaptive f detection
    detthresh=0.99
    # detection threshold
    dsegmin=5
    detmethod=fstat
    tb_prod=4000
    adaptivewlen=1200
    #length of window for AFD
    pthreshold=.01
    #p-value for time domain detection
    pfkid=0
    # pkfid for FK results
    corrthreshold=0.5
    # threshold for correlation values
    mineventlength
    # minimum event length in seconds
    fkresults=fk_res_brp
    # fk results to run detection on
    fdresults=fd_res_example_brp
    # where detection results are saved



    [AssocLocParams]
    network=YJ
    # network for association
    pfdetectid=0
    # detection ID from detection processing (all arrays for assoc must have same detect ID)
    pfkid=2
    # beamforming ID (all arrays must have same FK ID)
    beamwidth=10.0
    rangemax=1000.0
    clusterthresh=4.0
    trimthresh=None
    eventdetmin=3
    # minimum # of detections to form event
    eventarrmin=2
    # minimum number of arrays for event
    duration = 60
    # duration (minutes) for association windows

    fdtable_1=fd_res_example_brp
    #fdtable_2=fd_res_fsu
    #fdtable_3=fd_res_wmu
    # tables where detection results are stored
    resultstable = test_assoc
    # table where association results will be stored
