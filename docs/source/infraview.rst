.. _infraview:

=====================================
InfraView
=====================================

InfraPy includes an interactive graphical user interface (GUI) called InfraView.  This interface includes waveform analysis methods, beamforming and detection tools, and event identification and localization capabilities.  The InfraView GUI can be accessed from the command line by simply running,

    .. code-block:: bash

        infraview

--------------------
Array-Level Analyses
--------------------

- The newly opened InfraView GUI defaults to the waveform analysis tab as shown below. This part of the interface includes a waveform viewer (top left), data summary (bottom left), spectral viewer (top right), and filter controls (bottom right). 

    .. image:: _static/_images/infraview/wvfrms1.png
        :width: 900px
        :align: center

- Local waveform data can be ingested from :code:`File > Load Waveform(s)...`, and selecting the files in the system viewer.  Locate the infrapy installation directory and load the 4 :code:`examples/data/YJ.BRP*.SAC` files.  This will populate the data info in the bottom left window, visualize the waveform data in the upper left, and show the power spectral density in the upper right.  You can scroll to read through the trace and station info in the bottom left.

    .. image:: _static/_images/infraview/wvfrms2.png
        :width: 1200px
        :align: center

- A FDSN data downloader is supplied to enable one to retrieve data directly from IRIS or other FDSN compatible servers. It can be accessed from :code:`File > Import Waveform(s)...`

    .. image:: _static/_images/infraview/fdsn_loader.png
        :width: 300px
        :align: center

- The FDSN data downloader includes a station browser that allows you to search for stations by code, or can be used to identify stations that were active at a certain time, and/or within a certain range of a location.  This allow you to select appropriate stations and send them to the FDSN downloader window for easy viewing of its waveforms.

    .. image:: _static/_images/infraview/station_browser.png
        :width: 300px
        :align: center

- The blue and red windows in the waveform viewer can be moved and re-sized to identify some reference/noise segment (red) and an analysis segment of interest (blue).  The power spectral density (PSD) for each segment is shown in the upper right panel in red and blue, respectively. 

- Re-size the blue window to include all three high amplitude packets in the waveform set and set the red reference window either before or after to example the spectral content of the high amplitude arrivals.  Select the :code:`Apply Filter?` option on the right and adjust the gray window in the PSD to cover the frequencies where the blue line is notably above the red.  Then, click on the :code:`<-- Set Filter to -->` button to automatically adjust the filter frequencies.

    .. image:: _static/_images/infraview/wvfrms3.png
        :width: 1200px
        :align: center

- Once the analysis window and frequency range of interest are identified, switch to the Beamforming tab on the top of the GUI.

- Beamformer settings can be changed in the tab at the bottom of the window.  This allows you to set the beamforming algorithm, the window length and step, and the range and resolution of the back azimuth and trace velocities that are searched.

- Detector settings can also be set in the tabs at the bottom.  If the "Automatically calculate threshold" checkbox is checked, then the threshold will be calculated from the red window previously selected in the Waveforms analysis tab. The threshold can also be set manually.  The Back azimuth limit is the maximum spread in the back azimuth for which a detection can have. For a static explosion, this can be small, for a moving source, this could be set to 360 degrees.  The minimum peak with is the number of continuous points that must be above the threshold for something to be considered a detection.


    .. image:: _static/_images/infraview/beam1.png
        :width: 1200px
        :align: center

- Once you have the beamforming and detector settings that you want, you can click on "Run Beamforming" near the top.  The threshold will be calculated, then the beamformer will run.  You will see the window move through the waveform showing which part of the wave is being analyzed. Values for the F statistic, the back azimuth, and the trace velocity for each window will be plotted for each window step.

    .. image:: _static/_images/infraview/beam2.png
        :width: 1200px
        :align: center

-  When the beamforming run is complete, if detections were found, a window will pop up giving the analyst a chance to name the detection, add an event, and add a note to the detections.  Also, the analyst can discard a detection if not needed. When the Finish button is pressed, the detection(s) will be added to the table in the Detections tab, and will appear on the F-statistic plot with a grey box outlining the start and end of the detections.  These can be clicked on and moved around by the analyst if the auto-detector was deemed inaccurate.

    .. image:: _static/_images/infraview/detection1.png
        :width: 300px
        :align: center


- Once detections are recorded for a given array, clear out the waveforms, and repeat the process for another array...

    .. image:: _static/_images/infraview/fdsn_data.png
        :width: 1200px
        :align: center



----------------------
Network-Level Analyses
----------------------

- Load an example detection set for event ID...

    .. image:: _static/_images/infraview/load_detections.png
        :width: 1200px
        :align: center

- Load an example detection set for event ID...

    .. image:: _static/_images/infraview/event_id1.png
        :width: 1200px
        :align: center

- Select an event in the distance matrix...

    .. image:: _static/_images/infraview/event_id2.png
        :width: 1200px
        :align: center

- Run localization on the event

    .. image:: _static/_images/infraview/location1.png
        :width: 1200px
        :align: center

