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

- The blue and red windows in the waveform viewer can be moved and re-sized to identify some reference/noise segment (red) and an analysis segment of interest (blue).  The power spectral density (PSD) for each segment is shown in the upper right panel in red and blue, respectively. 

- Re-size the blue window to include all three high amplitude packets in the waveform set and set the red reference window either before or after to example the spectral content of the high amplitude arrivals.  Select the :code:`Apply Filter?` option on the right and adjust the gray window in the PSD to cover the frequencies where the blue line is notably above the red.  Then, click on the :code:`<-- Set Filter to -->` button to automatically adjust the filter frequencies.

    .. image:: _static/_images/infraview/wvfrms3.png
        :width: 1200px
        :align: center

- Once the analysis window and frequency range of interest are identified, switch to the Beamforming tab on the top of the GUI.

    .. image:: _static/_images/infraview/beam1.png
        :width: 1200px
        :align: center

- Run beamforming and identify detections...

    .. image:: _static/_images/infraview/beam2.png
        :width: 1200px
        :align: center

- Name and save detections...

    .. image:: _static/_images/infraview/detection1.png
        :width: 300px
        :align: center

- Download data from IRIS...using the station browser...

    .. image:: _static/_images/infraview/station_browser.png
        :width: 300px
        :align: center

- Loading data from IRIS...

    .. image:: _static/_images/infraview/fdsn_loader.png
        :width: 300px
        :align: center


- Continue analysis...

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

