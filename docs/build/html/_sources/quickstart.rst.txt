.. _quickstart:

=====================================
Quickstart
=====================================

A series of scripts illustrating how to use infrapy subroutines as stand-alone modules are found in the /examples folder.
The jupyter notebook documenting these steps is found in /tutorials/Quick Start.ipynb.  The notebook can be run by installing jupyter notebook via conda.

.. code-block:: python

    >> conda install jupyter notebook

Beamforming:

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
