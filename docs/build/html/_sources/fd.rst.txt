.. _fd:

===================================
Infrapy Detection (FD) Processing
===================================

Infrapy detects signals using an Adaptive F-Detector (AFD). The adaptive F-detector was developed by Arrowsmith et al., (2009) to account for both correlated and uncorrelated noise through modification of the conventional F-statistic. The detector accounts for temporal changes in noise by applying an adaptive window to update the detection distribution, which allows for the distinction between signal and correlated noise.

.. image:: _static/_images/AFD.png


__________________
Configuration file
__________________
Detection is run using the same configuration file as FK processing.  Detection requires input from FK processing results

.. code-block:: python

    >> infrapy run_fd --config_file BRPConfig.txt
