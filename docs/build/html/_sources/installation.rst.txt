.. _installation:

=====================================
Installation
=====================================

-------------------------------------
Operating Systems
-------------------------------------

Infrapy can currently be installed on machines running newer versions of Linux or Apple OSX.  A Windows-compatible version is in development.

-------------------------------------
Anaconda
-------------------------------------

The installation of infrapy currently depends on Anaconda to resolve and download the correct python libraries. So if you don't currently have anaconda installed
on your system, please do that first.

Anaconda can be downloaded from https://www.anaconda.com/distribution/. Either 3.x or 2.x will work since the numbers refer to the Python version of the default
environment.  Infrapy's installation will create a new environment and will install the version of Python that it needs into that environment.

-------------------------------------
Infrapy Installation
-------------------------------------

Once Anaconda is installed, you can install infrapy by navigating to the base directory of the infrapy package (there will be a file there
named infrapy_env.yml), and run:

.. code-block:: bash

    >> conda env create -f infrapy_env.yml

If this command executes correctly and finishes without errors, it should print out instructions on how to activate and deactivate the new environment:

To activate the environment, use:

.. code-block:: none

    >> conda activate infrapy_env

To deactivate an active environment, use

.. code-block:: none

    >> conda deactivate

-------------------------------------
Testing
-------------------------------------

Once the installation is complete, you can test that the InfraPy methods are set up and accessible by first activating the environment with:

.. code-block:: none

    >> conda activate infrapy_env

The InfraPy command line methods have usage summarizes that can be displayed via the :code:`--help` option.  On the command line, run:

.. code-block:: none

    infrapy --help

The usage information should be displayed:

.. code-block:: none

    Usage: infrapy [OPTIONS] COMMAND [ARGS]...

      infrapy - Python-based Infrasound Data Analysis Toolkit

      Command line interface (CLI) for running and visualizing infrasound analysis

    Options:
      -h, --help  Show this message and exit.

    Commands:
      dets              Visualize infrapy analysis results
      run_assoc         Associate detections into events
      run_fd            Identify detections from beamforming results
      run_fk            Run beamforming methods on waveform data
      run_fkd           Run beamforming and detection methods in sequence
      run_loc           Estimate source locations and times for events
      utils             Various utility functions for infrapy analysis

Each of the individual methods have usage information (e.g., :code:`infrapy run_fk --help`) that will be discussed in the :ref:`quickstart`

