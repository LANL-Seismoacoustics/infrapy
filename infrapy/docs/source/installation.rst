.. _installation:

=====================================
Installation
=====================================

------------
Dependencies
------------

**Operating Systems**

Infrapy can currently be installed on machines running newer versions of Linux or Apple OS X.  A Windows-compatible version is in development.

**Python Environment Management**

The installation of infrapy currently depends on `Anaconda <https://docs.anaconda.com/free/anaconda/install/index.html>`_, `Miniconda <https://www.anaconda.com/docs/getting-started/miniconda/system-requirements>`_, `Miniforge <https://github.com/conda-forge/miniforge>`_ or some similar Python package and environment manager to resolve and download the correct python libraries. So if you don't currently have anaconda installed
on your system, please do that first.  Infrapy's installation will create a new environment and will install the version of Python that it needs into that environment.

-----------------------------------------
Environment Construction and Installation 
-----------------------------------------

**Infrapy Installation**

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

**Python Geophysics Suite (PyGS) Installation**


Infrasound software tools developed by LANL SMEs have become increasing coupled in usage so that having them in a common Python environment is useful.  An in-development Python Geophysics Suite (PyGS) YML file is included in the InfraPy repository that will build an environment and install InfraPy, infraGA/GeoAc, and stochprop from GitHub.  It can be run using the same syntax as above,

.. code-block:: bash

    >> conda env create -f pygs_env.yml

All dependencies will be installed and the LANL Python libraries pulled from GitHub to complete the environment.  To finish setting up, activate the environment and compile the infraGA/GeoAc software,

.. code-block:: bash

    >> conda activate pygs
    >> infraga compile 

**Python Geophysics Suite (PyGS) Installation - Dev Version**


Because the PyGS YML file installs via GitHub cloning, it doesn't copy the examples/ directories from the various libraries for demonstration and also doesn't leave the source code easily accessible for any de-bugging or customization.  A separate developer version is also included that requires a few more steps.  Build an instance of the environment with just InfraPy included using the included YML file,

.. code-block:: bash

    >> conda env create -f pygs-dev_env.yml

Next, clone the other repositories if you don't have them,

.. code-block:: bash

    >> git clone https://github.com/LANL-Seismoacoustics/infraga.git
    >> git clone https://github.com/LANL-Seismoacoustics/stochprop.git

If you have SSH keys set up for GitHub, you can alternately clone as,

.. code-block:: bash
	
    >> git clone git@github.com:LANL-Seismoacoustics/infraga.git
    >> git clone git@github.com:LANL-Seismoacoustics/stochprop.git

Once the PyGS development environment is built, activate it and use pip with the :code:`-e` flag to install infraGA/GeoAc and stochprop.  As with the non-dev install, compile the infraGA/GeoAc ray tracing methods,

.. code-block:: bash

    >> conda activate pygs_dev

    >> cd /path/to/stochprop
    >> pip install -e .

    >> cd /path/to/infraga
    >> pip install -e .

    >> infraga compile 


This installation will clone the example directories with all relevant data and also allow you to interact with other :code:`git` branches to access any in-development features and capabilities.

-------
Testing
-------

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
      detect  Detect signatures in infrasound data
      doc     Open infrapy manual
      event   Build and analyse events
      plot    Visualize infrapy analysis results
      utils   Various utility functions for infrapy analysis

Each of the individual methods have usage information (e.g., :code:`infrapy detect beam --help`) that will be discussed in the :ref:`quickstart`

