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

Once the installation is complete, you can test some things by first activating the environment with:

.. code-block:: none

    >> conda activate infrapy_env

Then navigate to the /example directory located in the infrapy base directory, and run the test scripts via something like:

.. code-block:: none

    >> python test_beamforming.py

If infrapy was successfully installed, all of the test scripts should run and finish without any errors.

----------------------------------------
Running the InfraView GUI Application
----------------------------------------

Once installation is complete, and the new environment is activated, you can run the GUI with the command:

.. code-block:: none

    >> infraview
