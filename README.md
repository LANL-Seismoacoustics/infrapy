# Infrapy

Infrapy is a tool for processing infrasound and seismic array data. It
implements a database-centric approach for pipeline continuous near real-time
analysis. The pipeline includes analysis at station and network levels (using
beam-forming and clustering techniques, respectively) for the detection,
association and location of events.  The pipeline relies on the interaction of
the algorithms with a relational database structure to organize and store
waveform data, the parameters for the analysis, and results of both levels of
analysis. Our implementation can interact seamlessly with traditional (e.g.:
Oracle) and serverless (e.g.: SQLite) relational databases.


## Authorship
Infrapy was built upon previous similar (InfraMonitor) tools and
developed by the LANL Seismoacoustics (LANL-SA) Team.  

## Operating Systems

Infrapy can currently be installed on machines running newer versions of Linux or Apple OSX. A Windows-compatible version is in development. 

Additionally, it is assumed that you will have internet access for the installation.  If you don't then please contact us directly and we will help you get this installed.

## Anaconda

The installation of infrapy currently depends on Anaconda to resolve and download the correct python libraries. So if you don’t currently have anaconda installed on your system, please do that first.

Anaconda can be downloaded from https://www.anaconda.com/distribution/. Either 3.x or 2.x will work since the numbers refer to the Python version of the default environment. Infrapy’s installation will create a new environment and will install the version of Python that it needs into that environment.

## Installation

Once Anaconda is installed, the command below will create an environment named infrapy_env, install the necessary packages into it, and install infrapy.  Navigate to the base directory of the infrapy package (there will be a file there named infrapy_env.yml), and run:

    >> conda env create -f infrapy_env.yml

If this command executes correctly and finishes without errors, it should print out instructions on 
how to activate and deactivate the new environment:

    To activate the environment, use:

        >> conda activate infrapy_env

    To deactivate an active environment, use

        >> conda deactivate
 
## Testing

Once the installation is complete, you can test some things by first activating the environment with:

    >> conda activate infrapy_env

Then navigate to the /example directory located in the infrapy base directory, and run the test scripts via something like:

    >> python test_beamforming.py

If infrapy was successfully installed, all of the test scripts should run and finish without any errors.

## Supplemental Data

Some of the example scripts included with infrapy in the /scripts directory depend on some supplemental data.  This data can be found at [https://github.com/LANL-Seismoacoustics/infrapy-data](https://github.com/LANL-Seismoacoustics/infrapy-data)

Instructions found there will guide you in its installation.

## Infraview

We supply a GUI application to help with quick data, beamforming, and location analysis. Once installation is complete, and the new environment is activated, you can run the GUI with the command:

    >> infraview


## Errors/issues

If you have any errors or issues related to the installation or basic functionality, the best way to get them to us is by submitting a new issue in the Issues Tab above. 


