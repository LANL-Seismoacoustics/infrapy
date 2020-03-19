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

## Downloading

In a terminal, navigate to a directory that you would like to put infrapy in, then download the repository by either https:

    >> git clone  https://github.com/LANL-Seismoacoustics/infrapy.git
    
or by ssh:

    >> git clone git@github.com:LANL-Seismoacoustics/infrapy.git
    
This will create a folder named infrapy. This will be the base directory for your installation.

## Installation

With Anaconda installed and the repository cloned, you can now install infrapy. The command below will create an environment named infrapy_env, install the necessary packages into it, and install infrapy into that environment.  Navigate to the base directory of infrapy (there will be a file there named infrapy_env.yml), and run:

    >> conda env create -f infrapy_env.yml

If this command executes correctly and finishes without errors, it should print out instructions on 
how to activate and deactivate the new environment:

    To activate the environment, use:

        >> conda activate infrapy_env

    To deactivate an active environment, use

        >> conda deactivate
        
## Updating

Infrapy is in continued development.  Features are added, bugs are fixed, and documentation is improved fairly continuously. It's good practice to pull the latest updates on a regular basis.  To do this in a terminal, simply navigate into the infrapy directory and run:

    >> git pull
 
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

We supply a GUI application to help with quick data, beamforming, and location analysis. Once installation is complete, you can activate the infrapy_env and run the GUI with the commands:
    
    >> conda activate infrapy_env
    >> infraview

## Errors/issues

If you have any errors or issues related to the installation or basic functionality, the best way to get them to us is by submitting a new issue in the Issues Tab above. 

Questions and problems that might not rise to the level of an Issue can be directed to:
  
jwebster@lanl.gov (Installation and GUI questions)

fransiska@lanl.gov (Documentation, Tutorials, and Scripting questions)

pblom@lanl.gov (Algorithms and general science questions)
