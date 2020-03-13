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


## Installation
Infrapy requires anaconda to install.  The command below will create an environment named infrapy_env, install the necessary
packages into it, and install infrapy.  Navigate to the base directory of the infrapy 
package (there will be a file there named infrapy_env.yml), and run:

    >> conda env create -f infrapy_env.yml

If this command executes correctly and finishes without errors, it should print out instructions on 
how to activate and deactivate the new environment:

    To activate the environment, use:

        >> conda activate infrapy_env

    To deactivate an active environment, use

        >> conda deactivate

