.. _pipeline:

=======================================
Running Pipeline Processing in Infrapy
=======================================

The folder tutorials/cli contains all necessary data and configuration files to begin utilizing the pipeline processing methodologies in Infrapy.

Once installed, the steps to run pipeline processing in Infrapy are:

1. Either load local waveform data into a sqlite database or connect to a Oracle database following instructions in :ref:`pisces`.

2. Create a configuration file. See :ref:`config` for a detailed description of the parameters that need to be included. Two example configuration files, one for connecting to the provided sqlite file and one for connecting to an oracle database are provided.

3. Run the FK analysis for a specific array:

.. code-block:: python

    >> infrapy run_fk --config_file BRPConfig.txt

4. Run the FD analysis for a specific array that has already FK results, remember to locate the parameter id from your FK analysis (you can use the script read_pfk.py to find the correct id).

.. code-block:: python

    >> infrapy run_fd --config_file BRPConfig.txt

  4. Once you have run detection on 2+ array, run association processing:

  .. code-block:: python

      >> infrapy run_assoc --config_file BRPConfig.txt




.. toctree::
    :maxdepth: 5
    :titlesonly:

    config
    fk
    fd
    assoc
    scripts
