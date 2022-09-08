.. _yield:

===========================
Source Characterization
===========================

Source characteristization (e.g., yield estimation) can be performed using InfraPy's Spectral Yield Estimation (SpYE) methods.  Detections are analyzed to compute the spectral amplitude (magnitude of the Fourier transform) and compare with a reference noise sample at each array defined as pre- or post-detection (InfraPy can also use the beamed signal and residuals, but that's an in-development method and might produce unexpected results). Transmission loss statistics built external to InfraPy are utilized to estimate the infrasonic spectral amplitude at a location near the source from the various observations.  Finally, the near-source spectral amplitude estimate can be combined with a source model (currently the Kinney & Graham scaling laws and Friedlander blastwave model) to map the spectral amplitude to blastwave yield

- An example set of transmission loss models for summer in the western US are included in :code:`infrapy/propagation/priors/tloss/` for use in the example analysis for the Humming Roadrunner 5 explosion.

- Construction of models requires use of the `NCPAprop software <https://github.com/chetzer-ncpa/ncpaprop-release>`_ and methods in :code:`infrapy/propagation/infrasound.py` (a separate Python library named :code:`stochprop` is being developed that will provide more straight forward means of creating these models)


See the following for more references on yield estimation using regional infrasound:

- `Friedlander, 1946 <https://doi.org/10.1098/rspa.1946.0046>`_

- `Kinney & Graham, 1985 <https://doi.org/10.1007/978-3-642-86682-1>`_

- `Blom et al., 2018 <https://doi.org/10.1093/gji/ggy258>`_
