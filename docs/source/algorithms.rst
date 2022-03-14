.. _algorithms:

=====================================
Algorithms
=====================================

- Analyst methods are modular so that results from other processing tools (e.g., Bloodhound, Cardinal, PMCC) can be used in later analysis steps and vice versa.
- Algorithms are written to be data source agnostic so analysis can be performed regardless of data source once IO method is understood
    - Waveform analysis utilizes an ObsPy Stream object for data ingestion 
    - Detection information is stored in .json or ascii .dat files that can be ingested using functions in :code:`infrapy.likelihoods`

***************************
Station Level Processing
***************************
_______________________________________
:ref:`beamforming`
_______________________________________
- Beamforming estimates parameters of coherent signals
- Capabilities include methods to characterize transients as well as persistent signals
- Transient signals are identified using standard Bartlett beaming
- Persistent signals can be investigated using Minimum Variance Distortionless Response (MVDR) or MUltiple Signal Classification (MUSIC) algorithms

____________________________________
:ref:`afd`
____________________________________
- Adaptive Fisher statistics determine when to declare a detection

***************************
Network Level Processing
***************************

_________________________________
:ref:`association`
_________________________________
- Events are identified using a pair-based Bayesian algorithm that defines the association between detection pairs from their joint-likelihood and identifies events via hierarchical clustering analysis
- Current implementation utilizes only spatial and temporal coincidence, but additional detection information can further improve event identification
- Evaluation using a synthetic data set shows some mixing of spatially similar events poorly resolved by network geometry and occasional inclusion of "noise" detections in event clustering

__________________________________
:ref:`localization`
__________________________________
- Event analysis using the Bayesian Infrasonic Source Localization (BISL) methodology to estimate both the spatial location of the event as well as the origin time with quantified uncertainty
- Preliminary analysis of the back projections can be used to define the spatial region of interest or it can be specified by the analysis
- Analysis identifies the maximum posteriori solution
- The marginalized spatial distribution is approximated as 2d-normal to define 95 and 99% confidence ellipse bounds
- The marginalized temporal distribution is analyzed to identify the exact 95 and 99% confidence bounds
- Likelihood definitions relating detection parameters to spatial and temporal source characteristics are shared between association and localization analysis for consistency


    .. toctree::
        :maxdepth: 5
        :hidden:

        beamforming
        detection
        association
        localization
        yield
