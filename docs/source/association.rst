.. _association:

===========================
Association
===========================


InfraPy uses a pair-based, joint-likelihood method to identify those detection sets with high likelihood of originating from a common source.  Events are identified using a pair-based Bayesian algorithm that defines the association between detection pairs from their joint-likelihood and identifies events via hierarchical clustering analysis

- The current implementation utilizes only spatial and temporal coincidence, but additional detection information can further improve event identification

- Evaluation using a synthetic data set has shown some mixing of spatially similar events poorly resolved by network geometry and occasional inclusion of "noise" detections in event clustering


**Integration Region Analysis**

The joint-likelihood calculation used to quantify the probability that a pair of detections originated from a common source includes a pre-computation step to identify an integration region for the joint-likelihood.  The identification of this integration region can be complicated and is controlled by several parameters in InfraPy.  Consider the figure below showing a pair of detecting infrasound arrays (orange triangles).  The direction-of-arrival information from beamforming (fk) analysis is used to compute a back projection from each station out some maximum distance to identify intersections.  In InfraPy, this projection distance is controlled by :code:`range_max` parameter (:code:`--range-max` on the CLI).  In addition to this maximum projection range, the integration region analysis is dependent on the lighter blue projections that bound the central great circle path.  These bounding projections are defined by the :code:`back_az_width` parameter (:code:`--back-az-width` on the command line).  These values default to 2000 km and :math:`\pm 10^o` and are intended for regional analysis.  In the case of local or global scale analysis, the maximum range should be modified to be 10's of km or 1000's of km, respectively, depending on the likely source-receiver ranges.

    .. image:: _static/_images/beam-overlap.jpg
        :width: 400px
        :align: center

For a detailed summary of how the primary-primary intersection (green dot), primary-bounding intersections (dark-blue dots), and bounding-bounding intersections (light-blue dots) are used to identify the integration region center and radius as well as how special cases in which one array is within the beam of the other array, see Blom et al. (2020).

**Citation**

See the following for more references on infrasonic association:

- `Blom et al., 2020 <https://academic.oup.com/gji/advance-article-abstract/doi/10.1093/gji/ggaa105/5800992>`_





