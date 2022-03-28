.. _localization:

============
Localization
============

* Event analysis uses the Bayesian Infrasonic Source Localization (BISL) methodology to estimate both the spatial location of the event as well as the origin time with quantified uncertainty
* Preliminary analysis of the back projections can be used to define the spatial region of interest or it can be specified by the analysis
* Analysis identifies the maximum posteriori solution
* The marginalized spatial distribution is approximated as 2d-normal to define 95 and 99% confidence ellipse bounds
* The marginalized temporal distribution is analyzed to identify the exact 95 and 99% confidence bounds
* Likelihood definitions relating detection parameters to spatial and temporal source characteristics are shared between association and localization analysis for consistency

See the following for more references on BISL:
`Modrak et al., 2010 <https://academic.oup.com/gji/article/181/1/399/718964>`_
`Marcillo et al., 2013 <https://academic.oup.com/gji/article/196/1/375/586767>`_
`Blom et al., 2015 <https://academic.oup.com/gji/article/203/3/1682/2594791>`_
