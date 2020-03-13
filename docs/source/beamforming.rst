.. _beamforming:

===========================
Array Processing
===========================


The use of infrasonic arrays, specifically for CTBT applications, is preferable due to the inherent reduction in signal-to-noise ratios (SNR) originating from the summation of four or more recordings at each array. The nature of a decision-rule based detector requires data that has been pre-processed using beamforming methods. Beamforming, a form of array processing, is the first step in the analysis of data from infrasonic arrays.    Conventional beamforming methods (Bartlett, Capon) separate coherent and incoherent parts of a signal through the assumption of planar waves arriving at the array.  A signal backazimuth and slowness can be estimated as signals are shifted to account for travel time differentials across array elements, bringing the signal into phase across as the noise deconstructively cancels out.  In the classical, or Bartlett methodology , data records on each array element are time-shifted versions of the other with local noise,

See the following for more references on beamforming:
`Rost and Thomas 2002 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2000RG000100>`_
`Olson and Szuberla 2010 <https://link.springer.com/chapter/10.1007/978-0-387-30441-0_81>`_
`Costley 2013 <https://asa.scitation.org/doi/full/10.1121/1.4818940>`_
