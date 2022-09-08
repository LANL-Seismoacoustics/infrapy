.. _beamforming:

===========================
Array Processing
===========================


The use of infrasonic arrays, specifically for CTBT applications, is preferable due to the inherent reduction in signal-to-noise ratios (SNR) originating from the summation of four or more recordings at each array. The nature of a decision-rule based detector requires data that has been pre-processed using beamforming methods. Beamforming, a form of array processing, is the first step in the analysis of data from infrasonic arrays.    Conventional beamforming methods (Bartlett, Capon) separate coherent and incoherent parts of a signal through the assumption of planar waves arriving at the array.  A signal backazimuth and slowness can be estimated as signals are shifted to account for travel time differentials across array elements, bringing the signal into phase across as the noise deconstructively cancels out.  In the classical, or Bartlett methodology , data records on each array element are time-shifted versions of the other with local noise,


***************************
Mathematical Details
***************************

- Beamforming in InfraPy is completed in the frequency domain using a model in which the data on the :math:`m^\text{th}` sensor of an array is a combation of a time-shifted coherent signal and local noise,

    .. math::
        x_m \left( t \right) = g \left( t - \tau_m \right) + \eta_m \left( t \right)

    In vector form, this is,

    .. math::
        \vec{x} \left( t \right) = g \left( t - \vec{\tau} \right) + \vec{\eta} \left( t \right)

    Taking a Fourier transform of this relation gives the signal model form in the frequency domain,

    .. math::
        \vec{\mathcal{X}} \left( f \right) = \mathcal{G} \left( f \right) \vec{\Phi} + \vec{\mathcal{N}} \left( f \right), \quad \vec{\Phi} \left(f, \vec{\tau} \right) = e^{- 2 i \pi f \vec{\tau}}

    where :math:`\vec{\Phi}` is a vector of phasors often termed the "steering vector".

- The time delays in the steering vector can be defined for a plane wave assumption by characterizing the planewave by a direction-of-arrival azimuth, :math:`\varphi`, and a trace velocity at which the sound propagates across the ground surface, :math:`v_\text{tr} = \frac{c_0}{\cos \vartheta}`, that varies with the ambient sound speed, :math:`c_0`, and planewave inclination angle, :math:`\vartheta`.  These parameter define the slowness, :math:`\vec{s}` of the planewave, which can be combined with the locations of the array elements, :math:`\vec{z}_j`, to define the arrival time delays,

    .. math::
        \tau_j = \vec{s} \cdot \vec{z}_j, \quad \vec{s} = \left( \frac{\sin \varphi}{v_\text{tr}}, \frac{\cos \varphi}{v_\text{tr}} \right)

- The ordinary least squares estimate for the coherent signal is more commonly termed the Bartlett beam and can be defined as,

    .. math::
        \hat{\mathcal{G}}_\text{Bartlett} \left( f, \vec{s} \right) = \frac{\vec{\mathcal{X}}^\dagger \vec{\Phi}}{\vec{\Phi}^\dagger \vec{\Phi}}

    In the case that the background covariance is known, the generalized least squares (GLS) method can be applied,

    .. math::
        \hat{\mathcal{G}}_\text{GLS} \left( f, \vec{s}, \mathbf{S}_\mathcal{N} \right) = \frac{\vec{\mathcal{X}}^\dagger \mathbf{S}_\mathcal{N}^{-1} \vec{\Phi}}{\vec{\Phi}^\dagger \mathbf{S}_\mathcal{N}^{-1} \vec{\Phi}}

    The noise covariance can be estimated in the :ref:`infraview` interface via the red noise window and in the CLI :code:`run_fk` methods via the :code:`noise_start` and :code:`noise_end` parameters.  Applying the GLS method in an automated way requires adaptively estimating the noise covariance while analyzing data, which is an area of ongoing research.

    In both of these cases, the estimated beam power is often expressed as the square of the estimated signal amplitude (e.g., :math:`\mathcal{P}_\text{Bartlett} \left( f, \vec{s} \right) = \hat{\mathcal{G}}_\text{Bartlett}^2 \left( f, \vec{s} \right)`).


- Several covariance-based relations discussed in detail in Krim & Viberg (1996) are also implemented in the InfraPy beamforming methods.  For each method, a signal covariance matrix is computed from the outer product of the signal vector, :math:`\mathbf{S}_\mathcal{X} = \frac{1}{J} \sum_j{ \vec{\mathcal{X}}_j \vec{\mathcal{X}}_j ^\dagger}` and that matrix is used in analysis.  In order to utilize these methods, multiple sub-windows are needed to ensure that the estimated covariance matrix is full rank for inversion or eigen-decomposition; however, a whitening step is included to ensure these options are possible.  From the above definition, it's clear that the Bartlett beam can be defined from the signal covariance as,

    .. math::
        \hat{\mathcal{P}}_\text{Bartlett} \left( f, \vec{s} \right) = \left| \frac{\vec{\Phi}^\dagger \mathbf{S}_\mathcal{X} \vec{\Phi}}{\vec{\Phi}^\dagger \vec{\Phi}} \right|^\frac{1}{2}


    The Minimum Variance Distortionless Response (MVDR, also teremd the Capon beam) has the form,
  
    .. math::
        \hat{\mathcal{P}}_\text{MVDR} \left( f, \vec{s} \right) = \frac{1}{\vec{\Phi}^\dagger \mathbf{S}_\mathcal{X}^{-1} \vec{\Phi}}
        
    The MUltiple SIgnal Classification (MUSIC) algorithm takes an eigen-decomposition of the signal covariance and defines the noise sub-space from the :math:`M - q` eigenvectors associated with the lowest eigenvalues.  Denoting this subpace as :math:`\mathbf{\Pi}_q` where :math:`q` denotes the assumed number of signal eigenvectors, the MUSIC beam has a form similar to the MVDR,

    .. math::
        \hat{\mathcal{P}}_\text{MUSIC} \left( f, \vec{s}, q \right) = \frac{1}{\vec{\Phi}^\dagger \mathbf{\Pi}_q \vec{\Phi}}


See the following for more references on beamforming:

- `Rost and Thomas, 2002 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2000RG000100>`_

- `Olson and Szuberla, 2010 <https://link.springer.com/chapter/10.1007/978-0-387-30441-0_81>`_

- `Costley, 2013 <https://asa.scitation.org/doi/full/10.1121/1.4818940>`_

- `Krim and Viberg, 1996 <https://doi.org/10.1109/79.526899>`_