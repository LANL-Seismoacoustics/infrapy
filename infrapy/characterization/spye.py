# estimate_yield.py
#
# Functions to run the Bayesian Infrasonic Spectral Yield Estimation (BISYE) methods

import warnings

import numpy as np

import matplotlib.cm as cm
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline
from scipy.signal import savgol_filter

from ..detection import beamforming_new
from ..utils import prog_bar, confidence

from scipy.special import gamma

#########################
##   Kinney & Graham   ## 
##   Blastwave Model   ##
#########################
def kg_op(W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical"):
    """
        Kinney & Graham scaling law peak overpressure model
                
        Parameters
        ----------
        W: float
            Explosive yield of the source [kg eq. TNT]
        r: float
            Propagation distance [km]
        p_amb: float
            Ambient atmospheric pressure [kPa]
        T_amb: float
            Ambient atmospheric temperature [deg K]
        exp_type: string
            Type of explosion modeled, options are "chemical" or "nuclear"
        
        Returns
        -------
        p0: float
            Peak overpressure [Pa]    
    """
    
    fd = (p_amb / 101.325)**(1.0 / 3.0) * (T_amb / 288.15)**(1.0 / 3.0)
    sc_rng = fd / W**(1.0 / 3.0) * r * 1000.0
    
    if exp_type=="chemical":
        term1 = 1.0 + (sc_rng / 4.5)**2
        term2 = np.sqrt(1.0 + (sc_rng / 0.048)**2)
        term3 = np.sqrt(1.0 + (sc_rng / 0.32)**2)
        term4 = np.sqrt(1.0 + (sc_rng / 1.35)**2)
        
        result = 808.0 * term1 / (term2 * term3 * term4)
    else:
        term1 = (1.0 + sc_rng / 800.0)
        term2 = np.sqrt(1.0 + (sc_rng / 87.0)**2)
        
        result = 3.2e6 / sc_rng**3 * term1 * term2

    return p_amb * 1.0e3 * result


def kg_ppd(W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical"):
    """
        Kinney & Graham scaling law positive phase duration model
                
        Parameters
        ----------
        W : float
            Explosive yield of the source [kg eq. TNT]
        r : float
            Propagation distance [km]
        p_amb : float
            Ambient atmospheric pressure [kPa]
        T_amb : float
            Ambient atmospheric temperature [deg K]
        exp_type : string
            Type of explosion modeled, options are "chemical" or "nuclear"
            
        Returns
        -------
        t0 : float
            Positive phase duration [s]
    """

    fd = (p_amb / 101.325)**(1.0 / 3.0) * (T_amb / 288.15)**(1.0 / 3.0)
    sc_rng = fd / W**(1.0 / 3.0) * r * 1000.0
    
    if exp_type=="chemical":
        term1 = 1.0 + (sc_rng / 0.54)**10
        term2 = 1.0 + (sc_rng / 0.02)**3
        term3 = 1.0 + (sc_rng / 0.74)**6
        term4 = np.sqrt(1.0 + (sc_rng / 6.9)**2)
        
        result = 980.0 * term1 / (term2 * term3 * term4)
    else:
        term1 = np.sqrt(1.0 + (sc_rng / 100.0)**3)
        term2 = np.sqrt(1.0 + (sc_rng / 40.0))
        term3 = (1.0 + (sc_rng / 285.0)**5)**(1.0 / 6.0)
        term4 = (1.0 + (sc_rng / 50000.0))**(1.0 / 6.0)
        
        result = 180.0 * term1 / (term2 * term3 * term4)

    return result * W**(1.0 / 3.0) / 1e3


def blastwave(t, W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical", shaping_param=0.0):
    """ 
        Acoustic blastwave that can be used as a source
            for surface explosions out to several scale
            kilometers.
        
        Note: alpha = 0 produces the Friedlander
            blastwave model.
        
        Parameters
        ----------
        t: float
            Time [s]
        W: float
            Explosive yield of the source [kg eq. TNT]
        r: float
            Propagation distance [km]
        p_amb: float
            Ambient atmospheric pressure [kPa]
        T_amb: float
            Ambient atmospheric temperature [deg K]
        exp_type: string
            Type of explosion modeled, options are "chemical" or "nuclear"
        shaping_param: float
            Shaping parameter for waveform that controls high frequency trend (set to 0 for original Friedlander blastwave)

        Returns
        -------
        p : float
            Overpressure at time t            
    """

    p0 = kg_op(W, r, p_amb, T_amb, exp_type)
    t0 = kg_ppd(W, r, p_amb, T_amb, exp_type)
    
    x = t / t0 + (1.0 + shaping_param)
    x0 = (1.0 + shaping_param) - np.sqrt(1.0 + shaping_param)
    norm = x0**shaping_param * (1.0 - x0 / (1.0 + shaping_param)) * np.exp(-x0)
    
    result = np.empty_like(x)
    if len(np.atleast_1d(result)) == 1:
        if x >= 0.0:
            return p0 / norm * x**shaping_param * (1.0 - x / (1.0 + shaping_param)) * np.exp(-x)
        else:
            return 0.0
    else:
        result[x >= 0.0] = p0 / norm * x[x >= 0.0]**shaping_param * (1.0 - x[x >= 0.0] / (1.0 + shaping_param)) * np.exp(-x[x >= 0.0])
        result[x <  0.0] = 0.0
        return result


def blastwave_spectrum(f, W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical", shaping_param=0.0):
    """ 
        Fourier transform amplitude for the acoustic
            blastwave in sasm.acoustic.blastwave().
        
        Note: shaping_param = 0 produces the Friedlander
        blastwave model.
        
        Note: the peak of the spectrum occurs at
        f_0 = \frac{1}{2 \pi t_0} \frac{1}{\sqrt{\alpha + 1}
        and t0 corresponding to a given peak frequency is
        t_0 = \frac{1}{2 \pi f_0} \frac{1}{\sqrt{\alpha + 1}
        
        Parameters
        ----------
        f: float
            Frequency [Hz]
        W: float
            Explosive yield of the source [kg eq. TNT]
        r: float
            Propagation distance [km]
        p_amb: float
            Ambient atmospheric pressure [kPa]
        T_amb: float
            Ambient atmospheric temperature [deg K]
        exp_type: string
            Type of explosion modeled, options are "chemical" or "nuclear"
        shaping_param: float
            Shaping parameter for waveform that controls high frequency trend (set to 0 for original Friedlander blastwave)
        
        Returns
        -------
        P : float
            Spectral value at frequency f
    """

    p0 = kg_op(W, r, p_amb, T_amb, exp_type)
    t0 = kg_ppd(W, r, p_amb, T_amb, exp_type)

    omega = 2.0 * np.pi * f
    x0 = (1.0 + shaping_param) - np.sqrt(1.0 + shaping_param)
    norm = x0**shaping_param * (1.0 - x0 / (1.0 + shaping_param)) * np.exp(-x0)
    
    return p0 / norm * t0 * gamma(shaping_param + 1.0) * (omega * t0) / (1.0 + (omega**2 * t0**2))**(shaping_param / 2.0 + 1.0)


#################################
##        Spectral Yield       ##
##  Estimation (SpYE) Methods  ## 
#################################
def extract_spectra(det_list, st_list, win_buffer=0.25, ns_opt="pre"):
    """ 
        Extract spectra for a list of detections using a list of
            obspy streams with a defined window buffer and noise
            window option
        
        Parameters
        ----------
        det_list : :obj:`list` of :obj:`InfrasoundDetection`
            Iterable of detections with defined start and end times
        st_list : :obj:`list` of :obj:`obspy.Stream`
            Iterable of Obspy streams containing the array channels for each detection
        win_buffer : float
            Scaling factor defining the window buffer (a 20 second detection with a 
                buffer of 0.25 adds 5 seconds to the start and end of the window)
        ns_opt : string
            Option defining how the noise window is defined.  Options include "pre",
                "post", and "beam" to use the preceding, following, or beam residual

        Returns
        -------
        smn_spec : ndarray
            Spectral amplitude of the signal-minus-noise            
    """
    print("Computing detection spectra...")
    det_cnt = len(det_list)
    smn_spec = [0] * det_cnt

    for j in range(det_cnt):
        st = st_list[j]
        t_ref = (det_list[j].peakF_UTCtime - np.datetime64(st[j].stats.starttime)).astype('m8[s]').astype(float)
        det_len = det_list[j].end - det_list[j].start
        
        x, t, t0, geom = beamforming_new.stream_to_array_data(st)
        
        sig_t1 = t_ref + det_list[j].start - det_len * win_buffer
        sig_t2 = t_ref + det_list[j].end + det_len * win_buffer
        if ns_opt == "pre" or ns_opt == "post":
            X_sig, _, f_sig = beamforming_new.fft_array_data(x, t, window=[sig_t1, sig_t2], fft_window="tukey")
            X_sig = np.mean(abs(X_sig), axis=0)
            if ns_opt == "pre":
                ns_t1 = t_ref + det_list[j].start - det_len * win_buffer - det_len * (1.0 + 2.0 * win_buffer)
                ns_t2 = t_ref + det_list[j].start - det_len * win_buffer
            else:
                ns_t1 = t_ref + det_list[j].end + det_len * win_buffer
                ns_t2 = t_ref + det_list[j].end + det_len * win_buffer + det_len * (1.0 + 2.0 * win_buffer)
            X_ns, _, f_ns = beamforming_new.fft_array_data(x, t, window=[ns_t1, ns_t2], fft_window="tukey")
            X_ns = np.mean(abs(X_ns), axis=0)
    
        elif ns_opt == "beam":
            X_temp, _, f_sig = beamforming_new.fft_array_data(x, t, window=[sig_t1, sig_t2], fft_window="tukey")
            X_sig, residual = beamforming_new.extract_signal(X_temp, f_sig, np.array([det_list[j].back_azimuth, det_list[j].trace_velocity]), geom)
            
            X_sig = abs(X_sig)
            f_ns = f_sig
            X_ns = np.mean(abs(residual), axis=0)
        else:
            msg = "ERROR!  Invalid noise window option.  Options are 'pre', 'post', and 'beam'"
            raise ValueError(msg)

        smn_spec[j] = np.array([f_sig, savgol_filter(10.0 * np.log10(abs(X_sig - X_ns)), 49, 3)])

    return smn_spec


def run(det_list, smn_spec, src_loc, freq_band, tloss_models, resol=150, yld_rng=np.array([10.0, 10.0e3]), ref_src_rng=1.0, grnd_brst=True, p_amb=101.325, T_amb=288.15, exp_type="chemical"):
    """ 
        Run Spectral Yield Estimation (SpYE) methods to estimate explosive yield
        
        Parameters
        ----------
        det_list: :obj:`list` of :obj:`InfrasoundDetection`
            Iterable of detections with defined start and end times
        smn_spec: ndarray
            Spectral amplitude of the signal-minus-noise  
        freq_band: iterable
            List or tuple with minimum and maximum frequency (e.g.,  [f_min, f_max])
        tloss_models: list of frequencies and infrapy.propagation.TLossModel instances
            Propagation transmission loss model list with reference frequencies (see
            test/test_yield.py for construction example)
        resol: int
            Resolution used in analysis
        yld_rng: iterable
            List or tuple with minimum and maximum yields (e.g., [yld_min, yld_max])
        ref_src_rng: float
            Standoff distance used to define source model distance
        grnd_brst: boolean
            Boolean declaring whether the source is a ground burst (interaction with source
            doubles the effective yield for a ground burst vs. air burst)
        p_amb: float
            Ambient atmospheric pressure [kPa]
        T_amb: float
            Ambient atmospheric temperature [deg K]
        exp_type: string
            Type of explosion modeled, options are "chemical" or "nuclear"

        Returns
        -------
        yld_vals: 
            Values of explosive yield for the PDF                    
        yld_pdf: 
            Probability of the observed source having a given yield                    
        conf_bnds: ndarray
            Limits for the 68% and 95% confidence bounds on yield            
    """

    print("Estimating yield using spectral amplitudes...")
    det_cnt = len(det_list)
    
    freqs = np.logspace(np.log10(max(tloss_models[0][0], freq_band[0])), np.log10(min(tloss_models[0][-1], freq_band[1])), resol)
    pdf = np.empty((det_cnt, resol**2))
    
    obs_spec_ref = np.mean([max(smn_spec[j][1]) for j in range(det_cnt)])
    src_spec_vals = np.linspace(obs_spec_ref - 10.0, obs_spec_ref + 40.0, resol)

    # Compute the combined near-source spectral amplitude   
    for j in range(det_cnt):
        _, _, pdf[j] = det_list[j].src_spec_pdf(src_loc[0], src_loc[1], freqs, src_spec_vals, smn_spec[j], tloss_models)
        
    psd_fit = interp2d(freqs, src_spec_vals + 10.0 * np.log10(1.0 / ref_src_rng), np.product(pdf, axis=0).reshape((resol, resol)))
    yld_vals = np.logspace(np.log10(yld_rng[0]), np.log10(yld_rng[1]), resol)
    yld_pdf = np.empty_like(yld_vals)

    for n in range(len(yld_vals)):
        if grnd_brst:
            def temp(f):
                return psd_fit(f, 10.0 * np.log10(blastwave_spectrum(f, yld_vals[n] * 2.0, ref_src_rng, p_amb, T_amb, exp_type))) / f
        else:
            def temp(f):
                return psd_fit(f, 10.0 * np.log10(blastwave_spectrum(f, yld_vals[n] * 1.0, ref_src_rng, p_amb, T_amb, exp_type))) / f

        yld_pdf[n] = quad(temp, freqs[0], freqs[-1], limit=250, epsrel=5.0e-3)[0]

    yld_interp = interp1d(yld_vals, yld_pdf)

    conf_bnds = [0] * 2
    conf_bnds[0], _, temp = confidence.find_confidence(yld_interp, [yld_vals[0], yld_vals[-1]], 0.68)
    conf_bnds[1], _, temp = confidence.find_confidence(yld_interp, [yld_vals[0], yld_vals[-1]], 0.95)

    result = {'yld_vals': yld_vals, 'yld_pdf' : yld_pdf, 'conf_bnds': conf_bnds}

    return result 
