"""
infrapy.detection.beamforming_new.py

Methods for reading in time series data, analyzing
it, and identifying infrasonic signals using various
beamforming methods.

Author            Philip Blom (pblom@lanl.gov)

"""
import warnings

import numpy as np

from numba import jit, float64, complex128

from scipy import signal
from scipy import stats
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, root

from pyproj import Geod

wgs84_proj = Geod(ellps='sphere')

# ####################### #
#    Data manipulation    #
# ####################### #
def stream_to_array_data(stream, latlon=None, t_start=None, t_end=None):
    """Extract time series from ObsPy stream on common time samples and define the array geometry

        Extracts the time series from individual traces of an Obspy stream and identifies a
        common set of time samples where all are defined.  Interpolates the individual traces
        into a single numpy array (x) for which x[m] = x_m(t).  The geometry of the array is
        also extracted to enable beamforming analysis.

        Parameters
        ----------
        stream : ObsPy stream
            Obspy stream containing traces for all array elements
        latlon : 2darray
            (M x 2) 2darray containing the latitudes and longitudes of the array elements if they aren't in the stream

        Returns:
        ----------
        x : 2darray
            M x N matrix of array data, x[m][n] = x_m(t_n)
        t : 1darray
            Vector of N sampled points in time, t[n] = t_n
        t_ref : datetime64
            Datetime corresponding to t[0]
        dxdy : 2darray
            M x 2 matrix of slowness vectors

        """

    # define common time samples
    t0 = max(np.datetime64(tr.stats.starttime) + np.timedelta64(int(tr.times()[0] * 1e3), 'ms') for tr in stream)
    t1 = min(np.datetime64(tr.stats.endtime) + np.timedelta64(int(tr.times()[0] * 1e3), 'ms') for tr in stream)

    dt = max(1.0 / tr.stats.sampling_rate for tr in stream)
    t = np.arange(0.0, (t1 - t0).astype('m8[ms]').astype(float) * 1.0e-3, dt)

    # interpolate each channel and evaluate on common time samples
    x = np.empty((len(stream), len(t)))
    for m, tr in enumerate(stream):
        t_ref = (np.datetime64(tr.stats.starttime) - t0).astype('m8[ms]').astype(float) * 1.0e-3
        temp = interp1d(t_ref + tr.times(), signal.detrend(tr.data), kind='linear')
        x[m] = temp(t)

    # if start/end times are given, apply mask
    if t_start and t_end:
        mask = np.logical_and(t_start <= t, t <= t_end)
        t = t[mask]
        x = x[:, mask]

    # define the geometry of the array
    dxdy = np.zeros((len(stream), 2))
    if latlon is None:
        for m, tr in enumerate(stream):
            temp = wgs84_proj.inv(tr.stats.sac['stlo'], tr.stats.sac['stla'], stream[0].stats.sac['stlo'], stream[0].stats.sac['stla'])
            dxdy[m] = np.array((temp[2] * np.sin(np.radians(temp[0])), temp[2] * np.cos(np.radians(temp[0]))))
        # Or using 'Coordinate' from stats
        # ...
    else:
        for m in range(1, len(stream)):
            temp = wgs84_proj.inv(latlon[m][1], latlon[m][0], latlon[0][1], latlon[0][0])
            dxdy[m] = np.array((temp[2] * np.sin(np.radians(temp[0])), temp[2] * np.cos(np.radians(temp[0]))))

    return x, t, t0, dxdy


def fft_array_data(x, t, window=None, sub_window_len=None, sub_window_overlap=0.5, fft_window="hanning", normalize_windowing=False):
    """Compute the Fourier transform of the array data to perform analysis

        Compute the Fourier transform of the array data within an analysis window defined by window = [t1, t2]
        and potentially using subwindows to obtain a full rank covariance matrix for beamforming analyses
        requiring such data.  Multiple FFT window options are available and a normalization option scales to
        account for the amplitude loss at the window edges.

        Parameters
        ----------
        x : 2darray
            M x N matrix of array data, x[m][n] = x_m(t_n)
        t : 1darray
            Vector of N sampled points in time, t[n] = t_n
        window : float
            Start and end time of the window relative to times in t, [t_1, t_2]
        sub_window_len : float
            Duration of the subwindow in seconds
        sub_window_overlap : float
            Fraction of subwindow to overlap (limited range of 0.0 to 0.9)
        fft_window : str
            Fourier windowing method
        normalize_windowing : boolean
            Boolean to apply normalization of window scaling

        Returns:
        ----------
        X : 2darray
            M x N_f matrix of the FFT'd data, X[m][n] = X_m(f_n)
        S : 3darray
            M x M x N_f cube of the covariance matrices, S[m1][m2][n] = mean(X_{m1}(f_n) conj(X_{m2}(f_n)))
        f : 1darray
            Vector of N_f frequencies for the FFT'd data, f[n] = f_n

    """

    M, N = x.shape
    dt = t[1] - t[0]

    if window:
        mask = np.logical_and(window[0] <= t, t <= window[1]).astype(int)
        win_n1, win_N = mask.nonzero()[0][0], sum(mask)
    else:
        win_n1, win_N = 0, N

    if sub_window_len:
        if sub_window_overlap < 0.0 or sub_window_overlap > 0.9:
            msg = "Inappropriate value in subwindow overlap.  Value is expected to be a fraction of the window length between 0.0 and 0.9."
            raise ValueError(msg)

        sub_win_N = int(sub_window_len / dt)
        sub_win_step = sub_win_N * sub_window_overlap

        padded_N = 2**int(np.ceil(np.log2(sub_win_N)))
        N_f = int(padded_N / 2 + 1)

        f = (1.0 / dt) * (np.arange(float(N_f)) / padded_N)
        X = np.zeros((M, N_f), dtype=complex)
        S = np.zeros((M, M, N_f), dtype=complex)

        # add contributions to X(f) and S(f) for each subwindow
        window_cnt = 0
        for n in range(win_n1, win_n1 + win_N, int(sub_win_N * (1.0 - sub_window_overlap))):
            if n != win_n1 and n + sub_win_N > win_n1 + (win_N - 1):
                break

            temp = np.zeros((M, padded_N))
            temp[:, 0 : sub_win_N] = x[:, n : n + sub_win_N]
            if fft_window == "hanning":
                temp[:, 0 : sub_win_N] *= np.array([np.hanning(sub_win_N)] * M)
                if normalize_windowing:
                    temp /= np.mean(np.hanning(sub_win_N))
            elif fft_window == "bartlett":
                temp[:, 0 : sub_win_N] *= np.array([np.bartlett(sub_win_N)] * M)
                if normalize_windowing:
                    temp /= np.mean(np.bartlett(sub_win_N))
            elif fft_window == "blackman":
                temp[:, 0 : sub_win_N] *= np.array([np.blackman(sub_win_N)] * M)
                if normalize_windowing:
                    temp /= np.mean(np.blackman(sub_win_N))
            elif fft_window == "hamming":
                temp[:, 0 : sub_win_N] *= np.array([np.hamming(sub_win_N)] * M)
                if normalize_windowing:
                    temp /= np.mean(np.hamming(sub_win_N))
            elif fft_window == "tukey":
                temp[:, 0 : sub_win_N] *= np.array([signal.tukey(sub_win_N)] * M)
                if normalize_windowing:
                    temp /= np.mean(signal.tukey(sub_win_N))
            elif fft_window == "boxcar":
                temp[:, 0 : sub_win_N] *= 1.0
            else:
                msg = "Unrecognized method in fft_window.  Options are 'hanning', 'bartlett', 'blackman', 'hamming', 'tukey', or 'boxcar'."
                raise ValueError(msg)

            # fft the data and add their contribution to X and S
            fft = np.fft.rfft(temp, axis=1) * dt

            X += fft
            for nf in range(0, int(padded_N / 2 + 1)):
                S[:,:,nf] += np.outer(fft[:, nf], np.conj(fft[:, nf]))
            window_cnt += 1

        # scale by number of windows used
        X /= window_cnt
        S /= window_cnt

    else:
        padded_N = 2**int(np.ceil(np.log2(win_N)))
        N_f = int(padded_N / 2 + 1)

        f = (1.0 / dt) * (np.arange(float(N_f)) / padded_N)
        X = np.zeros((M, N_f), dtype=complex)
        S = np.zeros((M, M, N_f), dtype=complex)

        # window and zero pad the data
        temp = np.zeros((M, padded_N))
        temp[:, 0 : win_N] = x[:, win_n1 : win_n1 + win_N]
        if fft_window == "hanning":
            temp[:, 0 : win_N] *= np.array([np.hanning(win_N)] * M)
            if normalize_windowing:
                temp /= np.mean(np.hanning(win_N))
        elif fft_window == "bartlett":
            temp[:, 0 : win_N] *= np.array([np.bartlett(win_N)] * M)
            if normalize_windowing:
                temp /= np.mean(np.bartlett(win_N))
        elif fft_window == "blackman":
            temp[:, 0 : win_N] *= np.array([np.blackman(win_N)] * M)
            if normalize_windowing:
                temp /= np.mean(np.blackman(win_N))
        elif fft_window == "hamming":
            temp[:, 0 : win_N] *= np.array([np.hamming(win_N)] * M)
            if normalize_windowing:
                temp /= np.mean(np.hamming(win_N))
        elif fft_window == "tukey":
            temp[:, 0 : win_N] *= np.array([signal.tukey(win_N)] * M)
            if normalize_windowing:
                temp /= np.mean(signal.tukey(win_N))
        elif fft_window == "boxcar":
            temp[:, 0 : win_N] *= 1.0
        else:
            msg = "Unrecognized method in fft_window.  Options are 'hanning', 'bartlett', 'blackman', 'hamming', or 'boxcar'."
            temp[:, 0 : win_N] *= 1.0
            raise(ValueError)

        # fft the data and define X(f) and S(f)
        X = np.fft.rfft(temp, axis=1) * dt
        for nf in range(0, int(padded_N / 2) + 1):
            S[:,:,nf] = np.outer(X[:, nf], np.conj(X[:, nf]))

    return X, S, f


# ############################# #
#     Slowness and delays for   #
#  defining the steering vector #
# ############################# #
def build_slowness(back_azs, trc_vels):
    """Compute the slowness values for a polar grid

        Computes the slowness grid usingg a polar grid defined by a series of back azimuth
        values and trave velocity values.  Returns a grid specified such that grid[n] is
        the x and y component of the nth slowness vector.

        Parameters
        ----------
        back_azs : 1darray
            Back azimuth values for slowness grid, K_1 values
        trc_vels : 1darray
            Trace velocity values for slowness grid K_2 values

        Returns:
        ----------
        grid : 2darray
            (K_1 x K_2) by 2 array of slowness vectors
    """

    back_az_grid, trc_vel_grid = np.meshgrid(back_azs, trc_vels)
    slowness_grid = np.array([np.sin(np.radians(back_az_grid.flatten())) / trc_vel_grid.flatten(),
                              np.cos(np.radians(back_az_grid.flatten())) / trc_vel_grid.flatten()]).T

    return slowness_grid


def compute_delays(dxdy, param_grid, param_opt='planar', sph_vel=340.0, sph_src_ht=0.0):
    """Compute the delays for a planewave

        Computes the time delays for each pair in param_grid given the
        array geometry in dxdy.  For planar parameterization, grid
        specifies s_x and s_y of the slowness.  For spherical
        parameterization, it specifies the x,y location of the source
        and requires specification of the velocity of the wavefront.

        For the slowness grid, use the build_slowness function to
        convert back azimuth and trace velocity values into a grid.
        Use np.meshgrid and flatten to produce a grid for the
        spherical wavefront source grid.

        Parameters
        ----------
        dxdy : 2darray
            M x 2 matrix describing the array geometry
        param_grid : 2darray
            K x 2 matrix of parameterization vectors containing
            either the slowness components (for 'planar') or
            the source location (for 'spherical')

        Returns:
        ----------
        delays : 2darray
            K x M of time delays across the array for each slowness
    """

    delays = np.zeros((dxdy.shape[0], param_grid.shape[0]))
    for m in range(dxdy.shape[0]):
        if param_opt == 'planar':
            delays[m] = np.array([np.dot(param_grid[k], dxdy[m]) for k in range(param_grid.shape[0])])
        else:
            delays[m] = np.array([np.sqrt(np.linalg.norm(param_grid[k] - dxdy[m])**2 + sph_src_ht**2) / sph_vel for k in range(param_grid.shape[0])])

    return delays.T


# ########################## #
#     Linear algbera for     #
#   beampower calculations   #
# ########################## #

@jit(complex128[:](complex128[:,:], complex128[:]), nopython=True)
def project_Ab(A, b):
    """Project matrix of K vectors, a_k, onto a vector b

        Projects a vector, b, of length M onto a set of K vectors, a_k,
        each of length M producing a vector, c, of length K

        Parameters
        ----------
        A : 2darray
            K x M matrix representing a set of K vectors, a_k, each of length M
        b : 1darray
            Vector of length M to project onto each a_k

        Returns:
        ----------
        c : 1darray
            Vector c where each scalar c_k = a_k^\dagger b
    """

    K, M = A.shape

    result_real = np.zeros(K)
    result_imag = np.zeros(K)
    for k in range(K):
        for m in range(M):
            temp = np.conj(A[k][m]) * b[m]
            result_real[k] += temp.real
            result_imag[k] += temp.imag

    return result_real + 1.0j * result_imag


@jit(float64[:](complex128[:,:], complex128[:,:]), nopython=True)
def project_ABA(A, B):
    """
    Project matrix of K vectors, a_k, onto Hermitian matrix B

        Projects each of K vectors, a_k, in matrix, A, of dimension K x M
        onto a Hermitian matrix, B, of dimension M x M producing a
        vector of scalars, c, of length K

        Parameters
        ----------
        A : 2darray
            K x M matrix representing a set of K vectors, a_k, each of length M
        B : 2darray
            M x M Hermitian matrix

        Returns:
        ----------
        c : 1darray
            Vector c where each scalar c_k = a_k^\dagger B a_k
    """
    K, M = A.shape

    result = np.zeros(K)
    for k in range(K):
        for m1 in range(M):
            for m2 in range(M):
                result[k] += (np.conj(A[k][m1]) * B[m1][m2] * A[k][m2]).real

    return result


@jit(complex128[:](complex128[:,:], complex128[:,:], complex128[:]), nopython=True)
def project_ABc(A, B, c):
    """Project matrix of K vectors, a_k, through Hermitian matrix B and onto vector c

        Projects each of K vectors, a_k, in matrix, A, of dimension K x M
        through a Hermitian matrix, B, of dimension M x M and onto a vector,
        c, producing a vector of scalars, d, of length K

        Parameters
        ----------
        A : 2darray
            K x M matrix representing a set of K vectors, a_k, each of length M
        B : 2darray
            M x M Hermitian matrix
        c : 1darray
            Vector of length M

        Returns:
        ----------
        d : 1darray
            Vector d where each scalar d_k = a_k^\dagger B c
    """
    K, M = A.shape

    result_real = np.zeros(K)
    result_imag = np.zeros(K)
    for k in range(K):
        for m1 in range(M):
            for m2 in range(M):
                temp = np.conj(A[k][m1]) * B[m1][m2] * c[m2]
                result_real[k] += temp.real
                result_imag[k] += temp.imag

    return result_real + 1.0j * result_imag


# ####################### #
#           Run           #
#       Beamforming       #
# ####################### #
def compute_beam_power(data, steering, method="bartlett", ns_covar_inv=None, signal_cnt=1):
    """Compute the beampower for a specific frequency

        Cmoputes the beampower at a single frequency using either the FFT'd data, X(f),
        for Bartlett or GLS analysis or the covariance matrix, S(f), for Capon and MUSIC.

        Generalized Least Square (GLS) analysis requires a noise covariance for the
        background which must be M x M where M is the length of X(f).

        MUltiple SIgnal Classification (MUSIC) analysis requires knowledge of the number
        of coherent signals in the data specified as signal_cnt.

        Parameters
        ----------
        data : ndarray
            Vector of length M, X_m(f_n), for "bartlett" and "gls" or matrix of
            dimension M x M, S(f_n), for covariance based methods
        steering : 2darray
            Matrix representing K steering vectors each of length K
        method : str
            Beamforming method to be applied to the data (must match
            for of data)
        ns_covar_inv : 2darray
            Noise covariance used in "gls" beamforming method
        signal_cnt : int
            Number of signals assumed in MUSIC algorithm

        Returns:
        ----------
        beam_power : 1darray
            Beam power for each steering vector (length K)
        """


    # note: for Bartlett (bartlett_covar), Capon and MUSIC, "data" contains S(f), while
    # for Bartlett and GLS it contains X(f)
    if method == "bartlett" or method == "gls":
        if len(data.shape) != 1:
            msg = "Incompatible beamforming method and data format: {} requires data vector X(f).".format(method)
            warnings.warn(msg)
    elif method == "bartlett_covar" or method == "capon" or method == "music":
        if len(data.shape) != 2:
            msg = "Incompatible beamforming method and data format: {} requires covariance matrix S(f).".format(method)
            warnings.warn(msg)
    else:
        msg = "Invalid beamforming method: {}.".format(method)
        warnings.warn(msg)

    beam_power = np.empty(len(steering))
    if method == "bartlett":
            temp = project_Ab(steering, data)
            beam_power = (np.conj(temp) * temp).real

    elif method == "gls":
        if ns_covar_inv is None:
            # Note: generalized least squares with noise covariance of identity
            # is equivalent to Bartlett beam.
            temp = project_Ab(steering, data)
            beam_power = (np.conj(temp) * temp).real

            # beam_power = gls_beam(data, steering, np.eye(data.shape[0], dtype=np.complex))
        else:
            num = project_ABc(steering, ns_covar_inv, data)
            den = project_ABA(steering, ns_covar_inv)
            beam_power = (np.conj(num) * num).real / den**2

    elif method == "bartlett_covar":
        beam_power = project_ABA(steering, data)

    elif method == "capon":
        temp = data + 1.0e-3 * np.mean(np.diag(data)) * np.eye(data.shape[0])
        covariance_inverse = np.linalg.inv(temp)
        beam_power = 1.0 / project_ABA(steering, covariance_inverse)

    elif method == "music":
        temp = data + 1.0e-3 * np.mean(np.diag(data)) * np.eye(data.shape[0])
        _, eigenvectors = np.linalg.eigh(temp)
        eigenvectors = eigenvectors.T

        noise_subspace = np.dot(eigenvectors[:-signal_cnt].T, np.conj(eigenvectors[:-signal_cnt]))
        beam_power = 1.0 / project_ABA(steering, noise_subspace)

    else:
        msg = "Invalid beamforming method: {}.".format(method)
        warnings.warn(msg)
        beam_power = None

    return beam_power


def compute_beam_power_wrapper(args):
    return compute_beam_power(*args)


def run(X, S, f, dxdy, delays, freq_band, method="bartlett", ns_covar_inv=None, signal_cnt=1, normalize_beam=True, pool=None):
    """Run beamforming analysis over frequencies of interest

        Computes the beam at multiple frequencies within a specified band given data in X(f)
        and S(f) and frequencies f as produced by the fft_array_data function.

        Normalization of the beam returns coherence in the case of Bartlett and a normalized
        version of the Capon beam but does not alter the output of the MUSIC algorithm as
        its result is a mathematical projection onto a noise subspace.

        A multiprocessing pool can be used to accelerate calculation of different frequencies
        in parallel.

        Parameters
        ----------
        X : 2darray
            M x N_f matrix of the FFT'd data, X[m][n] = X_m(f_n)
        S : 3darray
            M x M x N_f cube of the covariance matrices, S[m1][m2][n] = mean(X_{m1}(f_n) conj(X_{m2}(f_n)))
        f : 1darray
            Frequencies
        delays : 1darray
            Set of delays for the parameterization (length K)
        freq_band : iterable
            List or tuple with minimum and maximum frequency (e.g.,  [f_min, f_max])
        method : str
            Beamforming method to be applied to the data (must match form of data)
        signal_cnt : int
            Number of signals assumed in MUSIC algorithm
        ns_covar_inv : 2darray
            Noise covariance used in "gls" beamforming method
        normalize_beam : boolean
            Option to normalize the beam and return coherence (value between 0 and 1)
        pool : multiprocessing pool
            Multiprocessing pool for accelerating caluclation (maps over frequency)
        param_opt : string
            Option for the solution parameterization: 'planar' or 'spherical'
        sph_vel : float
            Velocity of the wavefront in the 'spherical' param_opt method

        Returns:
        ----------
        bmpwr : 2darray
            Beam power for each steering vector at each frequency in the band (dimension K x N_f)
        """

    band_mask = np.logical_and(freq_band[0] <= f, f <= freq_band[1])
    X_msk = X[:, band_mask]
    S_msk = S[:, :, band_mask]
    f_msk = f[band_mask]

    f_cnt = f_msk.shape[0]
    if pool:
        if method == "bartlett_covar" or method == "capon" or method == "music":
            args = [(S_msk[:, :, nf], np.exp(-2.0j * np.pi * f_msk[nf] * delays) / np.sqrt(X_msk.shape[0]), method, None, signal_cnt) for nf in range(f_cnt)]
        else:
            if ns_covar_inv is not None:
                args = [(X_msk[:, nf], np.exp(-2.0j * np.pi * f_msk[nf] * delays) / np.sqrt(X_msk.shape[0]), method, (ns_covar_inv[:, :, band_mask])[:, :, nf], signal_cnt) for nf in range(f_cnt)]
            else:
                args = [(X_msk[:, nf], np.exp(-2.0j * np.pi * f_msk[nf] * delays) / np.sqrt(X_msk.shape[0]), method, None, signal_cnt) for nf in range(f_cnt)]
        beam_power = np.array(pool.map(compute_beam_power_wrapper, args))

    else:
        if method == "bartlett_covar" or method == "capon" or method == "music":
            beam_power = np.array([compute_beam_power(S_msk[:, :, nf], np.exp(-2.0j * np.pi * f_msk[nf] * delays) / np.sqrt(X_msk.shape[0]), method, None, signal_cnt) for nf in range(f_cnt)])
        else:
            if ns_covar_inv is not None:
                beam_power = np.array([compute_beam_power(X_msk[:, nf], np.exp(-2.0j * np.pi * f_msk[nf] * delays) / np.sqrt(X_msk.shape[0]), method, (ns_covar_inv[:, :, band_mask])[:, :, nf], signal_cnt) for nf in range(f_cnt)])
            else:
                beam_power = np.array([compute_beam_power(X_msk[:, nf], np.exp(-2.0j * np.pi * f_msk[nf] * delays) / np.sqrt(X_msk.shape[0]), method, None, signal_cnt) for nf in range(f_cnt)])

    if normalize_beam:
        if method == "bartlett" or method == "gls" or method == "bartlett_covar":
            beam_power = np.array([beam_power[nf] / (np.vdot(X_msk[:, nf], X_msk[:, nf])).real for nf in range(f_cnt)])
        elif method == "capon":
            beam_power = np.array([beam_power[nf] / (np.max(np.linalg.eigh(S_msk[:, :, nf])[0])) for nf in range(f_cnt)])

    return beam_power


# ####################### #
#         Analyze         #
#    Beamforming Result   #
# ####################### #
def pure_state_filter(S):
    """Compute the pure state filter applied to a Hermitian matrix, S(f)

        Computes the pure state filter for a matrix.  Here, the covariance matrix is utilized
        to measure the average coeherence across the entire array at a given frequency.

        Pure state filter value are useful for weighting a multi-frequency beam average.

        Parameters
        ----------
        S : 3darray
            Covariance matrix of data in analysis window for all frequencies,
            x(t) --> S(f) = mean(X(f) X^\dagger(f))

        Returns:
        ----------
        pure_state : 1darray
            Pure state filter value at each frequency
        """

    # convert to coherence matrix
    M = S.shape[0]
    coh = np.empty_like(S)
    for i in range(M):
        for j in range(M):
            coh[i][j] = S[i][j] / np.sqrt(abs(S[i][i]) * abs(S[j][j]))

    return (M * np.trace(np.matmul(coh, coh)) - np.trace(coh)**2) / ((M - 1) * np.trace(coh)**2)


def find_peaks(beam_power, slowness_vals1, slowness_vals2, signal_cnt=1, freq_weights=None):
    """Identify the peak(s) in the beampower defined over a slowness grid

        Finds the peaks of a distribution using a frequency averaged beamforming result
        over a defined slowness grid.

        Parameters
        ----------
        beam_power : 2darray
            Beam power for each steering vector at each frequency in the
            band (dimension K x N_f)
        slowness_vals1 : 1darray
            Slowness values along first axis (polar or Cartesian grid)
        slowness_vals2 : 1darray
            Slowness values along second axis (polar or Cartesian grid)
        signal_cnt : int
            Number of signals to identify in the slowness grid
        freq_weights : string or 1darray
            Weights or method to use in frequency averaging of the beam power

        Returns:
        ----------
        peaks : ndarray
            signal_cnt x 3 array of the peaks identified containing slowness
            value 1 (back azimuth), slowness value 2 (trace velocity),
            and beam value
        """

    # Average over frequency and reshape
    if freq_weights == "doa_proj":
        bm_cnt = beam_power.shape[0]

        slowness = build_slowness(slowness_vals1, slowness_vals2) # might need to simplify this part
        pk_slows = np.array([slowness[np.argmax(beam_power[j])] for j in range(bm_cnt)])
        c_ref = min(slowness_vals2)

        doa_x = np.array([c_ref * pk_slows[j][0] for j in range(bm_cnt)])
        doa_y = np.array([c_ref * pk_slows[j][1] for j in range(bm_cnt)])
        doa_z = np.array([np.sqrt(1.0 - 0.99 * c_ref**2 * (pk_slows[j][0]**2 + pk_slows[j][1]**2)) for j in range(bm_cnt)])

        doa_proj = doa_x * np.mean(doa_x) + doa_y * np.mean(doa_y) + doa_z * np.mean(doa_z)
        doa_std = np.sqrt(np.var(doa_x) + np.var(doa_y) + np.var(doa_z))

        wts = (1.0 + doa_proj) / 2.0 + doa_std
        avg_beam = np.average(beam_power, axis=0, weights=wts)
    else:
        avg_beam = np.average(beam_power, axis=0, weights=freq_weights)

    avg_beam = avg_beam.reshape(len(slowness_vals2), len(slowness_vals1))

    # Determine the number of peaks on the grid and compare with
    # the specified signal count
    peaks = []
    if signal_cnt == 1:
        x = np.argwhere(avg_beam == avg_beam.max())
        m = x[0][0]
        n = x[0][1]

        n_up, n_dn = min(n + 1, len(slowness_vals1) - 1), max(n - 1, 0)
        m_up, m_dn = min(m + 1, len(slowness_vals2) - 1), max(m - 1, 0)

        dPds1 = (avg_beam[m][n_up] - avg_beam[m][n_dn]) / (slowness_vals1[n_up] - slowness_vals1[n_dn])
        dPds2 = (avg_beam[m_up][n] - avg_beam[m_dn][n]) / (slowness_vals2[m_up] - slowness_vals2[m_dn])

        ddPds1s1 = (avg_beam[m][n_up] - 2.0 * avg_beam[m][n] + avg_beam[m][n_dn]) / ((slowness_vals1[n_up] - slowness_vals1[n_dn]) / 2.0)**2
        ddPds2s2 = (avg_beam[m_up][n] - 2.0 * avg_beam[m][n] + avg_beam[m_dn][n]) / ((slowness_vals2[m_up] - slowness_vals2[m_dn]) / 2.0)**2
        
        ddPds1s2 = (avg_beam[m_up][n_up] - avg_beam[m_up][n_dn] - avg_beam[m_dn][n_up] + avg_beam[m_dn][n_dn])
        ddPds1s2 = ddPds1s2 / ((slowness_vals1[n_up] - slowness_vals1[n_dn]) * (slowness_vals2[m_up] - slowness_vals2[m_dn]))
        
        if ddPds1s1 * ddPds2s2 - ddPds1s2**2 > 0.0:
            ds1 = - (ddPds2s2 * dPds1 - dPds2 * ddPds1s2) / (ddPds1s1 * ddPds2s2 - ddPds1s2**2)
            ds2 = - (ddPds1s1 * dPds2 - dPds1 * ddPds1s2) / (ddPds1s1 * ddPds2s2 - ddPds1s2**2)        
            dP = dPds1 * ds1 + dPds2 * ds2 + (ddPds1s1 / 2.0) * ds1**2 + (ddPds2s2 / 2.0) * ds2**2 + ddPds1s2 * ds1 * ds2
        else:
            ds1 = 0.0
            ds2 = 0.0
            dP = 0.0

        peaks.append([slowness_vals1[n] + ds1, slowness_vals2[m] + ds2, avg_beam[m][n] + dP, n, m])        
    else :
        for n in range(1, len(avg_beam[0, :-1])):
            if np.max(avg_beam[:, n - 1]) <= np.max(avg_beam[:, n]) >= np.max(avg_beam[:, n + 1]):
                m = np.argmax(avg_beam[1:-1, n])

                n_up, n_dn = min(n + 1, len(slowness_vals1) - 1), max(n - 1, 0)
                m_up, m_dn = min(m + 1, len(slowness_vals2) - 1), max(m - 1, 0)

                dPds1 = (avg_beam[m][n_up] - avg_beam[m][n_dn]) / (slowness_vals1[n_up] - slowness_vals1[n_dn])
                dPds2 = (avg_beam[m_up][n] - avg_beam[m_dn][n]) / (slowness_vals2[m_up] - slowness_vals2[m_dn])

                ddPds1s1 = (avg_beam[m][n_up] - 2.0 * avg_beam[m][n] + avg_beam[m][n_dn]) / ((slowness_vals1[n_up] - slowness_vals1[n_dn]) / 2.0)**2
                ddPds2s2 = (avg_beam[m_up][n] - 2.0 * avg_beam[m][n] + avg_beam[m_dn][n]) / ((slowness_vals2[m_up] - slowness_vals2[m_dn]) / 2.0)**2
        
                ddPds1s2 = (avg_beam[m_up][n_up] - avg_beam[m_up][n_dn] - avg_beam[m_dn][n_up] + avg_beam[m_dn][n_dn])
                ddPds1s2 = ddPds1s2 / ((slowness_vals1[n_up] - slowness_vals1[n_dn]) * (slowness_vals2[m_up] - slowness_vals2[m_dn]))
        
                if ddPds1s1 * ddPds2s2 - ddPds1s2**2 > 0.0:
                    ds1 = - (ddPds2s2 * dPds1 - dPds2 * ddPds1s2) / (ddPds1s1 * ddPds2s2 - ddPds1s2**2)
                    ds2 = - (ddPds1s1 * dPds2 - dPds1 * ddPds1s2) / (ddPds1s1 * ddPds2s2 - ddPds1s2**2)        
                    dP = dPds1 * ds1 + dPds2 * ds2 + (ddPds1s1 / 2.0) * ds1**2 + (ddPds2s2 / 2.0) * ds2**2 + ddPds1s2 * ds1 * ds2
                else:
                    ds1 = 0.0
                    ds2 = 0.0
                    dP = 0.0

                peaks.append([slowness_vals1[n] + ds1, slowness_vals2[m] + ds2, avg_beam[m][n] + dP, n, m]) 

    peaks = np.array(peaks)
    sorting = peaks[:, 2].argsort()[::-1]
    peaks = peaks[sorting]

    if len(peaks) < signal_cnt:
        print("WARNING. Only found " + str(len(peaks) + " local maxima in the grid."))
    peaks = peaks[:signal_cnt]

    return peaks[:, :3]


def project_beam(beam_power, back_az_vals, trc_vel_vals, freq_weights=None, method="max"):
    """Project polar slowness grid onto only azimuth

        Projects the polar slowness grid onto back azimuth and trace velocity in order to
        more easily view each.  The method can either use the maximum value to project or
        average to approximate the marginal distribution

        Parameters
        ----------
        beam_power : 2darray
            Beam power for each steering vector at each frequency in the band (K x N_f)
        back_az_vals : 1darray
            Back azimuth values defining polar slowness grid
        trc_vel_vals : 1darray
            Trace velocity values defining polar slowness grid
        freq_weights : 1darray
            Weights to use in frequency averaging of the beam power
        method : str
            Determines whether mean or maximum along trace velocity axis is used to
            define the projections

        Returns:
        ----------
        back_az_proj : 1darray
            Projection of the beam power onto the back azimuth axis
        trc_vel_proj : 1darray
            Projection of the beam power onto the trace velocity axis
        """

    # Average over frequency and reshape
    avg_beam = np.average(beam_power, axis=0, weights=freq_weights)
    avg_beam = avg_beam.reshape(len(trc_vel_vals), len(back_az_vals))

    # Project by taking maximum along perpendicular axis
    if method == "max":
        back_az_proj = np.amax(avg_beam, axis=0)
        trc_vel_proj = np.amax(avg_beam, axis=1)
    elif method == "mean":
        back_az_proj = np.mean(avg_beam, axis=0)
        trc_vel_proj = np.mean(avg_beam, axis=1)
    else:
        msg = "Invalue method in beam projection.  Options are 'max' and 'mean'."
        raise ValueError(msg)

    return back_az_proj, trc_vel_proj

def extract_signal(X, f, slowness, dxdy):
    """Extract the signal along the beam for a given slowness vector

        Extract the "best beam" signal from the array data for a given slowness pair and
        array geometry.  Returns both the extracted signal and the residual on each trace
        of the array

        Note: following Laslo's work, the frequency domain Fisher ratio can be computed as:
            F[nf] = abs(sig_est)**2 / np.mean(np.abs(residual), axis=1)**2 * (X.shape[1] - 1)

        Parameters
        ----------
        f : 1darray
            Frequencies
        X : 2darray
            FFT of data in analysis window, x(t) --> X(f)
        slowness : 2darray
            Slowness components (either back azimuth and trace velocity or
            s_x and s_y depending on slowness option)
        dxdy : 2darray
            Array geometry

        Returns:
        ----------
        sig_estimate : 1darray
            Extracted frequency domain signal along the beam
        residual : 2darray
            Residual across the array once beamed signal is extracted
        """

    delays = (dxdy[:, 0] * np.sin(np.radians(slowness[0])) + dxdy[:, 1] * np.cos(np.radians(slowness[0]))) / slowness[1]

    sig_estimate = np.empty_like(X[0])
    residual = np.empty_like(X)
    for nf in range(len(f)):
        steering = np.exp(-2.0j * np.pi * f[nf] * delays)
        sig_estimate[nf] = np.vdot(steering, X[:, nf]) / np.vdot(steering, steering)
        residual[:, nf] = X[:, nf] - sig_estimate[nf] * steering

    return sig_estimate, residual



# ###################### #
#        Identify        #
#       Detections       #
# ###################### #
def calc_det_thresh(fstat_vals, det_thresh, TB_prod, channel_cnt):
    kde = stats.gaussian_kde(fstat_vals)
    def temp_kde(f):
        return -kde.pdf(f)

    def temp_fstat(f):
        return -stats.f(TB_prod, TB_prod * (channel_cnt - 1)).pdf(f)

    c_scalar = minimize_scalar(temp_kde).x / minimize_scalar(temp_fstat).x

    return stats.f(TB_prod, TB_prod * (channel_cnt - 1)).ppf(det_thresh) * c_scalar[0]

def detect_signals(times, beam_results, win_len, TB_prod, channel_cnt, det_thresh=0.99, min_seq=5, back_az_lim=15, fixed_thresh=None):
    """Identify detections with beamforming results

        Identify detection in the beamforming results using either Kernel Density
        Estimate (KDE) fits to the f-statistic distribution or the adaptive
        F-detector methods developed by Arrowsmith.

        Parameters
        ----------
        times : 1darray
            Times of beamforming results as numpy datetime64's
        beam_results : 2darray
            Beamforming results consisting of back azimuth, trace velocity, and
            f-value at each time step. This is a 2D array with dimensions (len(times), 3), 
            where the first column has back azimuth values, the second has trace velocity 
            values, and the third has f-statistic values
        win_len : float
            Window length to define the adaptive fstat threshold
        TB_prod : int
            Time-bandwidth product needed to compute the Fisher statistic
        channel_cnt : int
            Number of channels on the array; needed to compute the Fisher statistic
        det_thresh : float
            Threshold for declaring a detection
        min_seq : int
            Threshold for the number of sequential above-threshold values to declare
            a detection
        back_az_lim : float
            Threshold below which the maximum separation of back azimuths must be
            in order to declare a detection
        fixed_thresh : float
            A fixed detection threshold for fstat values (overrides adaptive 
                threshold calculation)

        Returns:
        ----------
        dets : list
            List of identified detections including detection time, relative start
            and end times of the detection, back azimuth, trace velocity, and f-stat.
        """

    back_az_vals = beam_results[:, 0]
    trc_vel_vals = beam_results[:, 1]
    fstat_vals = beam_results[:, 2]

    det_mask = np.zeros_like(fstat_vals, dtype=bool)

    print("Running detection analysis...")
    if fixed_thresh:
        det_mask = (fstat_vals > fixed_thresh)
    else:
        # check first full window for detections using KDE of all f-values
        thresh = calc_det_thresh(fstat_vals, det_thresh, TB_prod, channel_cnt)
        init_win_mask = (times - times[0]).astype('m8[s]').astype(float) < win_len
        det_mask[init_win_mask] = (fstat_vals[init_win_mask] > thresh)

        for n, tn in enumerate(times):
            if (tn - times[0]).astype('m8[s]').astype(float) > win_len:
                fstat_vals_win = fstat_vals[np.logical_and(tn - np.timedelta64(int(win_len), 's') <= times, times < tn)]
                thresh = calc_det_thresh(fstat_vals_win, det_thresh, TB_prod, channel_cnt)
                det_mask[n] = (fstat_vals[n] > thresh)

    # Check for detections shorter than the minimum sequence 
    #   length and with too large of back azimuth deviations
    n, dets = 0, []
    while n < (len(det_mask) - min_seq):
        if np.all(det_mask[n:n + min_seq]):
            det_len = min_seq
            while np.all(det_mask[n:n + det_len]) and n + det_len < len(det_mask):
                det_len += 1

            back_az_min = np.argmin(back_az_vals[n:n + det_len])
            back_az_max = np.argmax(back_az_vals[n:n + det_len])

            if abs(back_az_max - back_az_min) < back_az_lim or abs(back_az_max - back_az_min) - 360.0 < back_az_lim:
                pk_index = np.argmax(fstat_vals[n:n + det_len]) 

                if n == 0:
                    warnings.warn("Detection is close to start of analysis.  Detection start time is set to beginning of data, but this might be incorrect")
                if n + det_len == len(times):
                    warnings.warn("Detection is close to end of analysis. Detection end time is set to end of data, but this might be incorrect")

                try:
                    det_time = times[n + pk_index]
                    det_start = (times[n] - times[n + pk_index]).astype('m8[s]').astype(float)
                    det_end = (times[n + det_len - 1] - times[n + pk_index]).astype('m8[s]').astype(float)

                    back_az = back_az_vals[n + pk_index]
                    trc_vel = trc_vel_vals[n + pk_index]
                    fstat = fstat_vals[n + pk_index]
                    dets = dets + [[det_time, det_start, det_end, back_az, trc_vel, fstat]]
                except Exception as ex1:
                    print('Issue with detection time ' + str(det_time), ex1)

            n += det_len
        else:
            n += 1

    return dets
