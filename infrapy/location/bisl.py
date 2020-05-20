# infrapy.location.bisl.py
#
# The bisl methods utilized to analyze the posterior pdf
# obtained by combining detection likelihood pdfs.  The
# methods here use global coordiantes and assume two
# spatial parameters, latitude and longitude, along with
# the source time, t.
#
# Author            Philip Blom (pblom@lanl.gov)

import sys
import datetime
import time
import itertools
import random
import copyreg

import numpy as np

from scipy.integrate import quad, nquad
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import minimize, bisect
from scipy.special import i0
from scipy.stats import chi2

from pyproj import Geod

from ..propagation import likelihoods as lklhds
from ..utils import latlon as ll

####################################
#### Set Integration Parameters ####
####    and spherical globe     ####
####################################
int_opts = {'limit': 100, 'epsrel': 1.0e-3}
sph_proj = Geod(ellps='sphere')

def set_region(det_list, bm_width=10.0, rng_max=np.pi / 2.0 * 6370.0, rad_min=100.0, rad_max=1000.0):
    """Defines the integration region for computation of the BISL probability distribution

        Projects finite width beams from each of the detecting arrays and looks for intersections
        of the primary (center) and secondary (edge) lines to define the integration region
        for computation of the localization distribution

        Parameters
        ----------
        det_list : iterable of InfrasoundDetection instances
            Detections attributed to the event
        bm_width : float
            Width of the projected beam [degrees]
        rng_max : float
            Maximmum range for beam projection [km]
        rad_min : float
            Minimum radius of the integration region [km]
        rad_max : float
            Maximum radius of the integration region [km]

        Returns:
        ----------
        Center : float
            Center of the integration region as latitude, longitude pair [degrees]
        Radius : float
            Radius of the integration region [km]

        """

    det_cnt = len(det_list)
    pair_cnt = det_cnt**2

    main_intersects = np.empty((pair_cnt, 2,))
    beam_intersects = np.empty((8 * pair_cnt, 2,))

    main_intersects[:] = np.nan
    beam_intersects[:] = np.nan

    for n1 in range(det_cnt):

        latlon1 = np.array([det_list[n1].latitude, det_list[n1].longitude], dtype=np.float)

        proj1 = ll.sphericalfwd(latlon1, 90, det_list[n1].back_azimuth)[0][0]
        proj1up = ll.sphericalfwd(latlon1, 90, det_list[n1].back_azimuth + bm_width)[0][0]
        proj1dn = ll.sphericalfwd(latlon1, 90, det_list[n1].back_azimuth - bm_width)[0][0]

        for n2 in range(det_cnt):
            pair_index = n1 * det_cnt + n2

            latlon2 = np.array([det_list[n2].latitude, det_list[n2].longitude], dtype=np.float)

            if not np.array_equal(latlon1,latlon2):
                proj2 = ll.sphericalfwd(latlon2, 90, det_list[n2].back_azimuth)[0][0]
                proj2up = ll.sphericalfwd(latlon2, 90, det_list[n2].back_azimuth + bm_width)[0][0]
                proj2dn = ll.sphericalfwd(latlon2, 90, det_list[n2].back_azimuth - bm_width)[0][0]

                main_intersects[pair_index] = ll.gcarc_intersect(latlon1, proj1, latlon2, proj2)[0]

                beam_intersects[pair_index * 8 + 0] = ll.gcarc_intersect(latlon1, proj1, latlon2, proj2up)[0]
                beam_intersects[pair_index * 8 + 1] = ll.gcarc_intersect(latlon1, proj1, latlon2, proj2dn)[0]
                beam_intersects[pair_index * 8 + 2] = ll.gcarc_intersect(latlon1, proj1up, latlon2, proj2)[0]
                beam_intersects[pair_index * 8 + 3] = ll.gcarc_intersect(latlon1, proj1dn, latlon2, proj2)[0]
                beam_intersects[pair_index * 8 + 4] = ll.gcarc_intersect(latlon1, proj1up, latlon2, proj2up)[0]
                beam_intersects[pair_index * 8 + 5] = ll.gcarc_intersect(latlon1, proj1dn, latlon2, proj2dn)[0]
                beam_intersects[pair_index * 8 + 6] = ll.gcarc_intersect(latlon1, proj1up, latlon2, proj2dn)[0]
                beam_intersects[pair_index * 8 + 7] = ll.gcarc_intersect(latlon1, proj1dn, latlon2, proj2up)[0]

    all_intersects = np.vstack((main_intersects, beam_intersects))

    # Compute geographic mean of the intersections to define the center of the region
    x, y, z = 0.0, 0.0, 0.0
    for n in range (9):
        if not np.isnan(all_intersects[n][0]):
            x += np.sin(np.radians(all_intersects[n][0])) * np.cos(np.radians(all_intersects[n][1]))
            y += np.sin(np.radians(all_intersects[n][0])) * np.sin(np.radians(all_intersects[n][1]))
            z += np.cos(np.radians(all_intersects[n][0]))

    norm = np.sqrt(x**2 + y**2 + z**2)
    x /= norm
    y /= norm
    z /= norm
    center = [np.degrees(np.arccos(z)), np.degrees(np.arctan2(y, x))]

    # Compute the distance to the non-central intersections to estimate the size of the region
    # Limit the size of the possible region (diameter of 100 - 1,000 km, can't include a detecting station)
    rngs = [np.nan] * len(all_intersects)
    for n in range (len(main_intersects),len(all_intersects)):
        if not np.isnan(all_intersects[n][0]): rngs[n] = sph_proj.inv(center[1], center[0], all_intersects[n][1], all_intersects[n][0], radians=False)[2] / 1000.0
    radius = np.nanmean(rngs)

    radius = max(radius, rad_min)
    radius = min(radius, rad_max)

    rngs = [0.0] * len(det_list)
    for n in range(det_cnt): rngs[n] = sph_proj.inv(center[1], center[0], det_list[n].longitude, det_list[n].latitude, radians=False)[2] / 1000.0 - 0.01
    radius = min(radius, min(rngs))

    if radius < rad_min:
        # If the region is too small due to a nearby array, shift the center and resize to rad_min
        index = rngs.index(min(rngs))
        array_az_to_center = sph_proj.inv(det_list[index].longitude, det_list[index].latitude, center[1], center[0])[0]
        new_center = sph_proj.fwd(det_list[index].longitude, det_list[index].latitude, array_az_to_center, rad_min * 1000.0 + 0.01)
        center = [new_center[1], new_center[0]]
        radius = rad_min

    return center, radius


def calc_conf_ellipse(means, variances, conf_lvl, pnts=100):
    """Compute the confidence ellipse around a latitude longitude point

        Computes the points on an ellipse centered at a latitude, longitude point with
        variances (E/W, N/S, and covariance) defiend in kilometers for a given confidence
        level.

        Parameters
        ----------
        means : float
            Latitude and longitude of the center
        variances : float
            East/West, North/South, and covariance for the ellipse
        conf_lvl : float
            Confidence level in percentage (e.g., 95.0 = 95%)
        pnts : int
            Number of points to return around the ellipse contour

        Returns:
        ----------
        latlon : float
            Latitutde and longitude points of the ellipse

        """

    width = chi2(2).ppf(conf_lvl / 100.0)

    lat_vals = np.linspace(means[0] - variances[0] * np.sqrt(width), means[0] + variances[0] * np.sqrt(width), pnts)
    ellps = np.array([lat_vals[0], means[1] + variances[1] * (variances[2] * (lat_vals[0] - means[0]) / variances[0] + np.sqrt(max(width - ((lat_vals[0] - means[0]) / variances[0])**2, 0.0) * (1.0 - variances[2]**2)))])

    for n in range(1, pnts):          ellps = np.vstack((ellps, [lat_vals[n], means[1] + variances[1] * (variances[2] * (lat_vals[n] - means[0]) / variances[0] + np.sqrt(max(width - ((lat_vals[n] - means[0]) / variances[0])**2, 0.0) * (1.0 - variances[2]**2)))]))
    for n in range(pnts - 1, -1, -1): ellps = np.vstack((ellps, [lat_vals[n], means[1] + variances[1] * (variances[2] * (lat_vals[n] - means[0]) / variances[0] - np.sqrt(max(width - ((lat_vals[n] - means[0]) / variances[0])**2, 0.0) * (1.0 - variances[2]**2)))]))
    return np.array(ellps[:, 0]), np.array(ellps[:, 1])


def find_confidence(func, lims, conf_lvl):
    """Computes the bounds for a function given a confidence level

        Identifies the points for which \int_{x_1}^{x_2}{f(x) dx} includes a given
        fraction of the overall integral such that f(x_1) = f(x_2)

        Parameters
        ----------
        func : function of single float parameter
            Function to be analyzed
        lims : float
            Limits for the integration of the norm and bounds for possible limits

        Returns:
        ----------
        bnds : float
            The values of x_1 and x_2 for the confidence bound
        conf : float
            Actual confidence value obtained
        thresh : float
            Value of the function at x_1 and x_2

        """

    if conf_lvl > 1.0:
        print('WARNING - find_confidence cannot use conf > 1.0')
        return [lim0, lim1]

    def conf_func(x, thresh):
        val = func(x)
        if val >= thresh: return val
        else: return 0.0

    resol = 300
    x_vals = np.linspace(lims[0], lims[1], resol)
    f_vals = func(x_vals)

    f_max = max(f_vals)
    thresh_vals = np.linspace(0.0, f_max * 0.5, resol)

    norm = quad(func, lims[0], lims[1], limit=100, epsrel=1.0e-3)[0]

    conf_prev=1.0
    bnds=[]
    for n in range(resol):
        conf = quad(conf_func, lims[0], lims[1], (thresh_vals[n],), limit=100, epsrel=1.0e-3)[0] / norm

        if conf < conf_lvl < conf_prev:
            thresh = thresh_vals[n - 1] - (thresh_vals[n-1] - thresh_vals[n]) / (conf_prev - conf) * (conf_lvl - conf)
            conf = quad(conf_func, lims[0], lims[1], (thresh,), limit=100, epsrel=1.0e-3)[0] / norm

            for n in range(resol - 1):
                if (f_vals[n] - thresh) * (f_vals[n+1] - thresh) < 0.0:
                    bnds.append(x_vals[n])
            break

        conf_prev=conf

    return bnds, conf, thresh


def run(det_list, path_geo_model=None, custom_region=None, resol=180, bm_width=10.0, rng_max=np.pi / 2.0 * 6370.0, rad_min=100.0, angle=[-180,180],rad_max=1000.0, MaP_mthd="grid"):
    """Run analysis of the posterior pdf for BISL

        Compute the marginal disribution...

        Parameters
        ----------
        det_list : iterable of InfrasoundDetection instances
            Detections attributed to the event
        path_geo_model : Propagation-based, stochastic path geometry model
            Optional path geometry model if available
        custom_region : float
            Latitude, Longitude, and radius of a specific region to conduct analysis (overrides use of set_region function)
        resol : int
            Number of radial and azimuthal points used in the polar projection of the spatial PDFs
        bm_width : float
            Width of the projected beam [degrees]
        rng_max : float
            Maximmum range for beam projection [km]
        rad_min : float
            Minimum radius of the integration region [km]
        rad_max : float
            Maximum radius of the integration region [km]
        angle : float
            Minimum and maximum angles for polar projection of the spatial PDFs
        MaP_method : string
            Method to use for searching for the maximum a posteriori solution ("grid" or "random")

        Returns:
        ----------
        result : dictionary
            Dictionary containing all localization results:
                'lat_mean': Mean latitude for the marginalized spatial distribution
                'lon_mean' : Mean longitude for the marginalized spatial distribution
                'EW_var': Variance of the marginalized spatial distribution in the east/west direction in km
                'NS_var': Variance of the marginalized spatial distribution in the north/south direction in km
                'covar': Relative covariance, \sigma_{xy}^2 / (\sigma_x \sigma_y)
                't_mean': Mean marginalized temporal distribution
                't_var': Variance of the marginalized temporal distribution
                't_min' : 95% confidence bound lower limit for the marginalized temporal distribution
                't_max' : 95% confidence bound upper limit for the marginalized temporal distribution
                'lat_MaP': Maximum a Posteriori latitude
                'lon_MaP': Maximum a Posteriori longitude
                't_MaP': Maximum a Posteriori origin time
                'MaP_val' : Maximum a Posteriori value
        """

    # Determine region of interest and define the polar <--> latlon grid definition
    if custom_region:
        center = custom_region[0]
        radius = custom_region[1]
    else:
        center, radius = set_region(det_list, bm_width=bm_width, rng_max=rng_max, rad_min=rad_min, rad_max=rad_max)
    resol = int(resol)

    rngs = np.linspace(0.0, radius, resol)
    angles = np.linspace(angle[0], angle[1]-1, resol)
    #angles = np.linspace(-180.0, 179.0, resol)
    R, ANG = np.meshgrid(rngs, angles)
    proj_rngs = R.flatten()
    prof_azs = ANG.flatten()

    proj_lons, proj_lats = sph_proj.fwd(np.array([center[1]] * resol**2), np.array([center[0]] * resol**2), prof_azs, proj_rngs * 1e3)[:2]

    # Project the marginal spacial posterior in the region of interest for analysis
    pdf = lklhds.marginal_spatial_pdf(proj_lats, proj_lons, det_list, path_geo_model=path_geo_model)
    spatial_pdf = np.vstack((proj_lons, proj_lats, pdf))
    
    pdf_fit = interp2d(rngs, angles, pdf.reshape(resol, resol), kind='cubic')

    # Compute spatial means and variances
    #   dxdy -> r dr daz and absorb degrees/radians into normalization
    #   x = r sin(az), y = r cos(az)
    def integrand(r, az):
        return pdf_fit(r, az)[0] * r
    norm = nquad(integrand, [[0.0, radius],[-180.0, 179.0]], opts=[int_opts, int_opts])[0]

    def integrand(r, az):
        return r * np.sin(np.radians(az)) * pdf_fit(r, az)[0] * r
    x_mean = nquad(integrand, [[0.0, radius],[-180.0, 179.0]], opts=[int_opts,int_opts])[0] / norm

    def integrand(r, az):
        return r * np.cos(np.radians(az)) * pdf_fit(r, az)[0] * r
    y_mean = nquad(integrand, [[0.0, radius],[-180.0, 179.0]], opts=[int_opts,int_opts])[0] / norm

    temp = sph_proj.fwd(center[1], center[0], np.degrees(np.arctan2(x_mean, y_mean)), np.sqrt(x_mean**2 + y_mean**2) * 1e3)
    lat_mean, lon_mean = temp[1], temp[0]

    def integrand(r, az):
        return (r * np.sin(np.radians(az)) - x_mean)**2 * pdf_fit(r, az)[0] * r
    x_var = np.sqrt(nquad(integrand, [[0.0, radius],[-180.0, 179.0]], opts=[int_opts,int_opts])[0] / norm)

    def integrand(r, az):
        return (r * np.cos(np.radians(az)) - y_mean)**2 * pdf_fit(r, az)[0] * r
    y_var = np.sqrt(nquad(integrand, [[0.0, radius],[-180.0, 179.0]], opts=[int_opts,int_opts])[0] / norm)

    def integrand(r, az):
        return (r * np.sin(np.radians(az)) - x_mean) * (r * np.cos(np.radians(az)) - y_mean) * pdf_fit(r, az)[0] * r
    covar = nquad(integrand, [[0.0, radius],[-180.0, 179.0]], opts=[int_opts,int_opts])[0] / (norm * x_var * y_var)

    # Temporal analysis
    # Use region edge to determine limits of possible source times with celerities between 0.2 and 0.4 km/s
    time_lims = [det_list[0].peakF_UTCtime - np.timedelta64(int(rng_max / 0.2 * 1e3), 'ms'), det_list[0].peakF_UTCtime - np.timedelta64(int(0.01 / 0.4 * 1e3), 'ms')]

    conf_x, conf_y = calc_conf_ellipse([x_mean, y_mean], [x_var, y_var, covar], 99.0)
    edge_lons, edge_lats = sph_proj.fwd(np.array([center[1]] * len(conf_x)), np.array([center[0]] * len(conf_x)), np.degrees(np.arctan2(conf_x, conf_y)), np.sqrt(conf_x**2 + conf_y**2) * 1e3)[:2]

    for det in det_list:
        det_rngs = sph_proj.inv(np.array([det.longitude] * len(edge_lats)), np.array([det.latitude] * len(edge_lats)), edge_lons, edge_lats, radians=False)[2] / 1000.0
        time_lims[0] = max(det.peakF_UTCtime - np.timedelta64(int(max(det_rngs) / 0.2 * 1e3), 'ms'), time_lims[0])
        time_lims[1] = min(det.peakF_UTCtime - np.timedelta64(int(min(det_rngs) / 0.4 * 1e3), 'ms'), time_lims[1])
    #prog_bar.increment(3)

    # For each time step, project the joint likelihood at constant time and integrate out r, az
    dts = np.linspace(0.0, (time_lims[1] - time_lims[0]).astype('m8[ms]').astype(float) / 1e3, resol)
    t_vals = np.array([time_lims[0]] * len(dts))
    for n in range(len(dts)):
        t_vals[n] = t_vals[n] + np.timedelta64(int(dts[n] * 1e3), 'ms')

    time_marg_pdf = np.empty_like(dts)
    for n in range(resol):
        skip = int(np.sqrt(resol * 4.0))
        pdf = lklhds.joint_pdf(proj_lats[::skip], proj_lons[::skip], np.array([t_vals[n]] * len(proj_lats[::skip])), det_list, path_geo_model=path_geo_model)
        time_marg_pdf[n] = np.mean(pdf) * (np.pi * radius**2)

    time_pdf_fit = interp1d(dts, time_marg_pdf, kind='cubic')

    def integrand(t):
        return time_pdf_fit(t)
    time_norm = quad(integrand, dts[0], dts[-1], epsrel=1.0e-3)[0]
    time_marg_pdf /= time_norm

    temporal_pdf = [t_vals, time_marg_pdf]

    time_pdf_fit = interp1d(dts, time_marg_pdf, kind='cubic')

    # Analyze the resulting pdf to identify mean, variance, and exact 95% and 99% confidence bounds
    time_bnds_90 = find_confidence(time_pdf_fit, [dts[0], dts[-1]], 0.90)
    time_bnds_95 = find_confidence(time_pdf_fit, [dts[0], dts[-1]], 0.95)
    time_bnds_99 = find_confidence(time_pdf_fit, [dts[0], dts[-1]], 0.99)

    def integrand(t):
        return t * time_pdf_fit(t)
    time_mean = quad(integrand, dts[0], dts[-1], epsrel=1.0e-3)[0]

    def integrand(t):
        return (t - time_mean)**2 * time_pdf_fit(t)
    time_var =  np.sqrt(quad(integrand, dts[0], dts[-1], epsrel=1.0e-3)[0])

    # Maximum a Posteriori analysis
    if MaP_mthd=='random':
        n_coarse = int(1e5)

        map_lons, map_lats = sph_proj.fwd(np.array([lon_mean] * n_coarse), np.array([lat_mean] * n_coarse), np.random.uniform(-180.0, 180.0, n_coarse), abs(np.random.normal(scale=np.sqrt(x_var**2 + y_var**2), size=n_coarse)))[:2]
        map_dts = np.random.normal(loc=time_mean, scale=time_var, size=n_coarse)
    else:
        grid_resol = 200

        R, ANG = np.meshgrid(np.linspace(0.0, 3.0 * np.sqrt(x_var**2 + y_var**2), int(grid_resol)), np.linspace(-180.0, 179.0, int(grid_resol)))
        map_lons, map_lats = sph_proj.fwd(np.array([lon_mean] * int(grid_resol)**2), np.array([lat_mean] * int(grid_resol)**2), ANG.flatten(), R.flatten() * 1e3)[:2]
        map_dts = np.linspace(max(time_mean - 3.0 * time_var, 0.0), min(time_mean + 3.0 * time_var, (time_lims[1] - time_lims[0]).astype('m8[ms]').astype(float) / 1e3), grid_resol**2)

    map_tms = np.array([time_lims[0]] * len(map_dts))
    for n in range(len(map_dts)):
        map_tms[n] += np.timedelta64(int(map_dts[n] * 1e3), 'ms')

    func_vals = lklhds.joint_pdf(map_lats, map_lons, map_tms, det_list, path_geo_model=path_geo_model)

    best_est = (map_lats[np.argmax(func_vals)], map_lons[np.argmax(func_vals)], map_dts[np.argmax(func_vals)])
    def f(X): return -lklhds.joint_pdf(X[0], X[1], time_lims[0] + np.timedelta64(int(X[2] * 1e3), 'ms'), det_list, path_geo_model)
    MaP = minimize(f, best_est, method='SLSQP', options={'maxiter':1000, 'disp':False})

    result = {'lat_mean': lat_mean, 'lon_mean' : lon_mean,
              'EW_var': x_var, 'NS_var': y_var,
              'covar': covar,
              't_mean': time_lims[0] + np.timedelta64(int(time_mean * 1e3), 'ms'),
              't_var': time_var,
              't_min' : time_lims[0] + np.timedelta64(int(min(time_bnds_90[0]) * 1e3), 'ms'),
              't_max' : time_lims[0] + np.timedelta64(int(max(time_bnds_90[0]) * 1e3), 'ms'),
              'lat_MaP': MaP.x[0],
              'lon_MaP': MaP.x[1],
              't_MaP': time_lims[0] + np.timedelta64(int(MaP.x[2] * 1e3), 'ms'),
              'MaP_val' : -MaP.fun / norm,
			  'spatial_pdf' : spatial_pdf,
			  'temporal_pdf' : temporal_pdf}

    return result,spatial_pdf

def summarize(result, confidence_level=95):
    """Outputs results of BISL analysis

        Prints all results to screen in a readable format

        Parameters
        ----------
        bisl_result : dictionary
            Dictionary output of the BISL analysis methods

        """

    summary = ('Maximum a posteriori analysis: \n'
               '\tSource location: {slat}, {slon} \n'
               '\tSource time: {stime} \n'

               'Source location analysis:\n'
               '\tLatitude (mean and variance): {slatmean} +/- {nsvar} km. \n'
               '\tLongitude (mean and variance): {slonmean} +/- {ewvar} km.\n'
               '\tCovariance: ' + str(round(result['covar'],3)) + '\n'
               '\tArea of {s_confidence}% confidence ellipse: {conf} square kilometers\n'

               'Source time analysis:\n'
               '\tMean and variance: {stmean} +/- {stvar} second\n'
               '\tExact 95% confidence bounds: [{ex95confmin}, {ex95confmax}]\n'
               )
    summary = summary.format(
        slat = round(result['lat_MaP'],3),
        slon = round(result['lon_MaP'],3),
        stime = result['t_MaP'],
        slatmean = round(result['lat_mean'],3),
        nsvar = round(result['NS_var'], 3),
        slonmean = round(result['lon_mean'], 3),
        ewvar = round(result['EW_var'], 3),
        s_confidence = str(confidence_level),
        conf = round(np.pi * result['NS_var'] * result['EW_var'] * chi2(2).ppf(confidence_level/100), 3),
        stmean = result['t_mean'],
        stvar = round(result['t_var'],3),
        ex95confmin = result['t_min'],
        ex95confmax = result['t_max']
        )

    return summary
