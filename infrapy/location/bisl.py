# infrapy.location.bisl.py
#
# The bisl methods utilized to analyze the posterior pdf
# obtained by combining detection likelihood pdfs.  The
# methods here use global coordiantes and assume two
# spatial parameters, latitude and longitude, along with
# the source time, t.
#
# Author            Philip Blom (pblom@lanl.gov)


import warnings

import numpy as np

from scipy.integrate import simps
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.stats import chi2

from obspy import UTCDateTime

from pyproj import Geod

from ..propagation import likelihoods as lklhds
from ..propagation import infrasound
from ..utils import latlon as ll
from ..utils import prog_bar

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

        latlon1 = np.array([det_list[n1].latitude, det_list[n1].longitude], dtype=np.float64)

        proj1 = ll.sphericalfwd(latlon1, 90, det_list[n1].back_azimuth)[0][0]
        proj1up = ll.sphericalfwd(latlon1, 90, det_list[n1].back_azimuth + bm_width)[0][0]
        proj1dn = ll.sphericalfwd(latlon1, 90, det_list[n1].back_azimuth - bm_width)[0][0]

        for n2 in range(det_cnt):
            pair_index = n1 * det_cnt + n2

            latlon2 = np.array([det_list[n2].latitude, det_list[n2].longitude], dtype=np.float64)

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
    if np.all(np.isnan(all_intersects)):
        raise ValueError('Detection set contains no beam intersects. Can not run BISL')
    
    # Compute geographic mean of the intersections to define the center of the region
    x, y, z = 0.0, 0.0, 0.0
    for n in range (9):
        if not np.isnan(all_intersects[n][0]):
            x += np.cos(np.radians(all_intersects[n][0])) * np.cos(np.radians(all_intersects[n][1]))
            y += np.cos(np.radians(all_intersects[n][0])) * np.sin(np.radians(all_intersects[n][1]))
            z += np.sin(np.radians(all_intersects[n][0]))

    norm = np.sqrt(x**2 + y**2 + z**2)
    x /= norm
    y /= norm
    z /= norm
    center = [np.degrees(np.arcsin(z)), np.degrees(np.arctan2(y, x))]

    # Compute the distance to the non-central intersections to estimate the size of the region
    # Limit the size of the possible region (diameter of 100 - 1,000 km, can't include a detecting station)
    rngs = [np.nan] * len(all_intersects)
    for n in range (len(main_intersects),len(all_intersects)):
        if not np.isnan(all_intersects[n][0]): rngs[n] = sph_proj.inv(center[1], center[0], all_intersects[n][1], all_intersects[n][0], return_back_azimuth=True, radians=False)[2] / 1000.0
    radius = np.nanmean(rngs)

    radius = max(radius, rad_min)
    radius = min(radius, rad_max)

    rngs = [0.0] * len(det_list)
    for n in range(det_cnt): rngs[n] = sph_proj.inv(center[1], center[0], det_list[n].longitude, det_list[n].latitude, return_back_azimuth=True, radians=False)[2] / 1000.0 - 0.01
    radius = min(radius, min(rngs))

    if radius < rad_min:
        # If the region is too small due to a nearby array, shift the center and resize to rad_min
        index = rngs.index(min(rngs))
        array_az_to_center = sph_proj.inv(det_list[index].longitude, det_list[index].latitude, center[1], center[0], return_back_azimuth=True)[0]
        new_center = sph_proj.fwd(det_list[index].longitude, det_list[index].latitude, array_az_to_center, rad_min * 1000.0 + 0.01, return_back_azimuth=True)
        center = [new_center[1], new_center[0]]
        radius = rad_min

    return center, radius


def calc_conf_ellipse(means, st_devs, conf_lvl, pnts=100):
    """Compute the confidence ellipse around a latitude longitude point

        Computes the points on an ellipse centered at a latitude, longitude point with
        standard deviations (E/W, N/S, and covariance) defiend in kilometers for a given confidence
        level.

        Parameters
        ----------
        means : float
            Latitude and longitude of the center
        st_devs : float
            East/West and North/South standard deviations, and covariance for the ellipse
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

    lat_vals = np.linspace(means[0] - st_devs[0] * np.sqrt(width), means[0] + st_devs[0] * np.sqrt(width), pnts)
    ellps = np.array([lat_vals[0], means[1] + st_devs[1] * (st_devs[2] * (lat_vals[0] - means[0]) / st_devs[0] + np.sqrt(max(width - ((lat_vals[0] - means[0]) / st_devs[0])**2, 0.0) * (1.0 - st_devs[2]**2)))])

    for n in range(1, pnts):          ellps = np.vstack((ellps, [lat_vals[n], means[1] + st_devs[1] * (st_devs[2] * (lat_vals[n] - means[0]) / st_devs[0] + np.sqrt(max(width - ((lat_vals[n] - means[0]) / st_devs[0])**2, 0.0) * (1.0 - st_devs[2]**2)))]))
    for n in range(pnts - 1, -1, -1): ellps = np.vstack((ellps, [lat_vals[n], means[1] + st_devs[1] * (st_devs[2] * (lat_vals[n] - means[0]) / st_devs[0] - np.sqrt(max(width - ((lat_vals[n] - means[0]) / st_devs[0])**2, 0.0) * (1.0 - st_devs[2]**2)))]))
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
        return [lims[0], lims[1]]

    def conf_func(x, thresh):
        val = func(x)
        if val >= thresh:
            return val
        else:
            return 0.0

    resol = 200
    x_vals = np.linspace(lims[0], lims[1], resol)
    f_vals = func(x_vals)

    f_max = max(f_vals)
    thresh_vals = np.linspace(0.0, f_max * 0.5, resol)

    norm = simps(f_vals, x_vals)

    conf_prev=1.0
    bnds=[]
    for n in range(resol):
        conf = simps(np.array([conf_func(xj, thresh_vals[n]) for xj in x_vals]), x_vals) / norm

        if conf < conf_lvl < conf_prev:
            thresh = thresh_vals[n - 1] - (thresh_vals[n-1] - thresh_vals[n]) / (conf_prev - conf) * (conf_lvl - conf)
            conf = simps(np.array([conf_func(xj, thresh_vals[n]) for xj in x_vals]), x_vals) / norm

            for n in range(resol - 1):
                if (f_vals[n] - thresh) * (f_vals[n+1] - thresh) < 0.0:
                    bnds.append(x_vals[n])
            break

        conf_prev=conf

    return bnds, conf, thresh


def run(det_list, path_geo_model=None, custom_region=None, latlon_resol=0.05, tm_resol = 60.0, bm_width=10.0, rng_max=np.pi / 2.0 * 6370.0, rad_min=100.0, angle=[-180,180],rad_max=1000.0, verbose=True):
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

        Returns:
        ----------
        result : dictionary
            Dictionary containing all localization results:
                'lat_mean': Mean latitude for the marginalized spatial distribution
                'lon_mean' : Mean longitude for the marginalized spatial distribution
                'EW_stdev': Standard deviation of the marginalized spatial distribution in the east/west direction in km
                'NS_stdev': Standard deviation of the marginalized spatial distribution in the north/south direction in km
                'covar': Relative covariance, \sigma_{xy}^2 / (\sigma_x \sigma_y)
                't_mean': Mean marginalized temporal distribution
                't_stdev': Standard deviation of the marginalized temporal distribution
                't_min' : 95% confidence bound lower limit for the marginalized temporal distribution
                't_max' : 95% confidence bound upper limit for the marginalized temporal distribution
                'lat_MaP': Maximum a Posteriori latitude
                'lon_MaP': Maximum a Posteriori longitude
                't_MaP': Maximum a Posteriori origin time
                'MaP_val' : Maximum a Posteriori value
        """

    if verbose:
        print("Running Bayesian Infrasonic Source Localization (BISL) Analysis...")

    try:
        # Determine region of interest and define the polar <--> latlon grid definition
        if custom_region:
            center = (custom_region[0], custom_region[1])
            radius = custom_region[2]
        else:
            az_cnt = sum(det.back_azimuth is not None for det in det_list)
            if az_cnt > 1:
                center, radius = set_region(det_list, bm_width=bm_width, rng_max=rng_max, rad_min=rad_min, rad_max=rad_max)
            else:
                msg = "Detection set doesn't include at least 3 direction-of-arrival detections.  Analysis requires source region definition: --src-est '(lat, lon, radius)'"
                raise ValueError(msg)
    except Exception as e:
        raise ValueError(str(e)) from e
    
    resol = 200
    if verbose:
        print('\t' + "Identifying integration region...")
    rngs = np.linspace(0.0, radius, resol)
    angles = np.linspace(angle[0], angle[1] - 1, resol)
    
    R, ANG = np.meshgrid(rngs, angles)
    proj_rngs = R.flatten()
    prof_azs = ANG.flatten()

    proj_lons, proj_lats = sph_proj.fwd(np.array([center[1]] * resol**2), np.array([center[0]] * resol**2), prof_azs, proj_rngs * 1e3, return_back_azimuth=True)[:2]

    lat_vals = np.arange(np.min(proj_lats), np.max(proj_lats), latlon_resol)
    lon_vals = np.arange(np.min(proj_lons), np.max(proj_lons), latlon_resol)
    
    if len([det for det in det_list if det.peakF_UTCtime < UTCDateTime("9999-01-01T00:00:00")]) > 0:
        if verbose:
            print('\t' + "Defining grid and evaluating likelihoods...")
        t_mins, t_maxs = [], []
        for det in det_list:
            rngs = np.array(sph_proj.inv(det.longitude * np.ones_like(proj_lons), det.latitude * np.ones_like(proj_lats), proj_lons, proj_lats, radians=False))[2] / 1000.0             
            t_mins = t_mins + [det.peakF_UTCtime - np.timedelta64(int(np.max(rngs) / 0.22 * 1.0e3), 'ms')]
            t_maxs = t_maxs + [det.peakF_UTCtime - np.timedelta64(int(np.min(rngs) / 0.34 * 1.0e3), 'ms')]
        tm_lims = [max(t_mins), min(t_maxs)]

        dt_vals = np.arange(0.0, (tm_lims[1] - tm_lims[0]).astype('m8[s]').astype(float), tm_resol)
        tm_vals = np.array([tm_lims[0] + np.timedelta64(int(dt_vals[n] * 1e3), 'ms') for n in range(len(dt_vals))])

        lat_grid, lon_grid, tm_grid = np.meshgrid(lat_vals, lon_vals, tm_vals, indexing='ij')     

        if verbose:
            print('\t\t Progress: ', end='')
            prog_bar.prep(5 * len(det_list))
            pdf = np.array([det.pdf(lat_grid, lon_grid, tm_grid, path_geo_model=path_geo_model, prog_step=5) for det in det_list if type(det) == lklhds.InfrasoundDetection]).prod(axis=0)
            prog_bar.close()
        else:
            pdf = np.array([det.pdf(lat_grid, lon_grid, tm_grid, path_geo_model=path_geo_model) for det in det_list if type(det) == lklhds.InfrasoundDetection]).prod(axis=0)

        if verbose:
            print('\t' + "Analyzing localization pdf...")
            print('\t\t' + "Normalizing and marginalizing...")
        spatial_pdf = simps(pdf, x=dt_vals)
        tm_pdf = simps(simps(pdf * np.cos(np.radians(lat_grid)), x=lon_vals, axis=1), x=lat_vals, axis=0)

        norm = simps(tm_pdf, dt_vals)
        spatial_pdf = spatial_pdf / norm
        tm_pdf = tm_pdf / norm

        if verbose:
            print('\t\t' + "Analyzing spatial PDF...")
        def simps_spatial(vals):
            result = simps(vals, x=lon_vals)
            result = simps(result * np.cos(np.radians(lat_vals)), x=lat_vals)
            return result
        
        lat_grid2, lon_grid2 = np.meshgrid(lat_vals, lon_vals, indexing='ij') 
        lat_mean, lon_mean = [simps_spatial(grid_vals * spatial_pdf) for grid_vals in [lat_grid2, lon_grid2]]

        temp = np.array(sph_proj.inv(lon_mean * np.ones_like(lon_grid2), lat_mean * np.ones_like(lat_grid2), lon_grid2, lat_grid2, return_back_azimuth=True, radians=False))
        dx, dy = temp[2] / 1000.0 * np.sin(np.radians(temp[0])), temp[2] / 1000.0 * np.cos(np.radians(temp[0]))

        x_stdev, y_stdev = [np.sqrt(simps_spatial(diff**2 * spatial_pdf)) for diff in [dx, dy]]
        covar = simps_spatial(dx * dy * spatial_pdf) / (x_stdev * y_stdev)

        if verbose:
            print('\t\t' + "Analyzing temporal PDF...")
        dt_mean = simps(dt_vals * tm_pdf, x=dt_vals)
        dt_stdev = np.sqrt(simps((dt_vals - dt_mean)**2 * tm_pdf, x=dt_vals))

        tm_mask = np.logical_and(dt_mean - 3.5 * dt_stdev < dt_vals, dt_vals < dt_mean + 3.5 * dt_stdev)
        time_bnds_90 = find_confidence(interp1d(dt_vals[tm_mask], tm_pdf[tm_mask], kind='cubic'), [dt_vals[tm_mask][0], dt_vals[tm_mask][-1]], 0.90)
        
        MaP_index = np.argmax(spatial_pdf.flatten())

        result = {'lat_mean': lat_mean, 'lon_mean' : lon_mean,
                    'EW_stdev': x_stdev, 'NS_stdev': y_stdev,
                    'covar': covar,
                    't_mean': tm_lims[0] + np.timedelta64(int(dt_mean * 1e3), 'ms'),
                    't_stdev': dt_stdev,
                    't_min' : tm_lims[0] + np.timedelta64(int(min(time_bnds_90[0]) * 1e3), 'ms'),
                    't_max' : tm_lims[0] + np.timedelta64(int(max(time_bnds_90[0]) * 1e3), 'ms'),
                    'lat_MaP': lat_grid.flatten()[MaP_index],
                    'lon_MaP': lon_grid.flatten()[MaP_index],
                    't_MaP': tm_grid.flatten()[MaP_index],
                    'MaP_val' : spatial_pdf.flatten()[MaP_index],
                    'spatial_pdf' : np.vstack((lon_grid2.flatten(), lat_grid2.flatten(), spatial_pdf.flatten())),
                    'temporal_pdf' : [tm_vals, tm_pdf]}

    else:
        if verbose:
            print('\t' + "Defining grid and evaluating likelihoods...")
        lat_grid, lon_grid = np.meshgrid(lat_vals, lon_vals, indexing='ij')

        if verbose:
            print('\t\t Progress: ', end='')
            prog_bar.prep(5 * len(det_list))
            pdf = np.array([det.az_pdf(lat_grid, lon_grid, path_geo_model=path_geo_model, prog_step=5) for det in det_list if type(det) == lklhds.InfrasoundDetection]).prod(axis=0)
            prog_bar.close()
        else:
            pdf = np.array([det.az_pdf(lat_grid, lon_grid, path_geo_model=path_geo_model) for det in det_list if type(det) == lklhds.InfrasoundDetection]).prod(axis=0)

        if verbose:
            print('\t' + "Analyzing localization pdf...")
        def simps_2dim(vals):
            result = simps(vals, x=lon_vals)
            result = simps(result, x=lat_vals)
            return result
        
        integrand = pdf * np.cos(np.radians(lat_grid))
        integrand = integrand / simps_2dim(integrand)

        # Spatial statistics
        lat_mean, lon_mean = [simps_2dim(grid_vals * integrand) for grid_vals in [lat_grid, lon_grid]]

        temp = np.array(sph_proj.inv(lon_mean * np.ones_like(lon_grid), lat_mean * np.ones_like(lat_grid), lon_grid, lat_grid, return_back_azimuth=True, radians=False))
        dx = temp[2] / 1000.0 * np.sin(np.radians(temp[0]))
        dy = temp[2] / 1000.0 * np.cos(np.radians(temp[0]))

        x_stdev = np.sqrt(simps_2dim(dx**2 * integrand))
        y_stdev = np.sqrt(simps_2dim(dy**2 * integrand))
        covar = simps_2dim(dx * dy * integrand) / (x_stdev * y_stdev)

        LAT, LON = np.meshgrid(np.unique(lat_grid), np.unique(lon_grid))
        spatial_pdf = np.vstack((LON.flatten(), LAT.flatten(), pdf.flatten()))

        MaP_index = np.argmax(spatial_pdf.flatten())

        result = {'lat_mean': lat_mean, 'lon_mean' : lon_mean,
                    'EW_stdev': x_stdev, 'NS_stdev': y_stdev,
                    'covar': covar,
                    'lat_MaP': lat_grid.flatten()[MaP_index],
                    'lon_MaP': lon_grid.flatten()[MaP_index],
                    'MaP_val' : spatial_pdf.flatten()[MaP_index],
                    'spatial_pdf' : spatial_pdf}
        
    return result 
    

def summarize(result, confidence_level=90):
    """Outputs results of BISL analysis

        Prints all results to screen in a readable format

        Parameters
        ----------
        bisl_result : dictionary
            Dictionary output of the BISL analysis methods

        """

    if 't_MaP' in result:
        if confidence_level != 90:
            dt_vals = np.array([(tm_val - result['temporal_pdf'][0][0]).astype('m8[ms]').astype(float) / 1.0e3 for tm_val in result['temporal_pdf'][0]])

            dt_mean = (result['t_mean'] - result['temporal_pdf'][0][0]).astype('m8[ms]').astype(float) / 1.0e3
            tm_mask = np.logical_and(dt_mean - 4.0 * result['t_stdev'] < dt_vals, dt_vals < dt_mean + 4.0 * result['t_stdev'])

            tm_conf = find_confidence(interp1d(dt_vals[tm_mask], np.array(result['temporal_pdf'][1])[tm_mask], kind='cubic'), [dt_vals[tm_mask][0], dt_vals[tm_mask][-1]], confidence_level / 100.0)            

            tm_min_val = result['temporal_pdf'][0][0] + np.timedelta64(int(min(tm_conf[0]) * 1e3), 'ms')
            tm_max_val = result['temporal_pdf'][0][0] + np.timedelta64(int(max(tm_conf[0]) * 1e3), 'ms') 
        else:
            tm_min_val = result['t_min']
            tm_max_val = result['t_max']

        summary = ('Maximum a posteriori analysis: \n'
                    '\tSource location: {slat}, {slon} \n'
                    '\tSource time: {stime} \n'

                    'Source location analysis:\n'
                    '\tLatitude (mean and standard deviation): {slatmean} +/- {nsvar} km. \n'
                    '\tLongitude (mean and standard deviation): {slonmean} +/- {ewvar} km.\n'
                    '\tCovariance: {covar}.\n'
                    '\tArea of {s_confidence}% confidence ellipse: {conf} square kilometers\n'

                   'Source time analysis:\n'
                   '\tMean and standard deviation: {stmean} +/- {stvar} second\n'
                   '\tExact {s_confidence}% confidence bounds: [{ex95confmin}, {ex95confmax}]\n'
                   )

        summary = summary.format(
            slat = np.round(result['lat_MaP'], 3),
            slon = np.round(result['lon_MaP'], 3),
            stime = result['t_MaP'],
            slatmean = np.round(result['lat_mean'], 3),
            nsvar = np.round(result['NS_stdev'], 3),
            slonmean = np.round(result['lon_mean'], 3),
            ewvar = np.round(result['EW_stdev'], 3),
            s_confidence = str(confidence_level),
            covar = np.round(result['covar'], 3),
            conf = np.round(np.pi * result['NS_stdev'] * result['EW_stdev'] * chi2(2).ppf(confidence_level / 100.0), 3),
            stmean = result['t_mean'],
            stvar = np.round(result['t_stdev'],3),
            ex95confmin = tm_min_val,
            ex95confmax = tm_max_val
            )
        
    else:
        summary = ('Maximum a posteriori analysis: \n'
                    '\tSource location: {slat}, {slon} \n\n'

                    'Source location analysis:\n'
                    '\tLatitude (mean and standard deviation): {slatmean} +/- {nsvar} km. \n'
                    '\tLongitude (mean and standard deviation): {slonmean} +/- {ewvar} km.\n'
                    '\tCovariance: {covar}.\n'
                    '\tArea of {s_confidence} confidence ellipse: {conf} square kilometers\n'
                   )

        summary = summary.format(
            slat = np.round(result['lat_MaP'], 3),
            slon = np.round(result['lon_MaP'], 3),
            slatmean = np.round(result['lat_mean'], 3),
            nsvar = np.round(result['NS_stdev'], 3),
            slonmean = np.round(result['lon_mean'], 3),
            ewvar = np.round(result['EW_stdev'], 3),
            s_confidence = str(confidence_level),
            covar = np.round(result['covar'], 3),
            conf = np.round(np.pi * result['NS_stdev'] * result['EW_stdev'] * chi2(2).ppf(confidence_level / 100.0), 3),
            )

    return summary
