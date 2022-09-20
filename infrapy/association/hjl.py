# infrapy.hjl.py
#
# Part of the infrapy methods which analyzes detection
# pairs to determine the normalization of their joint
# likelihood.  This value is used as the association
# parameter to identify (aggregate) detections into
# events.
#
# Author(s) Philip Blom (pblom@lanl.gov)
#           Garrett Euler (ggeuler@lanl.gov)

import time
import itertools
import random
import sys
import os

import numpy as np

import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.pyplot as plt

from datetime import datetime

from pyproj import Geod

from scipy.cluster import hierarchy
from scipy.integrate import quad, nquad
from scipy.interpolate import interp2d
from scipy.spatial.distance import squareform

from obspy import UTCDateTime

from ..propagation import likelihoods as lklhds
from ..propagation import infrasound as infrsnd

from ..utils import prog_bar
from ..utils import latlon as ll

# ################################ #
#    Set Integration Parameters    #
#     and Spherical Earth Model    #
# ################################ #
int_opts = {'limit': 50, 'epsrel': 1.0e-3}
sph_proj = Geod(ellps='sphere')

# ################################ #
#       Combining a Pair of        #
#      Detection Likelihoods       #
# ################################ #
def set_region(det1, det2, bm_width=10.0, rng_max=np.pi / 2.0 * 6370.0, rad_min=100.0, rad_max=1000.0, ):
    """Defines the integration region for computation of the joint-likelihood for a pair of detections

        Projects finite width beams from each of the detecting arrays and looks for intersections
        of the primary (center) and secondary (edge) lines to define the integration region
        for computation of the joint-likelihiood between the detection pair

        Parameters
        ----------
        det1 : InfrasoundDetection
            First detection (see infrapy.propagation.likelihoods)
        det2 : InfrasoundDetection
            Second detection (see infrapy.propagation.likelihoods)
        bm_width : float
            Width of the projected beam [degrees]
        rng_max : float
            Maximmum range for beam projection [km]
        rad_min : float
            Minimum radius of the integration region [km]
        rad_max : float
            Maximum radius of the integration region [km]

        Returns:
        Successs : boolean
            True if region was defined, False if not
        Center : float
            Center of the integration region as latitude, longitude pair [degrees]
        Radius : float
            Radius of the integration region [km]

        """

    latlon1 = np.array([det1.latitude, det1.longitude], dtype=np.float)
    latlon2 = np.array([det2.latitude, det2.longitude], dtype=np.float)

    if np.allclose(latlon1, latlon2):
        # if detections are on the same array, center is the array location and radius is the maximum value
        return True, latlon1, rad_max

    proj1 = ll.sphericalfwd(latlon1, rng_max / 6370.0 * 180.0 / np.pi, det1.back_azimuth)[0][0]
    proj1up = ll.sphericalfwd(latlon1, rng_max / 6370.0 * 180.0 / np.pi, det1.back_azimuth + bm_width)[0][0]
    proj1dn = ll.sphericalfwd(latlon1, rng_max / 6370.0 * 180.0 / np.pi, det1.back_azimuth - bm_width)[0][0]

    proj2 = ll.sphericalfwd(latlon2, rng_max / 6370.0 * 180.0 / np.pi, det2.back_azimuth)[0][0]
    proj2up = ll.sphericalfwd(latlon2, rng_max / 6370.0 * 180.0 / np.pi, det2.back_azimuth + bm_width)[0][0]
    proj2dn = ll.sphericalfwd(latlon2, rng_max / 6370.0 * 180.0 / np.pi, det2.back_azimuth - bm_width)[0][0]

    intersect = [0] * 9
    intersect[0] = ll.gcarc_intersect(latlon1, proj1, latlon2, proj2)[0]
    intersect[1] = ll.gcarc_intersect(latlon1, proj1, latlon2, proj2up)[0]
    intersect[2] = ll.gcarc_intersect(latlon1, proj1, latlon2, proj2dn)[0]
    intersect[3] = ll.gcarc_intersect(latlon1, proj1up, latlon2, proj2)[0]
    intersect[4] = ll.gcarc_intersect(latlon1, proj1dn, latlon2, proj2)[0]
    intersect[5] = ll.gcarc_intersect(latlon1, proj1up, latlon2, proj2up)[0]
    intersect[6] = ll.gcarc_intersect(latlon1, proj1dn, latlon2, proj2dn)[0]
    intersect[7] = ll.gcarc_intersect(latlon1, proj1up, latlon2, proj2dn)[0]
    intersect[8] = ll.gcarc_intersect(latlon1, proj1dn, latlon2, proj2up)[0]

    if np.count_nonzero(~np.isnan(intersect[:5])) < 1:
        return False, [0.0, 0.0], 0.0

    # Compute geographic mean of the intersections to define the center of the region
    x, y, z = 0.0, 0.0, 0.0
    for n in range(5):
        if not np.isnan(intersect[n][0]):
            x += np.cos(np.radians(intersect[n][0])) * np.cos(np.radians(intersect[n][1]))
            y += np.cos(np.radians(intersect[n][0])) * np.sin(np.radians(intersect[n][1]))
            z += np.sin(np.radians(intersect[n][0]))

    norm = np.sqrt(x**2 + y**2 + z**2)
    x /= norm
    y /= norm
    z /= norm
    center = np.array([np.degrees(np.arcsin(z)), np.degrees(np.arctan2(y, x))])

    # Check if either array is in the beam of the other array
    # Use the range to the nearer array as the initial radius of the region of interest
    temp = sph_proj.inv(det1.longitude, det1.latitude, det2.longitude, det2.latitude, radians=False)

    az_diff1 = det1.back_azimuth - temp[0]
    if az_diff1 > 180.0:
        az_diff1 -= 360.0
    if az_diff1 < -180.0:
        az_diff1 += 360.0

    az_diff2 = det2.back_azimuth - temp[1]
    if az_diff2 > 180.0:
        az_diff2 -= 360.0
    if az_diff2 < -180.0:
        az_diff2 += 360.0

    if abs(az_diff1) < bm_width and abs(az_diff2) < bm_width:
        # If both arrays are in the beams of the other, use mid-point
        # and distance for region definition

        center = np.array(np.flipud(sph_proj.fwd(det1.longitude, det1.latitude, temp[0], temp[2] / 2.0, radians=False)[:2]))
        reg_rad = (temp[2] / 1000.0) / 2.0 - 0.1
        reg_rad = min(reg_rad, rad_max)
    elif abs(az_diff1) < bm_width or abs(az_diff2) < bm_width:
        # If only one array is in the beam of the other, use
        # distance to nearer array as the radius
        rng1 = sph_proj.inv(center[1], center[0], det1.longitude, det1.latitude, radians=False)[2] / 1000.0 - 0.01
        rng2 = sph_proj.inv(center[1], center[0], det2.longitude, det2.latitude, radians=False)[2] / 1000.0 - 0.01

        reg_rad = min(rng1, rng2)
        reg_rad = min(reg_rad, rad_max)
    else:
        # Compute the distances to the non-central intersections to estimate the radius of the region of interest
        rngs = [np.nan] * 9
        for n in range(1, 9):
            if not np.isnan(intersect[n][0]):
                rngs[n] = sph_proj.inv(center[1], center[0], intersect[n][1], intersect[n][0], radians=False)[2] / 1000.0
        reg_rad = np.nanmean(rngs)

        reg_rad = max(reg_rad, rad_min)
        reg_rad = min(reg_rad, rad_max)

        reg_rad = min(reg_rad, sph_proj.inv(center[1], center[0], det1.longitude, det1.latitude)[2] / 1000.0 - 0.01)
        reg_rad = min(reg_rad, sph_proj.inv(center[1], center[0], det2.longitude, det2.latitude)[2] / 1000.0 - 0.01)

    if reg_rad < rad_min:
        # If the region is too small, shift the center and resize to rad_min
        rng1 = sph_proj.inv(center[1], center[0], det1.longitude, det1.latitude)[2] / 1000.0
        rng2 = sph_proj.inv(center[1], center[0], det2.longitude, det2.latitude)[2] / 1000.0

        if rng1 < rng2:
            array_az_to_center = sph_proj.inv(det1.longitude, det2.latitude, center[1], center[0])[0]
            new_center = sph_proj.fwd(det1.longitude, det1.latitude, array_az_to_center, rad_min * 1000.0 + 0.01)
        else:
            array_az_to_center = sph_proj.inv(det2.longitude, det2.latitude, center[1], center[0], radians=False)[0]
            new_center = sph_proj.fwd(det2.longitude, det2.latitude, array_az_to_center, rad_min * 1000.0 + 0.01)

        center = np.array([new_center[1], new_center[0]])
        reg_rad = rad_min

    return True, center, reg_rad


def compute_assoc_pair(det1, det2,  bm_width=10.0, rng_max=np.pi / 2.0 * 6370.0, rad_min=100.0, rad_max=1000.0, resol=180, prog_step=0):
    """Computes the joint-likelihiood for a pair of detections

        Projects finite width beams from each of the detecting arrays and looks for intersections
        of the primary (center) and secondary (edge) lines to define the integration region
        for computation of the joint-likelihiood between the detection pair

        Parameters
        ----------
        det1 : InfrasoundDetection
            First detection (see infrapy.propagation.likelihoods)
        det2 : InfrasoundDetection
            Second detection (see infrapy.propagation.likelihoods)
        bm_width : float
            Width of the projected beam [degrees]
        rng_max : float
            Maximmum range for beam projection [km]
        rad_min : float
            Minimum radius of the integration region [km]
        rad_max : float
            Maximum radius of the integration region [km]
        resol : int
            Number of radial and azimuthal points used in the polar projection of the likelihood PDFs
        prog_step : int
            Used to increment progress bar

        Returns:
        jntlklhd : float
            The joint-likelihood value for the pair of detections
    """

    array_sep = sph_proj.inv(det1.longitude, det1.latitude, det2.longitude, det2.latitude)[2] / 1000.0 
    if array_sep < 1.0:
        # Integration separates into a product of 1-dimensional
        # integrals for detections on the same array
        az_diff = det1.back_azimuth - det2.back_azimuth
        if az_diff > 180.0:
            az_diff -= 360.0
        elif az_diff < -180.0:
            az_diff += 360.0

        if abs(az_diff) > 2.0 * bm_width:
            jntlklhd = np.finfo(float).epsneg
        else:
            def jnt_az_pdf(az):
                arg = det1.kappa * np.cos(np.radians(az - det1.back_azimuth))
                arg += det2.kappa * np.cos(np.radians(az - det2.back_azimuth))
                return np.exp(arg) / (det1.vm_norm * det2.vm_norm)

            tms = np.array([0.0, (det2.peakF_UTCtime - det1.peakF_UTCtime).astype('m8[s]').astype(float)])
            def jnt_rng_pdf(rng):
                result = 0.0
                for indices in itertools.product(list(range(3)), repeat=2):
                    a, b, c = 0.0, 0.0, 0.0
                    coeff = 1.0

                    for n in range(2):
                        dt = tms[n] - rng * infrsnd.canon_rcel_mns[indices[n]]
                        sig = rng * infrsnd.canon_rcel_vrs[indices[n]]

                        a += 1.0 / sig**2
                        b += dt / sig**2
                        c += (dt / sig)**2
                        coeff *= infrsnd.canon_rcel_wts[indices[n]] / sig

                    result += coeff / np.sqrt(a) * np.exp(-1.0 / 2.0 * (c - b**2 / a))
                return result

            jntlklhd = 1.0 / np.sqrt(2.0 * np.pi) * quad(jnt_az_pdf, -180.0, 180.0)[0] * quad(jnt_rng_pdf, 0.0, rng_max)[0] / (np.pi * rng_max)
    else:
        # Compute integration region center and radius
        success, center, radius = set_region(det1, det2, bm_width=bm_width, rad_min=rad_min, rad_max=rad_max, rng_max=rng_max)
        if success:
            resol = int(resol)
            angles = np.linspace(-180.0, 180.0, resol)
            rngs = np.linspace(0.0, radius, resol)
            R, ANG = np.meshgrid(rngs, angles)

            R = R.flatten()
            ANG = ANG.flatten()
            temp = sph_proj.fwd(np.array([center[1]] * resol**2), np.array([center[0]] * resol**2), ANG, R * 1e3)
            pdf = lklhds.marginal_spatial_pdf(temp[1], temp[0], [det1, det2])
            pdf_fit = interp2d(rngs, angles, pdf.reshape(resol, resol), kind='cubic')

            def integrand(r, az):
                # dxdy -> r dr daz with angle changed to radians produces (r * pi / 180) factor
                return pdf_fit(r, az)[0] * r * np.pi / 180.0

            jntlklhd = nquad(integrand, [[0.0, radius],[-180.0, 180.0]], opts=[int_opts, int_opts])[0] / (np.pi * rng_max)
        else:
            jntlklhd = np.finfo(float).epsneg

    prog_bar.increment(prog_step)
    return jntlklhd


def compute_assoc_pair_wrapper(args):
    return compute_assoc_pair(*args)


def build_distance_matrix(det_list, bm_width=10.0, rng_max=np.pi / 2.0 * 6370.0, rad_min=100.0, rad_max=1000.0, resol=180,  pool=None, progress=False):
    """Computes the joint-likelihood for all pairs of detections to define the distance matrix

        Computes the joint-likelihood value for each unique pair of detections in a provided list and
        uses a negative-log-joint-likelihood to convert to non-Euclidean distance for clustering analysis

        Parameters
        ----------
        det_list : :obj:`list` of :obj:`InfrasoundDetection`
            List of detections (see infrapy.propagation.likelihoods)
        bm_width : float
            Width of the projected beam [degrees]
        rng_max : float
            Maximmum range for beam projection [km]
        rad_min : float
            Minimum radius of the integration region [km]
        rad_max : float
            Maximum radius of the integration region [km]
        resol : int
            Number of radial and azimuthal points used in the polar projection of the likelihood PDFs
        pool : pathos.multiprocessing.ProcessingPool
            Multiprocessing pool for accelerating calculations

        Returns:
        dist_matrix : 2darray
            Distance matrix describing joint-likelihood separations for all pairs

        """
    print('\tComputing joint-likelihoods...')
    det_cnt = len(det_list)
    n_tot, n_ref = det_cnt * (det_cnt - 1) / 2, 0

    if progress:
        print('\t\tProgress: \t', end='')
        prog_bar.prep(50)
    if pool:
        args, result, ids = [], [], []
        for n, m in itertools.combinations(list(range(det_cnt)), 2):
            if progress:
                step = int(np.floor((50.0 * (n_ref + 1)) / n_tot) - np.floor((50.0 * n_ref) / n_tot))
            else:
                step = 0        
            args.append([det_list[n], det_list[m], bm_width, rng_max, rad_min, rad_max, resol, step])
            ids.append([n, m])
            n_ref += 1
        out = pool.map(compute_assoc_pair_wrapper, args)
        for n in range(int(n_tot)):
            result.append([ids[n][0], ids[n][1], out[n]])
    else:
        result = []
        for n, m in itertools.combinations(list(range(det_cnt)), 2):
            if progress:
                step = int(np.floor((50.0 * (n_ref + 1)) / n_tot) - np.floor((50.0 * n_ref) / n_tot))
            else:
                step = 0
            assoc = compute_assoc_pair(det_list[n], det_list[m], bm_width=bm_width, rng_max=rng_max, rad_min=rad_min, rad_max=rad_max, resol=resol, prog_step=step)
            result.append([n, m, assoc])
            n_ref += 1
    if progress:
        prog_bar.close()

    dist_mat = np.zeros((det_cnt, det_cnt))
    for n in range(int(n_tot)):
        dist_mat[result[n][0]][result[n][1]] = -np.log10(max(result[n][2], np.finfo(float).epsneg))
        dist_mat[result[n][1]][result[n][0]] = -np.log10(max(result[n][2], np.finfo(float).epsneg))

    return dist_mat

def view_distance_matrix(distance_matrix, file_id=None, ordering=None):

    """View distance matrix used for clustering analysis

        Parameters
        ----------

        dist_matrix : 2darray
            Distance matrix describing joint-likelihood separations for all pairs

        """
    det_cnt = len(distance_matrix)

    plt.figure()
    plt.axes().set_aspect('equal')

    plt.xlim((0.0, det_cnt))
    plt.ylim((0.0, det_cnt))
    plt.xticks(np.arange(det_cnt) + 0.5, [])

    plt.ylabel("Detection ID")

    if ordering:
        plt.yticks(np.arange(det_cnt) + 0.5, ordering)
    else:
        plt.yticks(np.arange(det_cnt) + 0.5, np.arange(det_cnt))

    v1, v2 = np.meshgrid(np.arange(det_cnt), np.arange(det_cnt))
    sc_plot = plt.scatter(v1 + 0.5, v2 + 0.5, c=distance_matrix, s=3250.0 / det_cnt**1.2, cmap=cm.gnuplot2, marker="s", edgecolor='none')
    plt.colorbar(sc_plot)

    if file_id:
        plt.savefig(file_id + "-dist_mat.png", dpi=300)

    plt.show(block=False)
    plt.pause(2.0)
    plt.close()

def cluster(distance_matrix, threshold, linkage_method='weighted', show_result=False, file_id=None, den_label_size=9, mat_label_size=7, trim_indices=[]):
    """Computes the clustering solution for a distance matrix and threshold

        Computes the linkages for a distance matrix using SciPy's hierarchical (agglomerative) clustering
        methods and returns the event labels for the original detection list

        Parameters
        ----------
        distance_matrix : 2darray
            Two dimensional numpy array representing distance matrix
        threshold : float
            Threshold value defining linkage cutoff
        linkage_method : str
            Linkage method used by scipy.cluster.hierarchy.linkage
        show_results : boolean
            Boolean to plot dendrogram and sorted distance matrix to screeen
        file_id : str
            Prefix for output file is dendrogram and sorted distance matrix figure is saved
        den_label_size : float
            Font size off the dendrogram labels in the figure
        mat_label_size : float
            Font size of the distance matrix label in the figure
        trim_indices : :obj:`list` of :obj:`int` pairs
            Locations of linkages cut by trimming algorithm (only used in figure)

        Returns:
        links : scipy.cluster.hierarchy.linkage
            Output links from SciPy clustering analysis
        labels : :obj:`list` of :obj:`int`
            Labels of cluster memberships for each detection
        distance_matrix_sorted : 2darray
            The distance matrix sorted to put clusters together
        """
    det_cnt = len(distance_matrix)

    # Compute linkage and clusters and define labels
    # Note: subtract 1 so cluster indexing starts at 0
    links = hierarchy.linkage(squareform(distance_matrix), linkage_method)
    labels = hierarchy.fcluster(links, threshold, criterion='distance') - 1

    # Sort the distance matrix using the labels
    sorting = np.array([])
    for n in range(max(labels + 1)):
        sorting = np.concatenate((sorting, np.arange(det_cnt)[labels==n]))
    sorting = sorting.astype(int)

    distance_matrix_sorted = np.empty_like(distance_matrix)
    for n1 in range(det_cnt):
        for n2 in range(det_cnt):
            distance_matrix_sorted[n1][n2] = distance_matrix[sorting[n1], sorting[n2]]

    if show_result or file_id:
        f, (ax1, ax2) = plt.subplots(1, 2)
        plt.suptitle("Association Results")

        # hierarchy.set_link_color_palette(['0.5', 'b', '0.5','r', '0.5', 'g'])

        den = hierarchy.dendrogram(links, color_threshold=threshold, orientation="right", ax=ax1, leaf_font_size=den_label_size, above_threshold_color='0.5')
        ax1.axvline(x=threshold, linestyle='dashed', lw=2, color='k')

        ax2.set_aspect('equal', anchor='NE')
        ax2.set_xlim((0.0, det_cnt))
        ax2.set_ylim((0.0, det_cnt))
        plt.xticks(np.arange(det_cnt) + 0.5, [])
        plt.yticks(np.arange(det_cnt) + 0.5, sorting)
        ax2.xaxis.set_tick_params(length=1)
        ax2.yaxis.set_tick_params(length=1, labelsize=mat_label_size)

        ax1.set_ylabel("Detection ID")
        ax1.set_xlabel("Linkage Distance")

        v1, v2 = np.meshgrid(np.arange(det_cnt), np.arange(det_cnt))
        sc = ax2.scatter(v1 + 0.5, v2 + 0.5, c=distance_matrix_sorted, s=750.0 / det_cnt**1.2, cmap=cm.gnuplot2, marker="s", edgecolor='none')

        for indices in trim_indices:
            ax2.plot(np.where(sorting == indices[0])[0] + 0.5, np.where(sorting == indices[1])[0] + 0.5, color='0.5', marker='x', markersize=200.0 / det_cnt**1.2)
            ax2.plot(np.where(sorting == indices[1])[0] + 0.5, np.where(sorting == indices[0])[0] + 0.5, color='0.5', marker='x', markersize=200.0 / det_cnt**1.2)

        corner = 0
        for n in range(max(labels + 1)):
            ax2.add_patch(patches.Rectangle((corner, corner), len(np.arange(det_cnt)[labels==n]), len(np.arange(det_cnt)[labels==n]), alpha=0.2, color='k'))
            corner += len(np.arange(det_cnt)[labels==n])

        if file_id:
            plt.savefig(file_id + "-cluster_result.png", dpi=300)

        if show_result:
            plt.show(block=False)
            plt.pause(2.0)
            plt.close()

    return links, labels, distance_matrix_sorted

def trim_clusters(labels, distance_matrix, population_min=3, ratio_thresh=3.0):
    """Trims linkages in poorly shaped clusters

        Identifies poorly shaped clusters by the ratio of their mean inter-element
        distances and radius (maximum inter-element distance) and returns indices
        of the linkages that should be cut to improve clustering

        Parameters
        ----------
        labels : :obj:`list` of :obj:`int`
            Labels of cluster memberships for each detection
        distance_matrix : 2darray
            Two dimensional numpy array representing distance matrix
        population_min : int
            Minimum number of detections in a cluster to declare an event
        ratio_thresh : float
            Threshold for radius / mean inter-element distance to require trimming

        Returns:
        trim_indices : :obj:`list` of :obj:`int`
            Indices of linkages causing poor cluster shape
        """
    det_cnt = len(distance_matrix)

    trim_indices = []
    for n in range(max(labels + 1)):
        if len(np.arange(det_cnt)[labels==n]) >= population_min:
            indices = np.arange(det_cnt)[labels==n]

            distance_submatrix = np.empty((len(indices), len(indices)))
            distance_submatrix[:] = np.nan
            for n1 in range(len(indices)):
                for n2 in range(n1 + 1, len(indices)):
                    distance_submatrix[n1, n2] = distance_matrix[indices[n1]][indices[n2]]
                    distance_submatrix[n2, n1] = distance_matrix[indices[n1]][indices[n2]]
            max2median_ratio = np.nanmax(distance_submatrix) / np.nanmedian(distance_submatrix)

            if max2median_ratio > ratio_thresh:
                m1 = np.argmax([np.nanmean(col) for col in distance_submatrix])
                m2 = np.nanargmin(distance_submatrix[:, m1])
                trim_indices = trim_indices + [[indices[m1], indices[m2]]]
    
    return trim_indices


def summarize_clusters(labels, distance_matrix, population_min=3, show_result=False):
    """Prints summary of cluster association solution to screen

        Searches through labels to identify clusters with sufficient
        membership to declare events and summarizes cluster quality

        Parameters
        ----------
        labels : :obj:`list` of :obj:`int`
            Labels of cluster memberships for each detection
        distance_matrix : 2darray
            Two dimensional numpy array representing distance matrix
        population_min : int
            Minimum number of detections in a cluster to declare an event
        """
    clusters = []
    qualities = []

    det_cnt = len(distance_matrix)
    for n in range(max(labels + 1)):
        if len(np.arange(det_cnt)[labels==n]) >= population_min:
            avg_spacing, diam, cnt = 0.0, 0.0, 0
            spacing = np.zeros(len(np.arange(det_cnt)[labels==n]))

            for j, n1 in enumerate(np.arange(det_cnt)[labels==n]):
                for n2 in np.arange(det_cnt)[labels==n]:
                    avg_spacing += distance_matrix[n1][n2] / (len(np.arange(det_cnt)[labels==n]) * (len(np.arange(det_cnt)[labels==n]) - 1.0))
                    spacing[j] += distance_matrix[n1][n2] / (len(np.arange(det_cnt)[labels==n]) - 1.0)
                    diam = max(diam, distance_matrix[n1][n2])

            clusters += [list(np.arange(det_cnt)[labels==n])]
            qualities += [avg_spacing]

            if show_result:
                print('\tCluster Summary:')
                print('\tDetection IDs:', np.arange(det_cnt)[labels==n])
                print('\t\tDetection Average Spacing:', spacing)
                print('\t\tCluster Average Spacing:',  avg_spacing)
                print('\t\tCluster Diameter:',  diam, '\n')

    return clusters, qualities

def run(det_list, threshold, dist_max=10.0, bm_width=10.0, rng_max=np.pi / 2.0 * 6370.0, rad_min=100.0, rad_max=1000.0, resol=180, show_result=None, file_id=None, linkage_method='weighted', trimming_thresh=None, trim_thresh_scalar=1.0, prg_bar=False, pool=None):
    """Run the Hierarchical Joint-Likelihood (HJL) association analysis

        Runs the clustering analysis on a list of detecctions and returns the
        membership labels and sorted distance matrix for event identification

        Parameters
        ----------
        det_list : :obj:`list` of :obj:`InfrasoundDetection`
            List of detections (see infrapy.propagation.likelihoods)
        threshold : float
            Threshold value defining linkage cutoff
        dist_max: float
            Maximum value allowed in distance matrix
        bm_width : float
            Width of the projected beam [degrees]
        rng_max : float
            Maximmum range for beam projection [km]
        rad_min : float
            Minimum radius of the integration region [km]
        rad_max : float
            Maximum radius of the integration region [km]
        resol : int
            Number of radial and azimuthal points used in the polar projection of the likelihood PDFs
        show_results : boolean
            Boolean to plot dendrogram and sorted distance matrix to screeen
        file_id : str
            Prefix for output file is dendrogram and sorted distance matrix figure is saved
        linkage_method : str
            Linkage method used by scipy.cluster.hierarchy.linkage
        trim_thresh_scalar : float
            Scalar modifying the threshold value for linkage cutoff in the trimmed result
        pool : pathos.multiprocessing.ProcessingPool
            Multiprocessing pool for accelerating calculations

        Returns:
        labels : :obj:`list` of :obj:`int`
            Labels of cluster memberships for each detection
        distance_matrix_sorted : 2darray
            The distance matrix sorted to put clusters together
        """

    # Compute distance matrix via joint-likelihood values and adjust maximum distances
    if file_id:
        if os.path.isfile(file_id + "-dm.npy"):
            dists = np.load(file_id + "-dm.npy")
        else:
            dists = np.array(build_distance_matrix(det_list, bm_width=bm_width, rng_max=rng_max, rad_min=rad_min, rad_max=rad_max, resol=resol, pool=pool, progress=prg_bar))
            np.save(file_id + "-dm", dists)
        dists[dists > dist_max] = dist_max
        print('\t' + "Clustering detections into events...")
        _, labels, sorted_dists = cluster(dists, threshold, linkage_method=linkage_method, show_result=show_result, file_id=file_id + "-orig")
    else :
        dists = np.array(build_distance_matrix(det_list, bm_width=bm_width, rng_max=rng_max, rad_min=rad_min, rad_max=rad_max, resol=resol, pool=pool, progress=prg_bar))
        dists[dists > dist_max] = dist_max
        print('\t' + "Clustering detections into events...")
        _, labels, sorted_dists = cluster(dists, threshold, linkage_method=linkage_method, show_result=show_result)

    # Trim clusters with poor shape
    if trimming_thresh:
        print('\t' + "Trimming poor linkages and repeating clustering analysis...")
        dists_new = dists.copy()
        trim_indices = []
        while True:
            new_indicies = trim_clusters(labels, dists_new, ratio_thresh=trimming_thresh)
            if len(new_indicies) == 0:
                break
            else :
                trim_indices = trim_indices + new_indicies

            # Re-run clustering on trimmed matrix
            for indices in trim_indices:
                dists_new[indices[0], indices[1]] = dist_max
                dists_new[indices[1], indices[0]] = dist_max

            if file_id:
                _, labels, _ = cluster(dists_new, threshold * trim_thresh_scalar, linkage_method=linkage_method, show_result=show_result, file_id=file_id + "-trim", trim_indices=trim_indices)
            else:
                _, labels, _ = cluster(dists_new, threshold * trim_thresh_scalar, linkage_method=linkage_method, show_result=show_result, trim_indices=trim_indices)

    return labels, sorted_dists


def id_events(det_list, threshold, starttime=None, endtime=None, dist_max=10.0, bm_width=10.0, rng_max=2500.0, rad_min=100.0, rad_max=1000.0, 
    resol=180, linkage_method='weighted', trimming_thresh=3.8, cluster_det_population=3, cluster_array_population=2, prg_bar=True, pool=None):
    """
    
    
    """


    # Compute maximum propagation time and analysis window length [minutes]
    max_prop_time = int(rng_max / 0.22)
    analysis_window = int(max_prop_time * 0.5)

    # check if start and end times are defined and otherwise select from earliest and latest detections
    if starttime is None or endtime is None:
        starttime = np.min(np.array([det.peakF_UTCtime for det in det_list]))
        endtime = np.max(np.array([det.peakF_UTCtime for det in det_list]))

        starttime = UTCDateTime(starttime.astype(datetime)) - max_prop_time
        endtime = UTCDateTime(endtime.astype(datetime))
    else:
        starttime = UTCDateTime(starttime)
        endtime = UTCDateTime(endtime)
    duration = int(endtime - starttime)

    # run clustering analysis
    events, event_qls = [], []
    for dt in range(0, duration, analysis_window):
        window_start = starttime +  dt # np.timedelta64(dt, 'm')
        window_end = starttime + (dt + analysis_window + max_prop_time) # np.timedelta64(dt + int(analysis_window + max_prop_time), 'm')
        print('\n' + "Running event identification for:", window_start, "-", window_end)

        temp = [(n, det) for n, det in enumerate(det_list) if np.logical_and(window_start <= UTCDateTime(det.peakF_UTCtime.astype(datetime)), UTCDateTime(det.peakF_UTCtime.astype(datetime)) <= window_end)]
        key = [pair[0] for pair in temp]
        new_list = [pair[1] for pair in temp]

        if len(temp) >= cluster_det_population:
            labels, dists = run(new_list, threshold, dist_max=dist_max, bm_width=bm_width, rng_max=rng_max, rad_min=rad_min, rad_max=rad_max, resol=resol, 
                linkage_method=linkage_method, trimming_thresh=trimming_thresh, prg_bar=prg_bar, pool=pool)
            clusters, qualities = summarize_clusters(labels, dists, population_min=cluster_det_population)

            for n in range(len(clusters)):
                events += [[key[n] for n in clusters[n]]]
                event_qls += [10.0**(-qualities[n])]

    # clean up clusters
    print('\n' + "Cleaning up and merging clusters...")
    event_cnt = len(events)
    for n1 in range(event_cnt):
        for n2 in range(n1 + 1, event_cnt):
            if len(events[n1]) > 0 and len(events[n2]) > 0:
                set1, set2 = set(events[n1]), set(events[n2])
                rel_overlap = len(set1.intersection(set2)) / min(len(set1), len(set2))

                if rel_overlap > 0.5:
                    events[n1], events[n2] = list(set1.union(set2)), []
                    event_qls[n1], event_qls[n2] = max(event_qls[n1], event_qls[n2]), -1.0

    for n, ev_ids in enumerate(events):
        if len(ev_ids) > 0:
            locs = np.array([[det_list[j].latitude, det_list[j].longitude] for j in ev_ids])
            unique_cnt = max(len(np.unique(locs[:, 0])), len(np.unique(locs[:, 1])))
            if unique_cnt < cluster_array_population:
                events[n] = []
                event_qls[n] = -1.0

    events = [ei for ei in events if len(ei) > 0]
    event_qls = [eqi for eqi in event_qls if eqi > 0]

    return events, event_qls
