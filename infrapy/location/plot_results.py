# infrapy.location.plot_results.py
#
# A series of functions to plot results from bisl location methods utilizing basemap
# transition from basemap to cartopy to come
#
# @fkdd  fransiska at lanl dot gov

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
from infrapy.location import bisl
from ..propagation import likelihoods as lklhds
from ..utils import latlon as ll
from IPython import embed

import matplotlib.pyplot as plt



#from mpl_toolkits.basemap import Basemap
#from matplotlib import cm
from pyproj import Geod
wgs84 = Geod(ellps='WGS84')

int_opts = {'limit': 100, 'epsrel': 1.0e-3}
sph_proj = Geod(ellps='sphere')

import pandas as pd
from infrapy.utils import latlon
from IPython import embed

def plot_dets(det_list,bm_width, rng_max, rad_min, rad_max,display=True,save_fig=None):

    center, radius = bisl.set_region(det_list, bm_width=bm_width, rng_max=rng_max, rad_min=rad_min, rad_max=rad_max)

    det_cnt = len(det_list)
    fig1 = plt.figure(1)

    ax = fig1.add_subplot(1, 1, 1)


    array_lats, array_lons = [0] * len(det_list), [0] * len(det_list)
    for n in range(len(det_list)):
        array_lats[n], array_lons[n] = det_list[n]._InfrasoundDetection__lat, det_list[n]._InfrasoundDetection__lon

    fig_bnd_lat_min, fig_bnd_lat_max = min(array_lats) - 1.0, max(array_lats) + 1.0
    fig_bnd_lon_min, fig_bnd_lon_max = min(array_lons) - 2.0, max(array_lons) + 2.0
    fig_cntr_lat = 0.5 * (fig_bnd_lat_max + fig_bnd_lat_min)
    fig_cntr_lon = 0.5 * (fig_bnd_lon_max + fig_bnd_lon_min)

    rng = wgs84.inv(fig_bnd_lon_min, fig_cntr_lat, fig_bnd_lon_max, fig_cntr_lat, radians=False)[2] / 1000.0 / 5.0
    rng_scale = max(int(round(rng/50.0) * 50.0), 50.0)
    lon_scale = fig_bnd_lon_min + 0.125 * (fig_bnd_lon_max - fig_bnd_lon_min)
    lat_scale = fig_bnd_lat_min + 0.1 * (fig_bnd_lat_max - fig_bnd_lat_min)

    m1 = Basemap(llcrnrlon=fig_bnd_lon_min, llcrnrlat=fig_bnd_lat_min, urcrnrlon=fig_bnd_lon_max, urcrnrlat=fig_bnd_lat_max, resolution='l', projection='cass', lon_0=center[1], lat_0=center[0])
    # m1 = Basemap(projection='nsper',lon_0=center[1],lat_0=center[0], satellite_height=1.0e7, resolution='l')
    m1.drawmapboundary(linewidth=1.5)
    m1.drawcoastlines(linewidth=0.75)
    m1.drawcountries(linewidth=0.75)
    m1.drawstates(linewidth=0.75)
    m1.drawmapscale(lon_scale, lat_scale, fig_cntr_lon, fig_cntr_lat, rng_scale, barstyle='fancy')

    circle = np.array(wgs84.fwd(np.array([center[1]] * 181), np.array([center[0]] * 181), np.linspace(-180.0, 179.0, 181), np.array([radius * 1e3] * 181), radians=False))
    m1.plot(circle[0], circle[1], color='blue', linestyle='--',linewidth=2.0, latlon=True)

    beams = [0] * det_cnt
    for n, det in enumerate(det_list):
        proj = latlon.sphericalfwd([det._InfrasoundDetection__lat, det._InfrasoundDetection__lon], 30, det._InfrasoundDetection__back_az)[0][0]
        beams[n], = m1.drawgreatcircle(det._InfrasoundDetection__lon, det._InfrasoundDetection__lat, proj[1], proj[0], linewidth=0.5, color='Blue')

    m1.plot(array_lons, array_lats, 'y^', markersize=10.0, latlon=True)

    plt.title('Back Azimuth Projection')

    if display == True:
        plt.show(block=False)
        plt.pause(0.1)

    if save_fig:

        plt.figure(1)
        plt.savefig(str(save_fig) + '_detections_priors.png')


    if display == True:
        plt.close('all')

def plot_loc(det_list,result,pdf,bm_width, rng_max, rad_min, rad_max,angle,confidence,display=True,save_fig=None,grnd_trth=None):

    plt.figure(2)

    center, radius = bisl.set_region(det_list, bm_width=bm_width, rng_max=rng_max, rad_min=rad_min, rad_max=rad_max)
    #resol = 316
    #rngs = np.linspace(0.0, radius, resol)
    #angles = np.linspace(angle[0], angle[1]-1, resol)
    #angles = np.linspace(-180.0, 179.0, resol)
    #R, ANG = np.meshgrid(rngs, angles)
    #proj_rngs = R.flatten()
    #prof_azs = ANG.flatten()

    #proj_lons, proj_lats = sph_proj.fwd(np.array([center[1]] * resol**2), np.array([center[0]] * resol**2), prof_azs, proj_rngs * 1e3)[:2]
    #pdf = lklhds.marginal_spatial_pdf(proj_lats, proj_lons, det_list, path_geo_model=path_geo_model)
    m2 = Basemap(width=2.0 * radius * 1e3, height=2.0 * radius * 1e3, resolution='i', projection='stere', lat_ts=center[0], lat_0=center[0], lon_0=center[1])

    m2.drawmapboundary(linewidth=2.0)
    m2.drawcoastlines(linewidth=1.5)
    m2.drawcountries(linewidth=1.5)
    m2.drawstates(linewidth=1.5)
    circle = np.array(wgs84.fwd(np.array([center[1]] * 181), np.array([center[0]] * 181), np.linspace(-180.0, 179.0, 181), np.array([radius * 1e3] * 181), radians=False))
    rng_scale = max(int(round(radius / 5 / 50.0) * 50.0), 50.0)
    lat_scale = min(circle[1]) + 0.1 * (max(circle[1]) - min(circle[1]))
    lon_scale = min(circle[0]) + 0.1 * (max(circle[0]) - min(circle[0]))
    m2.drawmapscale(lon_scale, lat_scale, center[1], center[0], rng_scale, barstyle='fancy')

    plt.title('Spatial PDF and Confidence Bounds')


    m2.scatter(pdf[0], pdf[1], c=pdf[2], s = 5.0, latlon=True, cmap=cm.Blues, marker="s", alpha=0.5, edgecolor='none', vmin=0.0)
    for aa in range(len(confidence)):
        conf = confidence[aa]
        conf_x, conf_y = bisl.calc_conf_ellipse([0.0, 0.0],[result['EW_var'], result['NS_var'], result['covar']],conf)
        '''
        conf_lons, conf_lats = sph_proj.fwd(np.array([result['lon_mean']] * len(conf_dx)),
                                            np.array([result['lat_mean']] * len(conf_dy)),
                                            np.degrees(np.arctan2(conf_dx, conf_dy)),
                                            np.sqrt(conf_dx**2 + conf_dy**2) * 1e3)[:2]
        '''

        temp = wgs84.fwd(np.array([center[1]] * len(conf_x)), np.array([center[0]] * len(conf_x)), np.degrees(np.arctan2(conf_x, conf_y)), np.sqrt(conf_x**2 + conf_y**2) * 1e3)
        m2.plot(temp[0], temp[1], color='DarkGreen', linewidth=2.5, latlon=True)
    m2.plot(result['lon_MaP'], result['lat_MaP'], 'm.', markersize=12.5, latlon=True)

    if grnd_trth: m2.plot([grnd_trth[1]], [grnd_trth[0]], 'r*', markersize=10.0, latlon=True)


    if display == True:
        plt.show(block=False)
        plt.pause(0.1)

    if save_fig:

        plt.figure(2)
        plt.savefig(str(save_fig) + '_location_priors.png')
        #plt.close('all')
    if display == True:
        plt.close('all')

def plot_time(result,grnd_trth,display=False,save_fig=None):

     plt.figure(3)
     plt.title('Source Time Distribution')
     t_pdf = result['temporal_pdf']
     plt.plot(t_pdf[0], t_pdf[1], label='Temporal PDF')
     plt.axvline(x=result['t_MaP'], color='magenta', label='Max a posteriori time')
     plt.axvline(x=result['t_min'], color='green', ls='--', label='95% Conf. Bnd Limits')
     plt.axvline(x=result['t_max'], color='green', ls='--')

     if grnd_trth is not None:
       plt.axvline(x=grnd_trth[2], color='red', ls='-', lw=1.5, label='Ground Truth')


     if display==True:
         plt.show(block=False)
         plt.pause(0.1)

     if save_fig:

         plt.figure(3)
         plt.savefig(str(save_fig) + '_source_time_priors.png')
         plt.close('all')

     if display==True:
         plt.close('all')
