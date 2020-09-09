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
from ..utils import scale_bar
from IPython import embed

import matplotlib.pyplot as plt
from matplotlib import cm

int_opts = {'limit': 100, 'epsrel': 1.0e-3}


import pandas as pd
from infrapy.utils import latlon
from IPython import embed

import cartopy.crs as ccrs
import cartopy
sph_proj = Geod(ellps="WGS84")
pc_proj = ccrs.PlateCarree()
import cartopy.feature as cfeature


def plot_dets(det_list,display=True,save_fig=None):

    det_cnt = len(det_list)
    fig1 = plt.figure(1)

    ax = fig1.add_subplot(1, 1, 1,projection=pc_proj)
    ax.coastlines()
    #ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.STATES)

    array_lats, array_lons = [0] * len(det_list), [0] * len(det_list)
    for n in range(len(det_list)):
        array_lats[n], array_lons[n] = det_list[n]._InfrasoundDetection__lat, det_list[n]._InfrasoundDetection__lon


    fig_bnd_lat_min, fig_bnd_lat_max = min(array_lats) - 1.0, max(array_lats) + 1.0
    fig_bnd_lon_min, fig_bnd_lon_max = min(array_lons) - 2.0, max(array_lons) + 2.0

    ax.set_extent([fig_bnd_lon_min, fig_bnd_lon_max, fig_bnd_lat_min, fig_bnd_lat_max],crs=pc_proj)
    scale_bar.scale_bar(ax, (0.05, 0.05), 300)
    #embed()

    beams = [0] * det_cnt
    for n, det in enumerate(det_list):
        proj = latlon.sphericalfwd([det._InfrasoundDetection__lat, det._InfrasoundDetection__lon], 50, det._InfrasoundDetection__back_az)[0][0]
        #beams[n], = m1.drawgreatcircle(det._InfrasoundDetection__lon, det._InfrasoundDetection__lat, proj[1], proj[0], linewidth=0.5, color='Blue')
        plt.plot((det._InfrasoundDetection__lon, proj[1]),(det._InfrasoundDetection__lat, proj[0]), linewidth=0.5, color='Blue',transform=pc_proj)
        #fig1.lines.append(line)
    ax.plot(array_lons, array_lats, 'y^', markersize=12, transform=pc_proj)
    plt.title('Back Azimuth Projection')

    if display == True:
        plt.show(block=False)
        plt.pause(0.5)

    if save_fig:

        #plt.figure(1)
        fig1.savefig(str(save_fig) + '_detections.png')


    if display == True:
        plt.close('all')

def plot_loc(det_list,result,pdf,confidence,display=True,save_fig=None,grnd_trth=None, include_arrays=False):


    det_cnt = len(det_list)
    fig2 = plt.figure(2)

    ax = fig2.add_subplot(1, 1, 1,projection=pc_proj)
    ax.coastlines()
    #ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAKES)
    # ax.add_feature(cfeature.STATES)
    ax.add_feature(cfeature.BORDERS)

    array_lats, array_lons = [0] * len(det_list), [0] * len(det_list)
    for n in range(len(det_list)):
        array_lats[n], array_lons[n] = det_list[n]._InfrasoundDetection__lat, det_list[n]._InfrasoundDetection__lon


    if include_arrays:
        fig_bnd_lat_min, fig_bnd_lat_max = min(array_lats) - 1.0, max(array_lats) + 1.0
        fig_bnd_lon_min, fig_bnd_lon_max = min(array_lons) - 2.0, max(array_lons) + 2.0
    else:
        fig_bnd_lat_min, fig_bnd_lat_max = result['lat_mean'] - 7.5 * result['NS_stdev'] / 111.0, result['lat_mean'] + 7.5 * result['NS_stdev'] / 111.0
        fig_bnd_lon_min, fig_bnd_lon_max = result['lon_mean'] - 7.5 * result['EW_stdev'] / 111.0, result['lon_mean'] + 7.5 * result['EW_stdev'] / 111.0


    ax.set_extent([fig_bnd_lon_min, fig_bnd_lon_max, fig_bnd_lat_min, fig_bnd_lat_max],crs=pc_proj)
    scale_bar.scale_bar(ax, (0.05, 0.05), 250)

    plt.title('Spatial PDF and Confidence Bounds')

    ax.scatter(pdf[0], pdf[1], c=pdf[2], s = 5.0, cmap=cm.Blues, marker="s", alpha=0.5, edgecolor='none', vmin=0.0,transform=pc_proj)
    for aa in range(len(confidence)):
        conf = confidence[aa]
        conf_x, conf_y = bisl.calc_conf_ellipse([0.0, 0.0],[result['EW_stdev'], result['NS_stdev'], result['covar']],conf)

        conf_lons, conf_lats = sph_proj.fwd(np.array([result['lon_mean']] * len(conf_x)),
                                            np.array([result['lat_mean']] * len(conf_y)),
                                            np.degrees(np.arctan2(conf_x, conf_y)),
                                            np.sqrt(conf_x**2 + conf_y**2) * 1e3)[:2]

        ax.plot(conf_lons, conf_lats, color='green', transform=pc_proj)

    ax.plot(result['lon_MaP'], result['lat_MaP'], 'm.', markersize=12.5, transform=pc_proj)

    if grnd_trth: ax.plot([grnd_trth[1]], [grnd_trth[0]], 'r*', markersize=10.0, transform=pc_proj)


    if display == True:
        plt.show(block=False)
        plt.pause(0.1)

    if save_fig:

        plt.figure(2)
        plt.savefig(str(save_fig) + '_location.png')
        #plt.close('all')
    if display == True:
        plt.close('all')

def plot_time(result,display=False,save_fig=None,grnd_trth=None):

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
         plt.savefig(str(save_fig) + '_source_time.png')
         plt.close('all')

     if display==True:
         plt.close('all')
