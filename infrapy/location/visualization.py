# visualization.py
#
# Visualization methods for beamforming (fk) and detection (fd) results
#
# Philip Blom (pblom@lanl.gov)


import numpy as np

from pyproj import Geod, transform

import matplotlib.pyplot as plt 
from matplotlib import cm

import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.ticker as mticker

from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from infrapy.utils import data_io

sph_proj = Geod(ellps='sphere')


marker_size = 1.0
map_proj = crs.PlateCarree()
resol = '100m'  # use data at this scale (not working at the moment)

def plot_dets_on_map(det_list, range_max=1000.0, title=None, output_path=None, show_fig=True):
    '''
    Visualize detections on a Cartopy map

    '''   



    array_lats = np.array([det.latitude for det in det_list])
    array_lons = np.array([det.longitude for det in det_list])

    lat_min, lat_max = min(array_lats), max(array_lats)
    lon_min, lon_max = min(array_lons), max(array_lons)

    for det in det_list:
        gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=2, del_s=(range_max * 1.0e3 / 2))
        lat_min, lat_max = min(min(gc_path.lats), lat_min), max(max(gc_path.lats), lat_max)
        lon_min, lon_max = min(min(gc_path.lons), lon_min), max(max(gc_path.lons), lon_max)

    lat_min, lat_max = np.floor(lat_min), np.ceil(lat_max)
    lon_min, lon_max = np.floor(lon_min), np.ceil(lon_max)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=map_proj)

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)

    gl = ax.gridlines(crs=map_proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    lat_tick, lon_tick = int((lat_max - lat_min) / 4), int((lon_max - lon_min) / 4)
    gl.xlocator = mticker.FixedLocator(np.arange(lon_min - np.ceil(lon_tick / 2), lon_max + lon_tick, lon_tick))
    gl.ylocator = mticker.FixedLocator(np.arange(lat_min - np.ceil(lat_tick / 2), lat_max + lat_tick, lat_tick))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Add features (coast lines, borders)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    
    ax.plot(array_lons, array_lats, 'k^', markersize=7.5, transform=map_proj)

    for det in det_list:
        gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=100, del_s=(range_max * 1.0e3 / 100))
        ax.plot(list(gc_path.lons), list(gc_path.lats), '-b', linewidth=1.5, transform=map_proj)

    if title:
        plt.title(title)

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_loc(det_list, bisl_result, range_max=1000.0, zoom=False, title=None, output_path=None, show_fig=True):
    '''
    Visualize detections on a Cartopy map along with BISL analysis result
    
    '''   

    bisl_results = data_io.read_locs(bisl_result)
    spatial_pdf = np.array(bisl_results['spatial_pdf'])
    origin_times = np.array([np.datetime64(tn) for tn in bisl_results['temporal_pdf'][0]])
    origin_time_pdf = np.array(bisl_results['temporal_pdf'][1])

    if not zoom:
        array_lats = np.array([det.latitude for det in det_list])
        array_lons = np.array([det.longitude for det in det_list])

        lat_min, lat_max = min(array_lats), max(array_lats)
        lon_min, lon_max = min(array_lons), max(array_lons)

        for det in det_list:
            gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=2, del_s=(range_max * 1.0e3 / 2))
            lat_min, lat_max = min(min(gc_path.lats), lat_min), max(max(gc_path.lats), lat_max)
            lon_min, lon_max = min(min(gc_path.lons), lon_min), max(max(gc_path.lons), lon_max)
    else:
        lat_min, lat_max = 30.0, 40.0
        lon_min, lon_max = -100.0, 100.0

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=map_proj)

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)

    gl = ax.gridlines(crs=map_proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    lat_tick, lon_tick = int((lat_max - lat_min) / 4), int((lon_max - lon_min) / 4)
    gl.xlocator = mticker.FixedLocator(np.arange(lon_min - np.ceil(lon_tick / 2), lon_max + lon_tick, lon_tick))
    gl.ylocator = mticker.FixedLocator(np.arange(lat_min - np.ceil(lat_tick / 2), lat_max + lat_tick, lat_tick))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    if title:
        plt.set_title(title)

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_origin_time(bisl_results, title=None, output_path=None, show_fig=True):
    '''
    Visualize the origin time PDF computed in BISL
    
    '''  
        
    origin_times = np.array([np.datetime64(tn) for tn in bisl_results['temporal_pdf'][0]])
    origin_time_pdf = np.array(bisl_results['temporal_pdf'][1])

    mask = np.logical_and(np.datetime64(bisl_results['t_min']) <= origin_times,
                            origin_times <= np.datetime64(bisl_results['t_max']))

    plt.plot(origin_times, origin_time_pdf, '-k')
    plt.fill_between(origin_times[mask], 0.0, origin_time_pdf[mask], color='b', alpha=0.5)
    plt.xlabel("Origin Time")
    plt.ylabel("Probability")

    if title:
        plt.title(title)

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()

