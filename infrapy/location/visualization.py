# visualization.py
#
# Visualization methods for beamforming (fk) and detection (fd) results
#
# Philip Blom (pblom@lanl.gov)


import numpy as np

from pyproj import Geod

import matplotlib.pyplot as plt 
from matplotlib import cm

import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.ticker as mticker

import cartopy
import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from . import bisl

sph_proj = Geod(ellps='sphere')

marker_size = 1.0
map_proj = crs.PlateCarree()

# Need to decide on a color combination
back_az_color = "DarkRed"
conf_color = "r"
pdf_cm = cm.hot_r

def use_offline_maps(self, pre_existing_data_dir, turn_on=True):
    # call this function to initialize the use of offline maps.  turn_on will initialize the pre_existing_data_directory
    if turn_on:
        cartopy.config['pre_existing_data_dir'] = pre_existing_data_dir
    else:
        cartopy.config['pre_existing_data_dir'] = ""

def _setup_map(fig, latlon_bnds):
    lat_min, lat_max = latlon_bnds[0]
    lon_min, lon_max = latlon_bnds[1]

    ax = fig.add_subplot(1, 1, 1, projection=map_proj)
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)

    gl = ax.gridlines(crs=map_proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    lat_tick, lon_tick = max(1, int((lat_max - lat_min) / 4)), max(1, int((lon_max - lon_min) / 4))
    while len(np.arange(lat_min, lat_max, lat_tick)) < 3:
        lat_tick = lat_tick / 2.0
    while len(np.arange(lon_min, lon_max, lon_tick)) < 3:
        lon_tick = lon_tick / 2.0

    gl.xlocator = mticker.FixedLocator(np.arange(lon_min - np.ceil(lon_tick), lon_max + np.ceil(lon_tick), lon_tick))
    gl.ylocator = mticker.FixedLocator(np.arange(lat_min - np.ceil(lat_tick), lat_max + np.ceil(lat_tick), lat_tick))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Add features (coast lines, borders)
    ax.add_feature(cfeature.BORDERS, linewidth=1.0)
    ax.add_feature(cfeature.STATES, linewidth=1.0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.75)
    ax.add_feature(cfeature.LAKES, linewidth=0.5)
    ax.add_feature(cfeature.RIVERS, linewidth=0.5)

    return ax

def plot_dets_on_map(det_list, range_max=1000.0, title=None, output_path=None, show_fig=True):
    '''
    Visualize detections on a Cartopy map

    '''   

    array_lats = np.array([det.latitude for det in det_list])
    array_lons = np.array([det.longitude for det in det_list])

    lat_min, lat_max = min(array_lats), max(array_lats)
    lon_min, lon_max = min(array_lons), max(array_lons)

    for det in det_list:
        if det.back_azimuth is not None:
            gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=2, del_s=(range_max * 1.0e3 / 2))
            lat_min, lat_max = min(min(gc_path.lats), lat_min), max(max(gc_path.lats), lat_max)
            lon_min, lon_max = min(min(gc_path.lons), lon_min), max(max(gc_path.lons), lon_max)

    lat_min, lat_max = np.floor(lat_min), np.ceil(lat_max)
    lon_min, lon_max = np.floor(lon_min), np.ceil(lon_max)

    fig = plt.figure()
    ax = _setup_map(fig,[[lat_min, lat_max], [lon_min, lon_max]])

    for det in det_list:
        if det.back_azimuth is not None:
            gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=500, del_s=(range_max * 1.0e3 / 500))
            ax.plot(list(gc_path.lons), list(gc_path.lats), '.', color=back_az_color, markersize=1.5, transform=map_proj)
    ax.plot(array_lons, array_lats, 'k^', markersize=7.5, transform=map_proj)

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

    array_lats = np.array([det.latitude for det in det_list])
    array_lons = np.array([det.longitude for det in det_list])

    conf_x, conf_y = bisl.calc_conf_ellipse([0.0, 0.0],[bisl_result['EW_stdev'], bisl_result['NS_stdev'], bisl_result['covar']], 90)
    conf_latlon = sph_proj.fwd(np.array([bisl_result['lon_mean']] * len(conf_x)), np.array([bisl_result['lat_mean']] * len(conf_x)), np.degrees(np.arctan2(conf_x, conf_y)), np.sqrt(conf_x**2 + conf_y**2) * 1e3)

    if zoom:
        lat_min, lat_max = min(conf_latlon[1]), max(conf_latlon[1])
        lon_min, lon_max = min(conf_latlon[0]), max(conf_latlon[0])
    else:
        lat_min, lat_max = min(array_lats), max(array_lats)
        lon_min, lon_max = min(array_lons), max(array_lons)

        for det in det_list:
            if det.back_azimuth is not None:
                gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=2, del_s=(range_max * 1.0e3 / 2))
                lat_min, lat_max = min(min(gc_path.lats), lat_min), max(max(gc_path.lats), lat_max)
                lon_min, lon_max = min(min(gc_path.lons), lon_min), max(max(gc_path.lons), lon_max)

    lat_min, lat_max = np.floor(lat_min), np.ceil(lat_max)
    lon_min, lon_max = np.floor(lon_min), np.ceil(lon_max)

    fig = plt.figure()
    ax = _setup_map(fig, [[lat_min, lat_max], [lon_min, lon_max]])
    
    spatial_pdf = np.array(bisl_result['spatial_pdf'])
    ax.scatter(spatial_pdf[0].flatten(), spatial_pdf[1].flatten(), c=spatial_pdf[2].flatten(), marker="s", s=2.5, cmap=pdf_cm, transform=map_proj, alpha=0.5, edgecolor='none', vmin=0.0)
    ax.plot(conf_latlon[0], conf_latlon[1], color=conf_color, linewidth=1.5, transform=map_proj)

    if not zoom:
        for det in det_list:
            if det.back_azimuth is not None:
                gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=500, del_s=(range_max * 1.0e3 / 500))
                ax.plot(list(gc_path.lons), list(gc_path.lats), '.', color=back_az_color, markersize=1.5, transform=map_proj)
        ax.plot(array_lons, array_lats, 'k^', markersize=7.5, transform=map_proj)

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

    plt.figure()
    plt.plot(origin_times, origin_time_pdf, '-k')
    plt.fill_between(origin_times[mask], 0.0, origin_time_pdf[mask], color=conf_color, alpha=0.5)
    plt.xlabel("Origin Time")
    plt.ylabel("Probability")

    if title:
        plt.title(title)

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_spye(spye_result, title=None, output_path=None, show_fig=True):
    '''
    Visualize the yield estimate PDF from SpYE

    '''

    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3), gridspec_kw={'width_ratios': [2, 1]})

    ax1.set_xscale('log')
    ax1.plot(np.array(spye_result['yld_vals']), spye_result['yld_pdf'], '-k')
    ax1.fill_between(np.array(spye_result['yld_vals']), spye_result['yld_pdf'], where=np.logical_and(spye_result['conf_bnds'][0][0] <= np.array(spye_result['yld_vals']), np.array(spye_result['yld_vals']) <= spye_result['conf_bnds'][0][1]), color=back_az_color, alpha=0.25)
    ax1.fill_between(np.array(spye_result['yld_vals']), spye_result['yld_pdf'], where=np.logical_and(spye_result['conf_bnds'][1][0] <= np.array(spye_result['yld_vals']), np.array(spye_result['yld_vals']) <= spye_result['conf_bnds'][1][1]), color=conf_color, alpha=0.25)

    ax1.set_xlabel("Yield (eq. TNT) [tons]")
    ax1.set_ylabel("Probability")


    F, SP = np.meshgrid(spye_result['spec_freqs'], spye_result['spec_vals'])
    ax2.scatter(F.flatten(), SP.flatten(), c=np.array(spye_result['spec_pdf']).flatten(), cmap=pdf_cm)

    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Spectral Amplitude [Pa/Hz]")

    plt.tight_layout()

    if title:
        plt.title(title)
        
    if output_path:
        plt.savefig(output_path + ".spye.png", dpi=300) 
    
    if show_fig:
        plt.show()

