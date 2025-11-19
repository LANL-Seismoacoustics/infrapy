# visualization.py
#
# Visualization methods for beamforming (fk) and detection (fd) results
#
# Philip Blom (pblom@lanl.gov)


import numpy as np

from pyproj import Geod

from scipy.interpolate import interp1d
from scipy.stats import chi2

import matplotlib.pyplot as plt 
from matplotlib import cm

import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import matplotlib.dates as mdates

import cartopy
import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from . import bisl

sph_proj = Geod(ellps='sphere')

marker_size = 1.0
map_proj = crs.PlateCarree()

# Need to decide on a color combination
'''
back_az_color = "DarkRed"
conf_color = "r"
pdf_cm = cm.hot_r
'''

back_az_color = "darkgreen"
conf_color = "navy"
pdf_cm = cm.ocean_r


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
    ax.add_feature(cfeature.LAKES, linewidth=0.5, alpha=0.5)
    ax.add_feature(cfeature.RIVERS, linewidth=0.5)

    return ax


def _setup_map_ax(ax, latlon_bnds):
    lat_min, lat_max = latlon_bnds[0]
    lon_min, lon_max = latlon_bnds[1]

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
    ax.add_feature(cfeature.LAKES, linewidth=0.5, alpha=0.5)
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


def plot_ev_wvfrms(ev_data, use_loc=False, loc_index=0, use_gt=False):

    sta_list, sta_indices = np.unique([det['wvfrm_info'][0][0]['trace id'] for det in ev_data['det_info']], return_index=True)
    sta_cnt = len(sta_list)

    if use_loc:
        print("Plotting using localization result for ranges...")

        sta_locs = np.array([[ev_data['det_info'][k]['wvfrm_info'][0][0]['latitude'],  ev_data['det_info'][k]['wvfrm_info'][0][0]['longitude']] for k in sta_indices])
        src_loc = [ev_data['location'][loc_index]['result']['lat_mean'],ev_data['location'][loc_index]['result']['lon_mean']]
        sta_rngs = sph_proj.inv([src_loc[1]] * sta_cnt, [src_loc[0]] * sta_cnt, sta_locs[:, 1], sta_locs[:, 0], return_back_azimuth=True, radians=False)[2] / 1000.0

        fig = plt.figure(figsize=(8, 12))

        scaling = 25.0

        for j, sta in enumerate(sta_list):
            for det_k in ev_data['det_info']:
                if det_k['wvfrm_info'][0][0]['trace id'] == sta:
                    t0 = np.datetime64(det_k["peak f-stat time"])
                    l = np.argmax([max(fk_l['f-stat']) for fk_l in det_k['fk']])

                    t_vals = [t0 + np.timedelta64(int(dt * 1000.0), 'ms') for dt in det_k["beam"][l]["time"]]
                    plt.plot(t_vals, sta_rngs[j] + np.array(det_k["beam"][l]['signal']) * scaling, 'k', linewidth=0.5)

    elif use_gt:
        print("Plotting using ground truth location for ranges...")
    


    else:
        print("Plotting event waveforms station-by-station...")

        fig = plt.figure(figsize=(8, 1 + 1.5 * sta_cnt), layout="constrained")
        spec = fig.add_gridspec(sta_cnt, 1)

        ax0 = fig.add_subplot(spec[sta_cnt - 1])

        ax0.set_xlabel("")
        ax0.tick_params(axis='x', labelrotation=30)
        ax0.set_ylabel("Pressure [Pa]")

        t_min = min([np.datetime64(det_k["peak f-stat time"]) + np.timedelta64(int(det_k["beam"][0]["time"][0] * 1000.0), 'ms') for det_k in ev_data['det_info']])
        t_max = max([np.datetime64(det_k["peak f-stat time"]) + np.timedelta64(int(det_k["beam"][0]["time"][-1] * 1000.0), 'ms') for det_k in ev_data['det_info']])
        ax0.set_xlim((t_min, t_max))

        for det_k in ev_data['det_info']:
            if det_k['wvfrm_info'][0][0]['trace id'] == sta_list[0]:
                t0 = np.datetime64(det_k["peak f-stat time"])
                l = np.argmax([max(fk_l['f-stat']) for fk_l in det_k['fk']])

                t_vals = [t0 + np.timedelta64(int(dt * 1000.0), 'ms') for dt in det_k["beam"][l]["time"]]
                ax0.plot(t_vals, det_k["beam"][l]['signal'], 'k', linewidth=0.5)

        ax0.annotate(sta_list[0], (0.975, 0.95), xycoords='axes fraction', horizontalalignment = "right", verticalalignment="top")

        for k in range(1, sta_cnt):
            ax_k = fig.add_subplot(spec[sta_cnt - (k + 1)], sharex=ax0)
            ax_k.set_ylabel("Pressure [Pa]")
            plt.setp(ax_k.get_xticklabels(), visible=False)

            for det_k in ev_data['det_info']:
                if det_k['wvfrm_info'][0][0]['trace id'] == sta_list[k]:
                    t0 = np.datetime64(det_k["peak f-stat time"])
                    l = np.argmax([max(fk_l['f-stat']) for fk_l in det_k['fk']])

                    t_vals = [t0 + np.timedelta64(int(dt * 1000.0), 'ms') for dt in det_k["beam"][l]["time"]]
                    ax_k.plot(t_vals, det_k["beam"][l]['signal'], 'k', linewidth=0.5)

            ax_k.annotate(sta_list[k], (0.975, 0.95), xycoords='axes fraction', horizontalalignment = "right", verticalalignment="top")

    plt.show()





def plot_localization(det_list, loc_dict, grnd_truth_dict, confidence_level=90.0, range_max=None, output_path=None, show_fig=True):

    loc_params = loc_dict['params']
    loc_result = loc_dict['result']

    array_lats = np.array([det.latitude for det in det_list])
    array_lons = np.array([det.longitude for det in det_list])

    conf_x, conf_y = bisl.calc_conf_ellipse([0.0, 0.0],[loc_result['EW_stdev'], loc_result['NS_stdev'], loc_result['covar']], confidence_level)
    conf_latlon = sph_proj.fwd(np.array([loc_result['lon_mean']] * len(conf_x)), np.array([loc_result['lat_mean']] * len(conf_x)), np.degrees(np.arctan2(conf_x, conf_y)), np.sqrt(conf_x**2 + conf_y**2) * 1e3)

    lat_min, lat_max = min(array_lats), max(array_lats)
    lon_min, lon_max = min(array_lons), max(array_lons)

    for det in det_list:
        if det.back_azimuth is not None:
            gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=2, del_s=(range_max * 1.0e3 / 2), return_back_azimuth=False)
            lat_min, lat_max = min(min(gc_path.lats), lat_min), max(max(gc_path.lats), lat_max)
            lon_min, lon_max = min(min(gc_path.lons), lon_min), max(max(gc_path.lons), lon_max)

    lat_min, lat_max = np.floor(lat_min), np.ceil(lat_max)
    lon_min, lon_max = np.floor(lon_min), np.ceil(lon_max)

    lat_mean = (lat_max + lat_min) / 2.0
    ratio = (lat_max - lat_min) / (lon_max - lon_min) * np.cos(np.radians(lat_mean))

    fig = plt.figure(figsize=(10, 10 * ratio), layout="constrained")
    spec = fig.add_gridspec(4, 5)

    ax = fig.add_subplot(spec[:3, :3], projection=map_proj)
    ax_tm = fig.add_subplot(spec[3:, :3])
    ax_an = fig.add_subplot(spec[:2, 3:])
    ax_zm = fig.add_subplot(spec[2:, 3:], projection=map_proj)

    # Zoomed out map
    _setup_map_ax(ax, [[lat_min, lat_max], [lon_min, lon_max]])
    spatial_pdf = np.array(loc_result['spatial_pdf'])
    ax.scatter(spatial_pdf[0].flatten(), spatial_pdf[1].flatten(), c=spatial_pdf[2].flatten(), marker="s", s=5.0, cmap=pdf_cm, transform=map_proj, alpha=0.5, edgecolor='none', vmin=0.0)
    ax.plot(conf_latlon[0], conf_latlon[1], color=conf_color, linewidth=1.5, transform=map_proj)

    for det in det_list:
        if det.back_azimuth is not None:
            gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=500, del_s=(range_max * 1.0e3 / 500), return_back_azimuth=False)
            ax.plot(list(gc_path.lons), list(gc_path.lats), '.', color=back_az_color, markersize=1.5, transform=map_proj)
    ax.plot(array_lons, array_lats, 'k^', markersize=7.5, transform=map_proj)

    if 'latitude' in grnd_truth_dict.keys():
        ax.plot([float(grnd_truth_dict['longitude'])], [float(grnd_truth_dict['latitude'])], '*r', markersize=5.0, transform=map_proj)

    # Zoomed in map
    lat_min, lat_max = np.floor(min(conf_latlon[1])), np.ceil(max(conf_latlon[1]))
    lon_min, lon_max = np.floor(min(conf_latlon[0])), np.ceil(max(conf_latlon[0]))

    if 'latitude' in grnd_truth_dict.keys():
        lat_min, lat_max = min(lat_min, float(grnd_truth_dict['latitude'])), max(lat_max, float(grnd_truth_dict['latitude']))
        lon_min, lon_max = min(lon_min, float(grnd_truth_dict['longitude'])), max(lon_max, float(grnd_truth_dict['longitude']))

    _setup_map_ax(ax_zm, [[lat_min, lat_max], [lon_min, lon_max]])

    ax_zm.scatter(spatial_pdf[0].flatten(), spatial_pdf[1].flatten(), c=spatial_pdf[2].flatten(), marker="s", s=5.0, cmap=pdf_cm, transform=map_proj, alpha=0.5, edgecolor='none', vmin=0.0)
    ax_zm.plot(conf_latlon[0], conf_latlon[1], color=conf_color, linewidth=1.5, transform=map_proj)

    if 'latitude' in grnd_truth_dict.keys():
        ax_zm.plot([float(grnd_truth_dict['longitude'])], [float(grnd_truth_dict['latitude'])], '*r', markersize=10.0, transform=map_proj)

    # Origin time
    dt_vals = np.array([(np.datetime64(tm_val) - np.datetime64(loc_result['temporal_pdf'][0][0])).astype('m8[ms]').astype(float) / 1.0e3 for tm_val in loc_result['temporal_pdf'][0]])
    dt_mean = (np.datetime64(loc_result['t_mean']) - np.datetime64(loc_result['temporal_pdf'][0][0])).astype('m8[ms]').astype(float) / 1.0e3
    tm_mask = np.logical_and(dt_mean - 5.0 * loc_result['t_stdev'] < dt_vals, dt_vals < dt_mean + 5.0 * loc_result['t_stdev'])
    if confidence_level != 90.0:
            tm_conf = bisl.find_confidence(interp1d(dt_vals[tm_mask], np.array(loc_result['temporal_pdf'][1])[tm_mask], kind='cubic'), [dt_vals[tm_mask][0], dt_vals[tm_mask][-1]], confidence_level / 100.0)            
            tm_min_val = str(np.datetime64(loc_result['temporal_pdf'][0][0]) + np.timedelta64(int(min(tm_conf[0]) * 1e3), 'ms'))
            tm_max_val = str(np.datetime64(loc_result['temporal_pdf'][0][0]) + np.timedelta64(int(max(tm_conf[0]) * 1e3), 'ms'))
    else:
        tm_min_val = loc_result['t_min']
        tm_max_val = loc_result['t_max']

    origin_times = np.array([np.datetime64(tn) for tn in loc_result['temporal_pdf'][0]])
    origin_time_pdf = np.array(loc_result['temporal_pdf'][1])

    conf_mask = np.logical_and(np.datetime64(tm_min_val) <= origin_times, origin_times <= np.datetime64(tm_max_val))
    
    ax_tm.plot(origin_times[tm_mask], origin_time_pdf[tm_mask], '-k', linewidth=2.5)
    ax_tm.fill_between(origin_times[conf_mask], 0.0, origin_time_pdf[conf_mask], color=conf_color, alpha=0.5)
    if 'orig_tm' in grnd_truth_dict.keys():
        ax_tm.axvline(x = np.datetime64(grnd_truth_dict['orig_tm']), color='red')

    ax_tm.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax_tm.set_ylim(0)
    ax_tm.yaxis.set_ticklabels([])

    # annotation
    ax_an.axis('off')   
    conf_area = np.round(np.pi * loc_result['NS_stdev'] * loc_result['EW_stdev'] * chi2(2).ppf(confidence_level / 100.0), 2)

    def round_str(val):
        return str(np.round(val, 2))

    if loc_params['atmo_data'] is None: 
        loc_summary = 'BISL Summary\n' + '-' * 20 + '\n'
    else:
        loc_summary = 'TRIBL Summary\n' + '-' * 26 + '\n'

    loc_summary = loc_summary + str(np.round(loc_result['lat_mean'], 4)) + " deg N +/- " + round_str(loc_result["NS_stdev"]) + " km," + '\n'
    loc_summary = loc_summary + str(np.round(loc_result['lon_mean'], 4)) + " deg E +/- " + round_str(loc_result["EW_stdev"]) + " km," + '\n'
    loc_summary = loc_summary + loc_result['t_mean'] + " +/- " + round_str(loc_result['t_stdev']) + " s" + '\n\n'
    loc_summary = loc_summary + str(confidence_level) + "% confidence area: " + str(conf_area) + " sq km" + '\n'
    loc_summary = loc_summary + str(confidence_level) + "% confidence origin time:" + '\n  ' + tm_min_val + '\n  ' + tm_max_val + '\n'

    if loc_params['atmo_data'] is None: 
        if loc_params['pgm_file'] is None:
            loc_summary = loc_summary + '\n' + "Celerity model: " + loc_params['celerity_model']
        else:
            loc_summary = loc_summary + '\n' + "PGM file: " + loc_params['pgm_file']
    else:
        loc_summary = loc_summary + '\n' "atmo_data: " + loc_params['atmo_data']

    ax_an.text(0.0, 1.0, loc_summary, va="top", fontsize=10, bbox=dict(boxstyle="round, pad=0.2", fc="lightsteelblue", ec="black", lw=1))

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()
        

def plot_loc(det_list, bisl_result, range_max=1000.0, zoom=False, title=None, output_path=None, grnd_truth=None, show_fig=True):
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
                gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=2, del_s=(range_max * 1.0e3 / 2), return_back_azimuth=False)
                lat_min, lat_max = min(min(gc_path.lats), lat_min), max(max(gc_path.lats), lat_max)
                lon_min, lon_max = min(min(gc_path.lons), lon_min), max(max(gc_path.lons), lon_max)

    lat_min, lat_max = np.floor(lat_min), np.ceil(lat_max)
    lon_min, lon_max = np.floor(lon_min), np.ceil(lon_max)

    fig = plt.figure()
    ax = _setup_map(fig, [[lat_min, lat_max], [lon_min, lon_max]])
    
    spatial_pdf = np.array(bisl_result['spatial_pdf'])
    ax.scatter(spatial_pdf[0].flatten(), spatial_pdf[1].flatten(), c=spatial_pdf[2].flatten(), marker="s", s=5.0, cmap=pdf_cm, transform=map_proj, alpha=0.5, edgecolor='none', vmin=0.0)
    ax.plot(conf_latlon[0], conf_latlon[1], color=conf_color, linewidth=1.5, transform=map_proj)

    if not zoom:
        for det in det_list:
            if det.back_azimuth is not None:
                gc_path = sph_proj.fwd_intermediate(det.longitude, det.latitude, det.back_azimuth, npts=500, del_s=(range_max * 1.0e3 / 500), return_back_azimuth=False)
                ax.plot(list(gc_path.lons), list(gc_path.lats), '.', color=back_az_color, markersize=1.5, transform=map_proj)
        ax.plot(array_lons, array_lats, 'k^', markersize=7.5, transform=map_proj)


    if grnd_truth is not None:
        plt.plot([grnd_truth[1]], [grnd_truth[0]], '*m', markersize=7.5)

    if title:
        plt.set_title(title)

    if output_path:
        plt.savefig(output_path, dpi=300) 

    if show_fig:
        plt.show()


def plot_origin_time(bisl_results, title=None, output_path=None, show_fig=True, grnd_truth=None):
    '''
    Visualize the origin time PDF computed in BISL
    
    '''  
        
    origin_times = np.array([np.datetime64(tn) for tn in bisl_results['temporal_pdf'][0]])
    origin_time_pdf = np.array(bisl_results['temporal_pdf'][1])

    mask = np.logical_and(np.datetime64(bisl_results['t_min']) <= origin_times,
                            origin_times <= np.datetime64(bisl_results['t_max']))

    plt.figure(figsize=(7, 3))
    plt.plot(origin_times, origin_time_pdf, '-k', linewidth=2.5)
    plt.fill_between(origin_times[mask], 0.0, origin_time_pdf[mask], color=conf_color, alpha=0.5)
    plt.xlabel("Origin Time")
    plt.ylabel("Probability")
    plt.ylim(0)

    if grnd_truth is not None:
        print('\t' + "Including ground truth origin time...")
        plt.axvline(np.datetime64(grnd_truth), color='m', linewidth=3.33)
        
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

