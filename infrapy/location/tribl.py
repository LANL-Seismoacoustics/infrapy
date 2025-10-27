# infrapy.location.projection.py
#
# Back projection localization methods using the infraGA ray 
# tracing with auxiliary parameters to map direction-of-arrival
# (DOA) confidence into spatial and temporal confidence.
#
# Author            Philip Blom (pblom@lanl.gov)

import os
import fnmatch
import warnings
import subprocess
import tempfile
import json

from importlib.util import find_spec

if find_spec('infraga'):
    import wget
    from netCDF4 import Dataset

import numpy as np

from datetime import datetime
from pyproj import Geod

from scipy.integrate import simps
from scipy.interpolate import interp1d, interp2d

from . import bisl
from ..utils import prog_bar

sph_proj = Geod(ellps='sphere')
resol = '100m'  # use data at this scale (not working at the moment)

# ############################ #
#       Back Projection        #
#     Localization Methods     #
# ############################ #
def interp_etopo(ll_corner, ur_corner):
    etopo1 = Dataset(find_spec('infraga').submodule_search_locations[0] + "/resources/ETOPO1_Ice_g_gmt4.grd")

    grid_lons = etopo1.variables['x'][:]
    grid_lats = etopo1.variables['y'][:]
    grid_elev = etopo1.variables['z'][:]

    lat_mask = np.logical_and(ll_corner[0] - 2.0 <= grid_lats, grid_lats <= ur_corner[0] + 2.0).nonzero()[0]
    lon_mask = np.logical_and(ll_corner[1] - 2.0 <= grid_lons, grid_lons <= ur_corner[1] + 2.0).nonzero()[0]

    region_lat = grid_lats[lat_mask]
    region_lon = grid_lons[lon_mask]
    region_elev = grid_elev[lat_mask,:][:,lon_mask]

    # Change underwater values to sea surface
    region_elev[region_elev < 0.0] = 0.0

    return interp2d(region_lon, region_lat, region_elev / 1000.0, kind='linear')


def _compute_projections(det_list, atmo_file, temp_dest, grnd_snd_spd=None, latlon_bnds=None, bounces=100, cpu_cnt=None):
    lat_vals = [det.latitude for det in det_list]
    lon_vals = [det.longitude for det in det_list]

    if os.path.isfile(find_spec('infraga').submodule_search_locations[0] + "/resources/ETOPO1_Ice_g_gmt4.grd"):
        topo = interp_etopo([min(lat_vals), min(lon_vals)], [max(lat_vals), max(lon_vals)])
    else:
        print("Topography file not found.  Downloading from https://www.ngdc.noaa.gov/mgg/global/")
        download_url = "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"
        destination = find_spec('infraga').submodule_search_locations[0] + "/resources/ETOPO1_Ice_g_gmt4.grd.gz"
        try:
            if not os.path.isdir(os.path.split(destination)[0]):
                os.mkdir(os.path.split(destination)[0])
            
            print("Downloading ETOPO1 data...")
            wget.download(download_url, destination)
            print("Extracting...")
            os.system("gzip -d " + destination)
            print("ETOPO file successfully downloaded.")

            topo = interp_etopo([min(lat_vals), min(lon_vals)], [max(lat_vals), max(lon_vals)])

        except:
            print("Download failed.")
            print("Try manual download: " + download_url)
            print("Place extracted .grd file in " + find_spec('infraga').submodule_search_locations[0] + "/resources/")

    rcvr_elevs = np.array([topo(det.latitude, det.longitude)[0] for det in det_list])

    if grnd_snd_spd is None:
        # Compute sound speed from atmo file
        atmo = np.loadtxt(atmo_file)
        snd_spd = interp1d(atmo[:, 0], np.sqrt(0.14 * atmo[:, 5] / atmo[:, 4]))
        grnd_snd_spd = np.array([snd_spd(z_val) for z_val in rcvr_elevs])
    
    if len(grnd_snd_spd) == 1:
        grnd_snd_spd = [grnd_snd_spd[0]] * len(det_list)

    if len(grnd_snd_spd) != len(det_list):
        print('\t' + "Warning! Specificed grnd_snd_spd values don't match length of detections list.")
        return None

    command_list = []
    for n, det in enumerate(det_list):
        corners = np.meshgrid(latlon_bnds[0], latlon_bnds[1])
        _, _, temp = sph_proj.inv([det.longitude] * 4, [det.latitude] * 4, corners[1].flatten(), corners[0].flatten())
        max_rng = max(temp / 1000.0) * 1.2

        command = find_spec('infraga').submodule_search_locations[0] + "/bin/infraga-sph -back_proj " + atmo_file + " rcvr_lat=" + str(det.latitude) + " rcvr_lon=" + str(det.longitude)
        command = command + " azimuth=" + str(det.back_azimuth) + " inclination=" + str(np.degrees(np.arccos(min(grnd_snd_spd[n] / det.trace_velocity, 1.0))))
        command = command + " max_rng=" + str(max_rng) + " bounces=" + str(bounces) + " z_grnd=" + str(rcvr_elevs[n])
        command = command + " output_id=" + temp_dest + ".det-" + str(n) + " > /dev/null"
        
        command_list = command_list + [command]

    if cpu_cnt is not None:
        for j in range(0, len(command_list), cpu_cnt):
            procs_list = [subprocess.Popen(cmd, shell=True) for cmd in command_list[j:j + cpu_cnt]]
            for proc in procs_list:
                proc.communicate()
                proc.wait()
    else:
        procs_list = [subprocess.Popen(cmd, shell=True) for cmd in command_list]
        for proc in procs_list:
            proc.communicate()
            proc.wait()

    return grnd_snd_spd


class BackProjection(object):

    r_earth = 6370.0

    def __init__(self, detection, projection_file, det_time_std_dev=5.0, c0=340.0, c0_stdev=2.0, dt=1.0, az_limit=2.0):

        self.c0 = c0
        self.c0_stdev = c0_stdev

        # Develop method to estimate azimuth and inclination standard deviations from f-stat?
        v0 = detection.trace_velocity

        self.az_std_dev = np.degrees(1.0 / np.sqrt(2.0 * (detection.array_dim - 1.0) * detection.peakF_value))
        self.az_std_dev = max(self.az_std_dev, az_limit)


        self.tr_vel_std_dev = v0 * np.radians(self.az_std_dev)

        v0_up = v0 + v0 * np.radians(self.az_std_dev)
        v0_dn = v0 - v0 * np.radians(self.az_std_dev)

        incl_up = np.degrees(np.arccos(c0 / v0_up))
        incl_dn = np.degrees(np.arccos(min(1.0, c0 / v0_dn)))

        self.incl_std_dev = (incl_up - incl_dn) / 2.0

        incl_std_dev2 = self.c0_stdev**2 + (self.c0 / v0)**2 * self.tr_vel_std_dev**2
        incl_std_dev2 = self.incl_std_dev / (v0**2 - self.c0**2)
        incl_std_dev2 = np.degrees(np.sqrt(self.incl_std_dev))

        self.incl_std_dev = min(incl_std_dev2, self.incl_std_dev)

        self.det_time = np.datetime64(detection.peakF_UTCtime)

        # Read in projection and interpolate to resample
        projection = np.loadtxt(projection_file)

        self.tms = np.arange(projection[0][3], projection[-1][3], dt)

        self.lat = interp1d(projection[:, 3], projection[:, 0])(self.tms)
        self.lon = interp1d(projection[:, 3], projection[:, 1])(self.tms)
        self.alt = interp1d(projection[:, 3], projection[:, 2])(self.tms)

        self.sd_lat = interp1d(projection[:, 3], np.sqrt((projection[:, 4] * self.incl_std_dev)**2 + (projection[:, 8] * self.az_std_dev)**2))(self.tms)
        self.sd_lon = interp1d(projection[:, 3], np.sqrt((projection[:, 5] * self.incl_std_dev)**2 + (projection[:, 9] * self.az_std_dev)**2))(self.tms)
        self.sd_alt = interp1d(projection[:, 3], np.sqrt((projection[:, 6] * self.incl_std_dev)**2 + (projection[:, 10] * self.az_std_dev)**2))(self.tms)
        self.sd_tm = interp1d(projection[:, 3], np.sqrt((projection[:, 7] * self.incl_std_dev)**2 + (projection[:, 11] * self.az_std_dev)**2 + det_time_std_dev**2))(self.tms)

        self.norm = 1.0 / (4.0 * np.pi**2 * (self.sd_lat * self.sd_lon * self.sd_alt * self.sd_tm))

    def likelihood(self, lat0, lon0, alt0, t0, prog_step=0):
        t0 = np.atleast_1d(t0)
        dt = np.array([np.timedelta64(self.det_time - np.datetime64(tn)).astype('m8[ms]').astype(float) / 1.0e3 for tn in t0])

        result = np.array([self.norm[n] * np.exp(-1.0 / 2.0 * (((lat0 - self.lat[n]) / self.sd_lat[n])**2 + ((lon0 - self.lon[n]) / self.sd_lon[n])**2 
                                                                + ((alt0 - self.alt[n]) / self.sd_alt[n])**2 + ((dt - self.tms[n]) / self.sd_tm[n])**2)) for n in range(len(self.norm))])

        prog_bar.increment(n=prog_step)
        return np.sum(result, axis=0)
    

def build_projections(dets_list, atmo_file, projection_path, grnd_snd_spd=None, latlon_bnds=None, cpu_cnt=None, c0_stdev=2.5, det_time_std_dev=5.0, az_limit=2.0):

    c0 = _compute_projections(dets_list, atmo_file, temp_dest=projection_path, grnd_snd_spd=grnd_snd_spd, latlon_bnds=latlon_bnds, cpu_cnt=cpu_cnt)
    if c0 is not None:
        return [BackProjection(det, projection_path + ".det-" + str(n) + ".projection.dat", det_time_std_dev=det_time_std_dev, c0=c0[n], c0_stdev=c0_stdev, az_limit=az_limit) for n, det in enumerate(dets_list)]
    else:
        return None


def eval_on_grid(proj, lat_grid, lon_grid, alt_grid, tm_grid, prog_step):
    return proj.likelihood(lat_grid.flatten(), lon_grid.flatten(), alt_grid.flatten(), tm_grid.flatten(), prog_step=prog_step)

def eval_on_grid_wrapper(args):
    return eval_on_grid(*args)


def run(det_list, atmo_file, temp_path, bm_width=10.0, rng_max=2000.0, grid_resol=50, ll_corner=None, ur_corner=None, latlon_resol=None, tm_lims=None, tm_resol=None, alt_lims=None, alt_resol=1.0,
            grnd_snd_spd=340.0, c0_stdev=10.0, det_time_stdev=10.0, az_limit=2.0, verbose=True, show_prog=True, pool=None):

    if verbose:
        print("Running Time-Reversed Infrasonic Bayesian Localization (TRIBL) Analysis...")
        print('\t' + "Identifying integration region and building grid...")
    
    if alt_lims is None:
        alt_lims = [0.0, 0.0]
        alt_resol = 1.0

    lat_grid, lon_grid, alt_grid, tm_grid = bisl.build_grid(det_list, bm_width=bm_width, rng_max=rng_max, grid_resol=grid_resol, ll_corner=ll_corner, ur_corner=ur_corner,
                                                    latlon_resol=latlon_resol, include_tms=True, tm_lims=tm_lims, tm_resol=tm_resol, alt_lims=alt_lims, alt_resol=alt_resol)

    lat_vals = np.sort(np.unique(lat_grid))
    lon_vals = np.sort(np.unique(lon_grid))
    alt_vals = np.sqrt(np.unique(alt_grid))
    tm_vals = np.sort(np.unique(tm_grid))

    if verbose:
        print('\t' + "Computing back projections for detection list...")

    if pool is not None:
        cpu_cnt = pool._processes
    else:
        cpu_cnt = None

    projs = build_projections(det_list, atmo_file, temp_path, grnd_snd_spd=grnd_snd_spd, latlon_bnds=[[lat_vals[0], lat_vals[-1]], [lon_vals[0], lon_vals[-1]]], cpu_cnt=cpu_cnt, c0_stdev=c0_stdev, det_time_std_dev=det_time_stdev, az_limit=az_limit)

    if verbose:
        print('\t' + "Evaluating localization probability on grid...")
        print('\t\t Progress: ', end='')

    if show_prog or verbose:
        prog_bar.prep(5 * len(det_list))
        if pool:
            det_pdfs = pool.map(eval_on_grid_wrapper, [[proj, lat_grid, lon_grid, alt_grid, tm_grid, 5] for proj in projs])
        else:
            det_pdfs = np.array([eval_on_grid(proj, lat_grid, lon_grid, alt_grid, tm_grid, prog_step=5) for proj in projs])
        prog_bar.close()
    else:   
        if pool:
            det_pdfs = pool.map(eval_on_grid_wrapper, [[proj, lat_grid, lon_grid, alt_grid, tm_grid, 0] for proj in projs])
        else:
            det_pdfs = np.array([eval_on_grid(proj, lat_grid, lon_grid, alt_grid, tm_grid, prog_step=0) for proj in projs])

    pdf = np.prod(det_pdfs, axis=0)    
    pdf = pdf.reshape(lat_grid.shape)

    np.savez_compressed(temp_path + ".pdf", lat_vals=lat_vals, lon_vals=lon_vals, alt_vals=alt_vals, tm_vals=tm_vals, pdf=pdf)

    if np.max(pdf) > 0.0:
        result = bisl.analyze_pdf(pdf, lat_grid, lon_grid, tm_grid, verbose=verbose)
    else:
        result = {'norm' : 0.0}

        print('\nOne of the detections is driving the PDF to zero...')
        print('\tIndex\tmax(PDF)')
        for j, det in enumerate(det_pdfs):
            print('\t' + str(j) + '\t' + str(np.max(det)))

        print("Once it's working, try --det-mask ", np.where(np.max(det) > 0.0))

    return result



