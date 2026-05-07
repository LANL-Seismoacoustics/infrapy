#!/usr/bin/env python

import os 
import sys
import fnmatch
import click
import tempfile

from importlib.util import find_spec

import configparser as cnfg
import numpy as np

from obspy import Stream

from scipy.integrate import simps

from multiprocessing import Pool

from ..location import bisl, tribl
from ..propagation import infrasound
from ..characterization import spye

from ..utils import config
from ..utils import data_io

@click.command('run_loc', short_help="Estimate source locations and times for events", hidden=True)
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-detect-label", help="Detection path and pattern", default=None)
@click.option("--local-loc-label", help="Localization results path", default=None)
@click.option("--back-az-width", help="Width of beam projection (default: " + config.defaults['LOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Max source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)

@click.option("--grid-resol", help="Grid resolution (number of points) (default: " + config.defaults['LOC']['grid_resol'] + ")", default=None, type=int)

@click.option("--ll-corner", help="Lower left corner of region (lat, lon)", default=None)
@click.option("--ur-corner", help="Upper right corner of region (lat, lon)", default=None)
@click.option("--latlon-resol", help="Resolution of latitude/longitude grid (degrees)", default=None, type=float)

@click.option("--tm-min", help="Minimum origin time", default=None)
@click.option("--tm-max", help="Maximum origin time", default=None)
@click.option("--tm-resol", help="Resolution of origin time grid (seconds)", default=None, type=float)

@click.option("--celerity-model", help="Use included celerity model (default: '" + config.defaults['LOC']['celerity_model'], default=None, hidden=True)
@click.option("--rcel-wts", help="Custom reciprocal celerity model weights", default=None, hidden=True)
@click.option("--rcel-mns", help="Custom reciprocal celerity model means", default=None, hidden=True)
@click.option("--rcel-sds", help="Custom reciprocal celerity model standard deviations", default=None, hidden=True)
@click.option("--pgm-file", help="Path geometry model (PGM) file (optional)", default=None)

# add TRIBL options here...
@click.option("--atmo-data", help="Atmosphere data if using TRIBL", default=None)
@click.option("--alt-lims", help="Altitude limits if using TRIBL", default=None)
@click.option("--alt-resol", help="Altitude resolution if using TRIBL", default=None)
@click.option("--grnd-snd-spd", help="Sound speed at the ground", default=None)
@click.option("--c0-stdev", help="Sound speed uncertainty if using TRIBL", default=None)
@click.option("--det-tm-stdev", help="Detection time uncertainty", default=None)
@click.option("--local-temp-dir", help="Local temporary directory if using TRIBL", default=None)
@click.option("--cpu-cnt", help="CPU count for multithreading (default: None)", default=None, type=int)
def run_loc(config_file, detect_label, local_loc_label, back_az_width, range_max, grid_resol, ll_corner, ur_corner, latlon_resol, tm_min, tm_max, tm_resol, celerity_model, rcel_wts, rcel_mns, rcel_sds, pgm_file, atmo_data, alt_lims, alt_resol, grnd_snd_spd, c0_stdev, det_tm_stdev, local_temp_dir, cpu_cnt):
    '''
    Run Bayesian Infrasonic Source Localization (BISL) methods to estimate the source location and origin time for an event

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_loc --local-detect-label GJI_example-ev0  --local-loc-label GJI_example-ev0
    \tinfrapy run_loc --local-detect-label data/detection_set2.json --local-loc-label data/location2 --pgm-file ../infrapy/propagation/priors/UTTR_models/UTTR_06_1800UTC.pgm
    \tinfrapy run_loc --local-detect-label data/Blom_etal2024_GJI/UTTR  --local-loc-label data/Blom_etal2024_GJI/UTTR-BISL
    \tinfrapy run_loc --config-file config/tribl_example.config
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##      Localization Analysis      ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    # Data IO parameters
    detect_label = config.set_param(user_config, 'DETECTION IO', 'detect_label', detect_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')

    if ".loc.json" in local_loc_label:
       local_loc_label = local_loc_label[:-9] 
            
    click.echo('\n' + "Data summary:")
    click.echo("  detect_label: " + str(detect_label))
    click.echo("  local_loc_label: " + str(local_loc_label))

    # Algorithm parameters
    back_az_width = config.set_param(user_config, 'LOC', 'back_az_width', back_az_width, 'float')
    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    grid_resol = config.set_param(user_config, 'LOC', 'grid_resol', grid_resol, 'int')       

    ll_corner = config.set_param(user_config, 'LOC', 'll_corner', ll_corner, 'str')
    ur_corner = config.set_param(user_config, 'LOC', 'ur_corner', ur_corner, 'str')
    latlon_resol = config.set_param(user_config, 'LOC', 'latlon_resol', latlon_resol, 'float')

    tm_min = config.set_param(user_config, 'LOC', 'tm_min', tm_min, 'str')
    tm_max = config.set_param(user_config, 'LOC', 'tm_max', tm_max, 'str')
    tm_resol = config.set_param(user_config, 'LOC', 'tm_resol', tm_resol, 'float')

    celerity_model = config.set_param(user_config, 'LOC', 'celerity_model', celerity_model, 'str')
    rcel_wts = config.set_param(user_config, 'LOC', 'rcel_wts', rcel_wts, 'str')
    rcel_mns = config.set_param(user_config, 'LOC', 'rcel_mns', rcel_mns, 'str')
    rcel_sds = config.set_param(user_config, 'LOC', 'rcel_sds', rcel_sds, 'str')
    pgm_file = config.set_param(user_config, 'LOC', 'pgm_file', pgm_file, 'str')

    atmo_data = config.set_param(user_config, 'LOC', 'atmo_data', atmo_data, 'str')
    alt_lims = config.set_param(user_config, 'LOC', 'alt_lims', alt_lims, 'str')
    alt_resol = config.set_param(user_config, 'LOC', 'alt_resol', alt_resol, 'float')
    grnd_snd_spd = config.set_param(user_config, 'LOC', 'grnd_snd_spd', grnd_snd_spd, 'float')
    c0_stdev = config.set_param(user_config, 'LOC', 'c0_stdev', c0_stdev, 'float')
    det_tm_stdev = config.set_param(user_config, 'LOC', 'det_tm_stdev', det_tm_stdev, 'float')
    local_temp_dir = config.set_param(user_config, 'LOC', 'local_temp_dir', local_temp_dir, 'str')
    cpu_cnt = config.set_param(user_config, 'ASSOC', 'cpu_cnt', cpu_cnt, 'int')

    # Summarie parameters
    click.echo('\n' + "Parameter summary:")
    if ll_corner is None:
        click.echo("  back_az_width: " + str(back_az_width))
        click.echo("  range_max: " + str(range_max))
    else:
        click.echo("  ll_corner: " + str(ll_corner))
        click.echo("  ur_corner: " + str(ur_corner))

        ll_corner = np.array([float(val) for val in ll_corner.replace(" ","").split(",")])
        ur_corner = np.array([float(val) for val in ur_corner.replace(" ","").split(",")])

    if latlon_resol is not None:
        click.echo("  latlon_resol: " + str(latlon_resol))
    else: 
        click.echo("  grid_resol: " + str(grid_resol))

    if tm_min is not None:
        click.echo("  tm_min: " + str(tm_min))
        click.echo("  tm_max: " + str(tm_max))
        tm_lims = (np.datetime64(tm_min), np.datetime64(tm_max))
    else:
        tm_lims = None 

    if tm_resol is not None:
        click.echo("  tm_resol: " + str(tm_resol))
    else:
        if latlon_resol is not None:
            click.echo("  grid_resol: " + str(grid_resol))

    if atmo_data is not None:
        click.echo("  atmo_data: " + str(atmo_data))
        if grnd_snd_spd is not None:
            click.echo("  grnd_snd_spd: " + str(grnd_snd_spd))

        click.echo("  c0_stdev: " + str(c0_stdev))
        click.echo("  det_tm_stdev: " + str(det_tm_stdev))

        if alt_lims is not None:
            click.echo("  alt_lims: " + str(alt_lims))
            click.echo("  alt_resol: " + str(alt_resol))

        if local_temp_dir is not None:
            click.echo("  local_temp_dir: " + str(local_temp_dir))
    else:
        if pgm_file is not None:
            click.echo("  pgm_file: " + str(pgm_file))
            pgm = infrasound.PathGeometryModel()
            pgm.load(pgm_file)
        else:
            infrasound.set_celerity_model(celerity_model, rcel_wts=rcel_wts, rcel_mns=rcel_mns, rcel_sds=rcel_sds)
            pgm = None  

    if cpu_cnt is not None:
        click.echo("  cpu_cnt: " + str(cpu_cnt))
        pl = Pool(cpu_cnt)
    else:
        pl = None


    click.echo("")
    events = data_io.set_det_list(detect_label, merge=False)

    if type(events[0]) is list:
        # run localization analysis for multiple detection sets
        for j, det_list in enumerate(events):
            click.echo('\n' + "Running BISL on event " + str(j + 1) + " of " + str(len(events)))
            result = bisl.run(det_list, bm_width=back_az_width, rng_max=range_max, grid_resol=grid_resol, ll_corner=ll_corner, ur_corner=ur_corner, latlon_resol=latlon_resol, tm_lims=tm_lims, tm_resol=tm_resol, path_geo_model=pgm)

            # Determine output format for BISL results
            click.echo('\n' + "BISL Summary:")
            click.echo(bisl.summarize(result))

            click.echo("Writing localization result into " + local_loc_label + "_ev-" + str(j + 1) + "loc.json")
            data_io.write_json(result, local_loc_label + "_ev-" + str(j + 1) + "loc.json")

    else:
        # run a single localization analysis
        if atmo_data is None:           
            result = bisl.run(events, bm_width=back_az_width, rng_max=range_max, grid_resol=grid_resol, ll_corner=ll_corner, ur_corner=ur_corner, 
                                    latlon_resol=latlon_resol, tm_lims=tm_lims, tm_resol=tm_resol, path_geo_model=pgm)
        else:
            if find_spec('infraga'):
                if not os.path.isfile(find_spec('infraga').submodule_search_locations[0] + "/bin/infraga-sph"):
                    click.echo("InfraGA methods not compiled.  Run 'infraga compile' and try again.")
                else:                               
                    with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
                        if local_temp_dir is not None:
                            if not os.path.isdir(local_temp_dir):
                                os.mkdir(local_temp_dir)
                            tmpdirname = local_temp_dir

                        temp_path = tmpdirname + "/temp"

                        if "*" in atmo_data:
                            if len(os.path.dirname(atmo_data)) > 0:
                                file_path = os.path.dirname(atmo_data) + "/"
                            else:
                                file_path = ""

                            if "/" in atmo_data:
                                dir_files = os.listdir(os.path.dirname(atmo_data))
                            else:
                                dir_files = os.listdir(".")

                            file_list = []
                            for file in dir_files:
                                if fnmatch.fnmatch(file, os.path.basename(atmo_data)):
                                    file_list += [file]
                            file_list = np.sort(file_list)

                            click.echo('\n' + "Computing localization using atmosphere ensemble:")
                            norms = []
                            for k, file_name in enumerate(file_list):                            
                                print('\t' + str(k + 1) + '/' + str(len(file_list)) + '\t' + file_path + file_name + '\t', end='')
                                temp = tribl.run(events, file_path + file_name, temp_path + "-" + str(k), bm_width=back_az_width, rng_max=range_max, grid_resol=grid_resol, ll_corner=ll_corner, ur_corner=ur_corner,
                                                latlon_resol=latlon_resol, tm_lims=tm_lims, tm_resol=tm_resol, alt_lims=alt_lims, alt_resol=alt_resol, grnd_snd_spd=grnd_snd_spd, c0_stdev=c0_stdev, det_time_stdev=det_tm_stdev, verbose=False, show_prog=True, pool=pl) 
                                norms = norms + [temp['norm']]

                            norms = norms / np.sum(norms)

                            print('\n' + "Merging PDFs across the ensemble...")
                            tmp_0 = np.load(temp_path + "-0.pdf.npz")

                            lat_grid, lon_grid, alt_grid, tm_grid = np.meshgrid(tmp_0['lat_vals'], tmp_0['lon_vals'], tmp_0['alt_vals'], tmp_0['tm_vals'], indexing='ij')
                            lat_grid = np.squeeze(lat_grid)
                            lon_grid = np.squeeze(lon_grid)
                            alt_grid = np.squeeze(alt_grid)
                            tm_grid = np.squeeze(tm_grid)

                            pdf = tmp_0['pdf']

                            print('\t1/' + str(len(file_list)) + '\t' + file_path + file_list[0] + '\t' + str(norms[0]))

                            for k, file_name in enumerate(file_list[1:]):                            
                                tmp_k = np.load(temp_path + "-" + str(k + 1) + ".pdf.npz")
                                pdf = pdf + tmp_k['pdf']
                                print('\t' + str(k + 2) + '/' + str(len(file_list)) + '\t' + file_path + file_name + '\t' + str(norms[k + 1]))
                        
                            click.echo('\n' + "Analyzing combined localization PDF...")
                            result = bisl.analyze_pdf(pdf, lat_grid, lon_grid, tm_grid, verbose=True)

                        else:              
                            result = tribl.run(events, atmo_data, temp_path, bm_width=back_az_width, rng_max=range_max, grid_resol=grid_resol, ll_corner=ll_corner, ur_corner=ur_corner,
                                                latlon_resol=latlon_resol, tm_lims=tm_lims, tm_resol=tm_resol, alt_lims=alt_lims, alt_resol=alt_resol, grnd_snd_spd=grnd_snd_spd, c0_stdev=c0_stdev, det_time_stdev=det_tm_stdev, verbose=True, pool=pl)
            else:
                click.echo('\n' + "Can't run TRIBL methods without infraGA installed for ray tracing")
                return

        # Determine output format for BISL results
        click.echo('\n' + "Localization Summary:")
        click.echo(bisl.summarize(result))

        if ".loc.json" not in local_loc_label:
            local_loc_label = local_loc_label + ".loc.json"
        click.echo("Writing localization result into " + local_loc_label)
        data_io.write_json(result, local_loc_label)

    if pl is not None:
        pl.terminate()
        pl.close()




@click.command('regional', short_help="Run analysis using a single set of TLMs")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--local-detect-label", help="Detection results path", default=None)
@click.option("--local-loc-label", help="Localization results path", default=None)
@click.option("--local-yld-label", help="Output file for results", default=None)
@click.option("--tlm-label", help="Transmission loss model (TLM) path", default=None)
@click.option("--src-lat", help="Source latitude (if no loc result file)", default=None)
@click.option("--src-lon", help="Source longitude (if no loc result file)", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['YIELD']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['YIELD']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--yld-min", help="Minimum yield (default: " + config.defaults['YIELD']['yld_min'] + " [tons eq. TNT])", default=1.0, type=float)
@click.option("--yld-max", help="Maximum yield (default: " + config.defaults['YIELD']['yld_max'] + " [tons eq. TNT])", default=1000.0, type=float)
@click.option("--ref-rng", help="Reference range for blastwave model (default " + config.defaults['YIELD']['ref_rng'] + " km)", default=1.0, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['YIELD']['resolution'] + ")", default=None, type=int)
@click.option("--noise-option", help="Noise option ('pre', 'post', or 'beam')", default=None)
@click.option("--window-buffer", help="Window buffer scaling (default: " + config.defaults['YIELD']['window_buffer'] + ")", default=None, type=float)
@click.option("--amb-press", help="Ambient pressure (default: " + config.defaults['YIELD']['amb_press'] + " [Pa])", default=None, type=float)
@click.option("--amb-temp", help="Ambient temperature (default: " + config.defaults['YIELD']['amb_temp'] + " [K])", default=None, type=float)
@click.option("--grnd-burst", help="Ground burst assumption (default: " + config.defaults['YIELD']['grnd_burst'] + " [Hz])", default=None, type=bool)
@click.option("--exp-type", help="Explosion type ('chemical' or 'nuclear')", default=None)

def regional(config_file, local_wvfrms, fdsn, db_config, detect_label, local_loc_label, local_yld_label, tlm_label, 
                src_lat, src_lon, freq_min, freq_max, yld_min, yld_max, ref_rng, resolution, noise_option, window_buffer,
                amb_press, amb_temp, grnd_burst, exp_type):
    '''
    Run Spectral Yield Estimation (SpYE) methods to estimate the equivalent TNT yield of an above-ground explosion using a single set of transmission loss models (TLMs)

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_spye regional --local-wvfrms '../infrapy-data/hrr-5/*/*.sac' --local-detect-label data/HRR-5.dets.json --src-lat 33.5377 --src-lon -106.333961 --tlm-label "../infrapy/propagation/priors/tloss/2007_08-" --local-yld-label "HRR-5"
    '''

    click.echo("")
    click.echo("########################################")
    click.echo("##                                    ##")
    click.echo("##              InfraPy               ##")
    click.echo("##  Spectral Yield Estimation (SpYE)  ##")
    click.echo("##   Regional (Single TLM) Analysis   ##")
    click.echo("##                                    ##")
    click.echo("########################################")
    click.echo("")    

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None    

    # Waveform info
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')  
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None
    
    # Data IO parameters
    detect_label = config.set_param(user_config, 'DETECTION IO', 'detect_label', detect_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')
    src_lat = config.set_param(user_config, 'YIELD', 'src_lat', src_lat, 'float')
    src_lon = config.set_param(user_config, 'YIELD', 'src_lon', src_lon, 'float')

    click.echo('\n' + "Data parameters:")
    click.echo("  detect_label: " + str(detect_label))
    click.echo("  tlm_label: " + str(tlm_label))
    click.echo("  local_loc_label: " + str(local_loc_label))
    if local_loc_label is not None:
        src_loc = [local_loc_label['lat_mean'], local_loc_label['lon_mean']]
        # src_loc = [local_loc_label['lat_MaP'], local_loc_label['lon_MaP']]
        click.echo("    src_loc (from [...].loc.json file):", src_loc)
    elif src_lat is not None and src_lon is not None:
        src_loc = [src_lat, src_lon]
        click.echo("    src_lat: " + str(src_lat))
        click.echo("    src_lon: " + str(src_lon))
    else:
        print("Error.  Yield analysis requires either localization analysis file or specified source location")
        return 0

    if local_wvfrms is not None:
        click.echo("  local_wvfrms: " + str(local_wvfrms))
    elif fdsn is not None:
        click.echo("  fdsn: " + str(fdsn))
    elif db_config is not None:
        db_info = cnfg.ConfigParser()
        db_info.read(db_config)
        click.echo("  db_config: " + str(db_config))
    else:
        click.echo("Invalid data parameters.  Analysis requires 1 of:")
        click.echo("  local_wvfrms")
        click.echo("  fdsn")
        click.echo("  db_url (and other database info)")
        
    freq_min = config.set_param(user_config, 'YIELD', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'YIELD', 'freq_max', freq_max, 'float')
    yld_min = config.set_param(user_config, 'YIELD', 'yld_min', yld_min, 'float')
    yld_max = config.set_param(user_config, 'YIELD', 'yld_max', yld_max, 'float')
    ref_rng = config.set_param(user_config, 'YIELD', 'ref_rng', ref_rng, 'float')
    resolution = config.set_param(user_config, 'YIELD', 'resolution', resolution, 'int')

    noise_option = config.set_param(user_config, 'YIELD', 'noise_option', noise_option, 'str')
    window_buffer = config.set_param(user_config, 'YIELD', 'window_buffer', window_buffer, 'float')
    amb_press = config.set_param(user_config, 'YIELD', 'amb_press', amb_press, 'float')
    amb_temp = config.set_param(user_config, 'YIELD', 'amb_temp', amb_temp, 'float')
    grnd_burst = config.set_param(user_config, 'YIELD', 'grnd_burst', grnd_burst, 'bool')
    exp_type = config.set_param(user_config, 'YIELD', 'exp_type', exp_type, 'str')

    click.echo('\n' + "Algorithm parameters:")
    click.echo("  freq_min: " + str(freq_min))
    click.echo("  freq_max: " + str(freq_max))
    click.echo("  yld_min: " + str(yld_min))
    click.echo("  yld_max: " + str(yld_max))
    click.echo("  ref_rng: " + str(ref_rng))
    click.echo("  resolution: " + str(resolution))
    click.echo("  noise_option: " + str(noise_option))
    click.echo("  window_buffer: " + str(window_buffer))
    click.echo("  amb_press: " + str(amb_press))
    click.echo("  amb_temp: " + str(amb_temp))
    click.echo("  grnd_burst: " + str(grnd_burst))
    click.echo("  exp_type: " + str(exp_type))


    det_list = data_io.json_to_detection_list(detect_label)
    if local_wvfrms is not None:
        stream, _ = data_io.set_stream(local_wvfrms, None, None)
    else:
        click.echo("Non-local waveform ingestion for yield estimation isn't set up yet...")
        return

    click.echo("Collecting waveform data for each detection...")
    st_list = [Stream([tr for tr in stream if det.station in tr.stats.station]) for det in det_list]
    for n, det in enumerate(det_list):
        click.echo('\n' + "Detection network.station: " + det.network + "." + det.station)
        print(st_list[n])
    click.echo('')

    smn_specs = spye.extract_spectra(det_list, st_list, win_buffer=window_buffer, ns_opt=noise_option)
    
    # ######################### #
    #     Load TLoss Models     #
    # ######################### #
    click.echo("Loading transmission loss statistics...")
    tlm_dir = os.path.dirname(tlm_label)
    tlm_pattern = tlm_pattern = tlm_label.split("/")[-1]
    tlm_files = [file_name for file_name in np.sort(os.listdir(tlm_dir)) if fnmatch.fnmatch(file_name, tlm_pattern + "*")]

    models = [0] * 2
    models[0] = [float(file_name.split("Hz")[0][len(tlm_pattern):]) for file_name in tlm_files]
    models[1] = [0] * len(tlm_files)
    for n in range(len(tlm_files)):
        models[1][n] = infrasound.TLossModel()
        models[1][n].load(tlm_dir + "/" + tlm_files[n])

    # ######################## #
    #         Run Yield        #
    #    Estimation Methods    #
    # ######################## #
    spye_result = spye.run(det_list, smn_specs, src_loc, np.array([freq_min, freq_max]), models, 
                            yld_rng=np.array([yld_min * 1.0e3, yld_max * 1.0e3]), ref_src_rng=ref_rng, 
                            resol=resolution, grnd_brst=grnd_burst, p_amb=amb_press, T_amb=amb_temp, exp_type=exp_type)

    if ".yld.json" not in local_yld_label:
        local_yld_label = local_yld_label + ".yld.json"
    click.echo("Writing yield estimate result into " + local_yld_label)
    data_io.write_json(spye_result, local_yld_label)

    click.echo('\n' + 'Results Summary (tons eq. TNT):')
    click.echo('\t' + "Maximum a Posteriori Yield: " + str(spye_result['yld_vals'][np.argmax(spye_result['yld_pdf'])]))
    click.echo('\t' + "68% Confidence Bounds: " + str(spye_result['conf_bnds'][0]))
    click.echo('\t' + "95% Confidence Bounds: " + str(spye_result['conf_bnds'][1]))
    click.echo('')


@click.command('single-station', short_help="Estimate the near-source spectral amplitude from a single station")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-wvfrms", help="Local waveform data files", default=None)
@click.option("--fdsn", help="FDSN source for waveform data files", default=None)
@click.option("--db-config", help="Database configuration file", default=None)
@click.option("--local-detect-label", help="Detection results path", default=None)
@click.option("--local-loc-label", help="Localization results path", default=None)
@click.option("--local-pdf-label", help="Output file for results", default=None)
@click.option("--tlm-label", help="Transmission loss model (TLM) path", default=None)
@click.option("--det-index", help="Index of detection in file", default=0, type=int)
@click.option("--src-lat", help="Source latitude (if no loc result file)", default=None)
@click.option("--src-lon", help="Source longitude (if no loc result file)", default=None)
@click.option("--freq-min", help="Minimum frequency (default: " + config.defaults['YIELD']['freq_min'] + " [Hz])", default=None, type=float)
@click.option("--freq-max", help="Maximum frequency (default: " + config.defaults['YIELD']['freq_max'] + " [Hz])", default=None, type=float)
@click.option("--ref-rng", help="Reference range for blastwave model (default " + config.defaults['YIELD']['ref_rng'] + " km)", default=1.0, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['YIELD']['resolution'] + ")", default=None, type=int)
@click.option("--noise-option", help="Noise option ('pre', 'post', or 'beam')", default=None)
@click.option("--window-buffer", help="Window buffer scaling (default: " + config.defaults['YIELD']['window_buffer'] + ")", default=None, type=float)
def single_station(config_file, local_wvfrms, fdsn, db_config, detect_label, local_loc_label, local_pdf_label, tlm_label, 
                    det_index, src_lat, src_lon, freq_min, freq_max, ref_rng, resolution, noise_option, window_buffer):
    '''
    Run Spectral Yield Estimation (SpYE) methods to estimate the near-source acoustic spectral amplitude for a single detecting station

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_spye single-station --local-wvfrms '../infrapy-data/hrr-5/W220/*.sac' --local-detect-label data/HRR-5.dets.json --det-index 0 --src-lat 33.5377 --src-lon -106.333961 --tlm-label "../infrapy/propagation/priors/tloss/2007_08-" --local-pdf-label "HRR-5_W220"
    \tinfrapy run_spye single-station --local-wvfrms '../infrapy-data/hrr-5/W240/*.sac' --local-detect-label data/HRR-5.dets.json --det-index 1 --src-lat 33.5377 --src-lon -106.333961 --tlm-label "../infrapy/propagation/priors/tloss/2007_08-" --local-pdf-label "HRR-5_W240"
    \tinfrapy run_spye single-station --local-wvfrms '../infrapy-data/hrr-5/W340/*.sac' --local-detect-label data/HRR-5.dets.json --det-index 2 --src-lat 33.5377 --src-lon -106.333961 --tlm-label "../infrapy/propagation/priors/tloss/2007_08-" --local-pdf-label "HRR-5_W340"
    \tinfrapy run_spye single-station --local-wvfrms '../infrapy-data/hrr-5/W420/*.sac' --local-detect-label data/HRR-5.dets.json --det-index 3 --src-lat 33.5377 --src-lon -106.333961 --tlm-label "../infrapy/propagation/priors/tloss/2007_08-" --local-pdf-label "HRR-5_W420"
    \tinfrapy run_spye single-station --local-wvfrms '../infrapy-data/hrr-5/W460/*.sac' --local-detect-label data/HRR-5.dets.json --det-index 4 --src-lat 33.5377 --src-lon -106.333961 --tlm-label "../infrapy/propagation/priors/tloss/2007_08-" --local-pdf-label "HRR-5_W460"
    '''

    click.echo("")
    click.echo("########################################")
    click.echo("##                                    ##")
    click.echo("##              InfraPy               ##")
    click.echo("##  Spectral Yield Estimation (SpYE)  ##")
    click.echo("##       Single Station Analysis      ##")
    click.echo("##                                    ##")
    click.echo("########################################")
    click.echo("")     

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None    

    # Waveform info
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')  
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None
    
    # Data IO parameters
    detect_label = config.set_param(user_config, 'DETECTION IO', 'detect_label', detect_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')
    src_lat = config.set_param(user_config, 'YIELD', 'src_lat', src_lat, 'float')
    src_lon = config.set_param(user_config, 'YIELD', 'src_lon', src_lon, 'float')

    click.echo('\n' + "Data parameters:")
    click.echo("  detect_label: " + str(detect_label))
    click.echo("  det_index: " + str(det_index))
    click.echo("  tlm_label: " + str(tlm_label))
    click.echo("  local_pdf_label: " + str(local_pdf_label))
    click.echo("  local_loc_label: " + str(local_loc_label))
    if local_loc_label is not None:
        src_loc = [local_loc_label['lat_mean'], local_loc_label['lon_mean']]
        # src_loc = [local_loc_label['lat_MaP'], local_loc_label['lon_MaP']]
        click.echo("    src_loc (from [...].loc.json file):", src_loc)
    elif src_lat is not None and src_lon is not None:
        src_loc = [src_lat, src_lon]
        click.echo("    src_lat: " + str(src_lat))
        click.echo("    src_lon: " + str(src_lon))
    else:
        print("Error.  SpYE analysis requires either localization analysis file or specified source location")
        return 0

    if local_wvfrms is not None:
        click.echo("  local_wvfrms: " + str(local_wvfrms))
    elif fdsn is not None:
        click.echo("  fdsn: " + str(fdsn))
    elif db_config is not None:
        db_info = cnfg.ConfigParser()
        db_info.read(db_config)
        click.echo("  db_config: " + str(db_config))
    else:
        click.echo("Invalid data parameters.  Analysis requires 1 of:")
        click.echo("  local_wvfrms")
        click.echo("  fdsn")
        click.echo("  db_url (and other database info)")
        
    freq_min = config.set_param(user_config, 'YIELD', 'freq_min', freq_min, 'float')
    freq_max = config.set_param(user_config, 'YIELD', 'freq_max', freq_max, 'float')
    ref_rng = config.set_param(user_config, 'YIELD', 'ref_rng', ref_rng, 'float')
    resolution = config.set_param(user_config, 'YIELD', 'resolution', resolution, 'int')
    noise_option = config.set_param(user_config, 'YIELD', 'noise_option', noise_option, 'str')
    window_buffer = config.set_param(user_config, 'YIELD', 'window_buffer', window_buffer, 'float')

    click.echo('\n' + "Algorithm parameters:")
    click.echo("  freq_min: " + str(freq_min))
    click.echo("  freq_max: " + str(freq_max))
    click.echo("  ref_rng: " + str(ref_rng))
    click.echo("  resolution: " + str(resolution))
    click.echo("  noise_option: " + str(noise_option))
    click.echo("  window_buffer: " + str(window_buffer))

    # Load detection and stream info
    det_list = data_io.json_to_detection_list(detect_label)
    if local_wvfrms is not None:
        stream, _ = data_io.set_stream(local_wvfrms, None, None)
    else:
        click.echo("Non-local waveform ingestion for yield estimation isn't set up yet...")
        return

    temp = spye.extract_spectra([det_list[det_index]], [stream], win_buffer=window_buffer, ns_opt=noise_option)[0]
    f = temp[0]
    spec_amp = temp[1]

    # ######################### #
    #     Load TLoss Models     #
    # ######################### #
    click.echo("Loading transmission loss statistics...")
    tlm_dir = os.path.dirname(tlm_label)
    tlm_pattern = tlm_pattern = tlm_label.split("/")[-1]
    tlm_files = [file_name for file_name in np.sort(os.listdir(tlm_dir)) if fnmatch.fnmatch(file_name, tlm_pattern + "*")]

    tlms = [0] * 2
    tlms[0] = [float(file_name.split("Hz")[0][len(tlm_pattern):]) for file_name in tlm_files]
    tlms[1] = [0] * len(tlm_files)
    for n in range(len(tlm_files)):
        tlms[1][n] = infrasound.TLossModel()
        tlms[1][n].load(tlm_dir + "/" + tlm_files[n])

    # Define grid and estimate near-source spectral amplitude
    click.echo("Computing near-source spectral amplitude PDF...")
    f_grid, spec_grid, pdf = spye._single_station(det_list[det_index], np.vstack((f, spec_amp)), [src_lat, src_lon], tlms, [freq_min, freq_max], resolution, ref_rng)

    click.echo("Saving PDF info to " + local_pdf_label + ".spye_pdf.npz")
    np.savez(local_pdf_label + ".spye_pdf.npz", f_grid, spec_grid, pdf)


@click.command('combine', short_help="Combine near-source spectral amplitude PDFs from stations")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-pdf-label", help="Output file for results", default=None)
@click.option("--local-yld-label", help="Output file for results", default=None)
@click.option("--yld-min", help="Minimum yield (default: " + config.defaults['YIELD']['yld_min'] + " [tons eq. TNT])", default=1.0, type=float)
@click.option("--yld-max", help="Maximum yield (default: " + config.defaults['YIELD']['yld_max'] + " [tons eq. TNT])", default=1000.0, type=float)
@click.option("--ref-rng", help="Reference range for blastwave model (default " + config.defaults['YIELD']['ref_rng'] + " km)", default=1.0, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['YIELD']['resolution'] + ")", default=None, type=int)
@click.option("--amb-press", help="Ambient pressure (default: " + config.defaults['YIELD']['amb_press'] + " [Pa])", default=None, type=float)
@click.option("--amb-temp", help="Ambient temperature (default: " + config.defaults['YIELD']['amb_temp'] + " [K])", default=None, type=float)
@click.option("--grnd-burst", help="Ground burst assumption (default: " + config.defaults['YIELD']['grnd_burst'] + " [Hz])", default=None, type=bool)
@click.option("--exp-type", help="Explosion type ('chemical' or 'nuclear')", default=None)

def combine(config_file, local_pdf_label, local_yld_label, yld_min, yld_max, ref_rng, resolution, amb_press, amb_temp, grnd_burst, exp_type):
    '''
    Run Spectral Yield Estimation (SpYE) methods to combine near-source acoustic spectral amplitude from single station analyses

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_spye combine --local-pdf-label 'HRR-5*.npz' --local-yld-label HRR_5-separate
    '''

    click.echo("")
    click.echo("########################################")
    click.echo("##                                    ##")
    click.echo("##              InfraPy               ##")
    click.echo("##  Spectral Yield Estimation (SpYE)  ##")
    click.echo("##   Combine Single Station Results   ##")
    click.echo("##                                    ##")
    click.echo("########################################")
    click.echo("")     

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo("Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None  


    click.echo('\n' + "Data parameters:")
    click.echo("  local_pdf_label: " + str(local_pdf_label))
    click.echo("  local_yld_label: " + str(local_yld_label))

    yld_min = config.set_param(user_config, 'YIELD', 'yld_min', yld_min, 'float')
    yld_max = config.set_param(user_config, 'YIELD', 'yld_max', yld_max, 'float')
    ref_rng = config.set_param(user_config, 'YIELD', 'ref_rng', ref_rng, 'float')
    resolution = config.set_param(user_config, 'YIELD', 'resolution', resolution, 'int')
    amb_press = config.set_param(user_config, 'YIELD', 'amb_press', amb_press, 'float')
    amb_temp = config.set_param(user_config, 'YIELD', 'amb_temp', amb_temp, 'float')
    grnd_burst = config.set_param(user_config, 'YIELD', 'grnd_burst', grnd_burst, 'bool')
    exp_type = config.set_param(user_config, 'YIELD', 'exp_type', exp_type, 'str')

    click.echo('\n' + "Algorithm parameters:")
    click.echo("  yld_min: " + str(yld_min))
    click.echo("  yld_max: " + str(yld_max))
    click.echo("  ref_rng: " + str(ref_rng))
    click.echo("  resolution: " + str(resolution))
    click.echo("  amb_press: " + str(amb_press))
    click.echo("  amb_temp: " + str(amb_temp))
    click.echo("  grnd_burst: " + str(grnd_burst))
    click.echo("  exp_type: " + str(exp_type))


    # Set file list for ingestion of all single-station results
    file_list = []
    if "/" in local_pdf_label:
        dir_files = os.listdir(os.path.dirname(local_pdf_label))
    else:
        dir_files = os.listdir(".")

    for file in np.sort(dir_files):
        if fnmatch.fnmatch(file, os.path.basename(local_pdf_label)):
            file_list += [file]

    # Load single-station PDFs
    click.echo('\n' + "Loading and interpolating near-source spectral estimates...")

    # combine and project onto blastwave model
    spye_result = spye._combine(file_list, np.array([yld_min * 1.0e3, yld_max * 1.0e3]), ref_rng, resolution, amb_press, amb_temp, grnd_burst, exp_type)

    if ".yld.json" not in local_yld_label:
        local_yld_label = local_yld_label + ".yld.json"
    click.echo("Writing yield estimate result into " + local_yld_label)
    data_io.write_json(spye_result, local_yld_label)

    click.echo('\n' + 'Results Summary (tons eq. TNT):')
    click.echo('\t' + "Maximum a Posteriori Yield: " + str(spye_result['yld_vals'][np.argmax(spye_result['yld_pdf'])]))
    click.echo('\t' + "68% Confidence Bounds: " + str(spye_result['conf_bnds'][0]))
    click.echo('\t' + "95% Confidence Bounds: " + str(spye_result['conf_bnds'][1]))
    click.echo('')
