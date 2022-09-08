#!/usr/bin/env python

from email.policy import default
import os 
import sys
import fnmatch
from threading import local 
import click
import configparser as cnfg
import numpy as np

from obspy import Stream

from ..characterization import spye
from ..location import bisl
from ..propagation import infrasound

from ..utils import config
from ..utils import data_io

@click.command('run_loc', short_help="Estimate source locations and times for events")
@click.option("--config-file", help="Configuration file", default=None)
@click.option("--local-detect-label", help="Detection path and pattern", default=None)
@click.option("--local-loc-label", help="Localization results path", default=None)
@click.option("--back-az-width", help="Width of beam projection (default: " + config.defaults['LOC']['back_az_width'] + " [deg])", default=None, type=float)
@click.option("--range-max", help="Maximum source-receiver range (default: " + config.defaults['LOC']['range_max'] + " [km])", default=None, type=float)
@click.option("--resolution", help="Number of points/dimension for numerical sampling (default: " + config.defaults['LOC']['resolution'] + ")", default=None, type=int)
@click.option("--src-est", help="Estimated source location and radius of region to consider (default: None)", default=None)
@click.option("--pgm-file", help="Path geometry model (PGM) file (default: None)", default=None)
def run_loc(config_file, local_detect_label, local_loc_label, back_az_width, range_max, resolution, src_est, pgm_file):
    '''
    Run Bayesian Infrasonic Source Localization (BISL) methods to estimate the source location and origin time for an event

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_loc --local-detect-label GJI_example-ev0  --local-loc-label GJI_example-ev0
    \tinfrapy run_loc --local-detect-label data/detection_set2.json --local-loc-label data/location2 --pgm-file ../infrapy/propagation/priors/UTTR_models/UTTR_06_1800UTC.pgm
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
    local_detect_label = config.set_param(user_config, 'DETECTION IO', 'local_detect_label', local_detect_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')

    if ".loc.json" in local_loc_label:
       local_loc_label = local_loc_label[:-9] 
            
    click.echo('\n' + "Data summary:")
    click.echo("  local_detect_label: " + str(local_detect_label))
    click.echo("  local_loc_label: " + str(local_loc_label))

    # Algorithm parameters
    back_az_width = config.set_param(user_config, 'LOC', 'back_az_width', back_az_width, 'float')
    range_max = config.set_param(user_config, 'LOC', 'range_max', range_max, 'float')
    resolution = config.set_param(user_config, 'LOC', 'resolution', resolution, 'int')
    src_est = config.set_param(user_config, 'LOC', 'src_est', src_est, 'string')
    pgm_file = config.set_param(user_config, 'LOC', 'pgm_file', pgm_file, 'str')

    if src_est is not None:
        src_est = [float(x.strip('[( )]')) for x in src_est.split(',')]

    click.echo('\n' + "Parameter summary:")
    click.echo("  back_az_width: " + str(back_az_width))
    click.echo("  range_max: " + str(range_max))
    click.echo("  resolution: " + str(resolution))
    click.echo("  src_est: " + str(src_est))
    click.echo("  pgm_file: " + str(pgm_file))

    if pgm_file is not None:
        click.echo("")
        pgm = infrasound.PathGeometryModel()
        pgm.load(pgm_file)
    else:
        pgm = None

    click.echo("")
    events = data_io.set_det_list(local_detect_label, merge=False)
    if type(events[0]) is list:
        # run localization analysis for multiple detection sets
        for j, det_list in enumerate(events):
            click.echo('\n' + "Running BISL on event " + str(j + 1) + " of " + str(len(events)))
            result = bisl.run(det_list, path_geo_model=pgm, custom_region=src_est, resol=resolution, bm_width=back_az_width, rng_max=range_max, rad_min=100.0, rad_max=range_max/4.0)

            # Determine output format for BISL results
            click.echo('\n' + "BISL Summary:")
            click.echo(bisl.summarize(result))

            click.echo("Writing localization result into " + local_loc_label + "_ev-" + str(j + 1) + "loc.json")
            data_io.write_json(result, local_loc_label + "_ev-" + str(j + 1) + "loc.json")

    else:
        # run a single localization analysis
        click.echo("")
        result = bisl.run(events, path_geo_model=pgm, custom_region=src_est, resol=resolution, bm_width=back_az_width, rng_max=range_max, rad_min=100.0, rad_max=range_max/4.0)

        # Determine output format for BISL results
        click.echo('\n' + "BISL Summary:")
        click.echo(bisl.summarize(result))


        if ".loc.json" not in local_loc_label:
            local_loc_label = local_loc_label + ".loc.json"
        click.echo("Writing localization result into " + local_loc_label)
        data_io.write_json(result, local_loc_label)


@click.command('run_yield', short_help="Estimate explosive yield from a surface explosion")
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

def run_yield(config_file, local_wvfrms, fdsn, db_config, local_detect_label, local_loc_label, local_yld_label, tlm_label, 
                src_lat, src_lon, freq_min, freq_max, yld_min, yld_max, ref_rng, resolution, noise_option, window_buffer,
                amb_press, amb_temp, grnd_burst, exp_type):
    '''
    Run Spectral Yield Estimation (SpYE) methods to estimate the equivalent TNT yield of an above-ground explosion

    \b
    Example usage (run from infrapy/examples directory):
    \tinfrapy run_yield --local-wvfrms '../infrapy-data/hrr-5/*/*.sac' --local-detect-label data/HRR-5.dets.json --src-lat 33.5377 --src-lon -106.333961 --tlm-label "../infrapy/propagation/priors/tloss/2007_08-" --local-yld-label "HRR-5"
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##             InfraPy             ##")
    click.echo("##    Yield Estimation Analysis    ##")
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

    # Waveform info
    local_wvfrms = config.set_param(user_config, 'WAVEFORM IO', 'local_wvfrms', local_wvfrms, 'string')
    fdsn = config.set_param(user_config, 'WAVEFORM IO', 'fdsn', fdsn, 'string')  
    db_config = config.set_param(user_config, 'WAVEFORM IO', 'db_config', db_config, 'string')
    db_info = None
    
    # Data IO parameters
    local_detect_label = config.set_param(user_config, 'DETECTION IO', 'local_detect_label', local_detect_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')
    local_loc_label = config.set_param(user_config, 'DETECTION IO', 'local_loc_label', local_loc_label, 'string')
    src_lat = config.set_param(user_config, 'YIELD', 'src_lat', src_lat, 'float')
    src_lon = config.set_param(user_config, 'YIELD', 'src_lon', src_lon, 'float')

    click.echo('\n' + "Data parameters:")
    click.echo("  local_detect_label: " + str(local_detect_label))
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


    det_list = data_io.json_to_detection_list(local_detect_label)
    if local_wvfrms is not None:
        stream, _ = data_io.set_stream(local_wvfrms, None, None)
    else:
        click.echo("Non-local waveform ingestion for yield estimation isn't set up yet...")
        return

    st_list = [Stream([tr for tr in stream if det.station in tr.stats.station]) for det in det_list]
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

