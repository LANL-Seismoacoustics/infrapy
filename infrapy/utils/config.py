#!/usr/bin/env python

import os
import configparser as cnfg
import numpy as np

from pathlib import Path

from ..propagation import infrasound

from importlib.util import find_spec
if find_spec('stochprop'):
    import stochprop.propagation as sp_prop

# Set up default configuation
defaults = cnfg.ConfigParser()
defaults.read(str(Path(__file__).parent.parent / "resources" / "default.config"))

def set_param(user_config, section, param, cli_val, format='float'):
    if cli_val is not None:
        # return CLI value if entered
        return cli_val
    elif user_config is not None:
        # check if parameter is in user config and return default if it's not
        try:
            if user_config[section][param] == "None":
                return None
            else:
                if format == 'float':
                    return float(user_config[section][param])
                elif format == 'int':
                    return int(user_config[section][param])
                elif format == 'bool':
                    return user_config[section].getboolean(param)
                else:
                    return user_config[section][param].strip('"\'')
        except:
            try:
                if defaults[section][param] == "None":
                    return None
                else:
                    if format == 'float':
                        return float(defaults[section][param])
                    elif format == 'int':
                        return int(defaults[section][param])
                    elif format == 'bool':
                        return defaults[section].getboolean(param)
                    else:
                        return defaults[section][param].strip('"\'')
            except:
                return None
    else:
        # use default values if no CLI and no user config
        try:
            if defaults[section][param] == "None":
                return None
            else:
                if format == 'float':
                    return float(defaults[section][param])
                elif format == 'int':
                    return int(defaults[section][param])
                elif format == 'bool':
                    return defaults[section].getboolean(param)
                else:
                    return defaults[section][param].strip('"\'')
        except:
            return None


###################################
##  Configure Without CLI Input  ##
###################################
def set_bm_config(cnfg_file):
    bm_cnfg = cnfg.ConfigParser()
    if os.path.isfile(cnfg_file):
        bm_cnfg.read(cnfg_file)
    else:
        print('\nInvalid configuration file (' + cnfg_file + ' not found).\nDefault parameters will be used.\n')

    bm_params = {}
    bm_params['freq_min'] = set_param(bm_cnfg, 'FK', 'freq_min', None, 'float')
    bm_params['freq_max'] = set_param(bm_cnfg, 'FK', 'freq_max', None, 'float')
    bm_params['back_az_min'] = set_param(bm_cnfg, 'FK', 'back_az_min', None, 'float')
    bm_params['back_az_max'] = set_param(bm_cnfg, 'FK', 'back_az_max', None, 'float')
    bm_params['back_az_step'] = set_param(bm_cnfg, 'FK', 'back_az_step', None, 'float')
    bm_params['trace_vel_min'] = set_param(bm_cnfg, 'FK', 'trace_vel_min', None, 'float')
    bm_params['trace_vel_max'] = set_param(bm_cnfg, 'FK', 'trace_vel_max', None, 'float')
    bm_params['trace_vel_step'] = set_param(bm_cnfg, 'FK', 'trace_vel_step', None, 'float')
    bm_params['method'] = set_param(bm_cnfg, 'FK', 'method', None, 'string')
    bm_params['signal_start'] = set_param(bm_cnfg, 'FK', 'signal_start', None, 'string')
    bm_params['signal_end'] = set_param(bm_cnfg, 'FK', 'signal_end', None, 'string')
    bm_params['noise_start'] = set_param(bm_cnfg, 'FK', 'noise_start', None, 'string')
    bm_params['noise_end'] = set_param(bm_cnfg, 'FK', 'noise_end', None, 'string')
    bm_params['window_len'] = set_param(bm_cnfg, 'FK', 'window_len', None, 'float')
    bm_params['sub_window_len'] = set_param(bm_cnfg, 'FK', 'sub_window_len', None, 'float')
    bm_params['window_step'] = set_param(bm_cnfg, 'FK', 'window_step', None, 'float')

    det_params = {}
    det_params['window_len'] = set_param(bm_cnfg, 'FD', 'window_len', None, 'float')
    det_params['p_value'] = set_param(bm_cnfg, 'FD', 'p_value', None, 'float')
    det_params['min_duration'] = set_param(bm_cnfg, 'FD', 'min_duration', None, 'float')
    det_params['back_az_width'] = set_param(bm_cnfg, 'FD', 'back_az_width', None, 'float')
    det_params['fixed_thresh'] = set_param(bm_cnfg, 'FD', 'fixed_thresh', None, 'float')
    det_params['thresh_ceil'] = set_param(bm_cnfg, 'FD', 'thresh_ceil', None, 'float')
    det_params['merge_dets'] = set_param(bm_cnfg, 'FD', 'merge_dets', None, 'bool')

    return bm_params, det_params


def set_ev_bld_config(cnfg_file):
    hjl_cnfg = cnfg.ConfigParser()
    if os.path.isfile(cnfg_file):
        hjl_cnfg.read(cnfg_file)
    else:
        print('\nInvalid configuration file (' + cnfg_file + ' not found).\nDefault parameters will be used.\n')

    hjl_params = {}
    hjl_params['celerity_model'] = set_param(hjl_cnfg, 'ASSOC', 'celerity_model', None, 'string')
    hjl_params['back_az_width'] = set_param(hjl_cnfg, 'ASSOC', 'back_az_width', None, 'float')
    hjl_params['range_max'] = set_param(hjl_cnfg, 'ASSOC', 'range_max', None, 'float')
    hjl_params['resolution'] = set_param(hjl_cnfg, 'ASSOC', 'resolution', None, 'int')
    hjl_params['distance_matrix_max'] = set_param(hjl_cnfg, 'ASSOC', 'distance_matrix_max', None, 'float')
    hjl_params['cluster_linkage'] = set_param(hjl_cnfg, 'ASSOC', 'cluster_linkage', None, 'string')
    hjl_params['cluster_threshold'] = set_param(hjl_cnfg, 'ASSOC', 'cluster_threshold', None, 'float')
    hjl_params['trimming_threshold'] = set_param(hjl_cnfg, 'ASSOC', 'trimming_threshold', None, 'float')
    hjl_params['ev_population_min'] = set_param(hjl_cnfg, 'ASSOC', 'ev_population_min', None, 'int')
    hjl_params['ev_station_min'] = set_param(hjl_cnfg, 'ASSOC', 'ev_station_min', None, 'int')

    return hjl_params


def set_ev_loc_config(cnfg_file):
    loc_cnfg = cnfg.ConfigParser()
    if os.path.isfile(cnfg_file):
        loc_cnfg.read(cnfg_file)
    else:
        print('\nInvalid configuration file (' + cnfg_file + ' not found).\nDefault parameters will be used.\n')

    loc_params = {}
    loc_params['back_az_width'] = set_param(loc_cnfg, 'LOC', 'back_az_width', None, 'float')
    loc_params['range_max'] = set_param(loc_cnfg, 'LOC', 'range_max', None, 'float')
    loc_params['grid_resol'] = set_param(loc_cnfg, 'LOC', 'grid_resol', None, 'int')
    loc_params['ll_corner'] = set_param(loc_cnfg, 'LOC', 'll_corner', None, 'str')
    loc_params['ur_corner'] = set_param(loc_cnfg, 'LOC', 'ur_corner', None, 'str')
    loc_params['latlon_resol'] = set_param(loc_cnfg, 'LOC', 'latlon_resol', None, 'float')
    loc_params['tm_min'] = set_param(loc_cnfg, 'LOC', 'tm_min', None, 'str')
    loc_params['tm_max'] = set_param(loc_cnfg, 'LOC', 'tm_max', None, 'str')
    loc_params['tm_resol'] = set_param(loc_cnfg, 'LOC', 'tm_resol', None, 'float')
    loc_params['celerity_model'] = set_param(loc_cnfg, 'LOC', 'celerity_model', None, 'str').lower()
    loc_params['pgm_file'] = set_param(loc_cnfg, 'LOC', 'pgm_file', None, 'str')
    loc_params['atmo_data'] = set_param(loc_cnfg, 'LOC', 'atmo_data', None, 'str')
    loc_params['alt_lims'] = set_param(loc_cnfg, 'LOC', 'alt_lims', None, 'str')
    loc_params['alt_resol'] = set_param(loc_cnfg, 'LOC', 'alt_resol', None, 'float')
    loc_params['grnd_snd_spd'] = set_param(loc_cnfg, 'LOC', 'grnd_snd_spd', None, 'float')
    loc_params['c0_stdev'] = set_param(loc_cnfg, 'LOC', 'c0_stdev', None, 'float')
    loc_params['det_tm_stdev'] = set_param(loc_cnfg, 'LOC', 'det_tm_stdev', None, 'float')
    loc_params['az_limit'] = set_param(loc_cnfg, 'LOC', 'az_limit', None, 'float')
    loc_params['local_temp_dir'] = set_param(loc_cnfg, 'LOC', 'local_temp_dir', None, 'str')

    if loc_params['ll_corner'] is not None:
        loc_params['ll_corner'] = np.array([float(val) for val in loc_params['ll_corner'].replace(" ","").split(",")])
        loc_params['ur_corner'] = np.array([float(val) for val in loc_params['ur_corner'].replace(" ","").split(",")])
        tm_lims = (np.datetime64(loc_params['tm_min']), np.datetime64(loc_params['tm_max']))

        if loc_params['alt_lims'] is not None:
            loc_params['alt_lims'] = np.array([float(val) for val in loc_params['alt_lims'].replace(" ","").split(",")])

        loc_params['back_az_width'] = None
        loc_params['range_max'] = None
        loc_params['grid_resol'] = None
    else:
        loc_params['alt_resol'] = None
        loc_params['c0_stdev'] = None
        loc_params['az_limit'] = None
        tm_lims = None

    if loc_params['atmo_data'] is None:
        if loc_params['pgm_file'] is not None:
            if find_spec('stochprop'):
                print.echo("  pgm_file: " + str(loc_params['pgm_file']))
                pgm = sp_prop.PathGeometryModel()
                pgm.load(loc_params['pgm_file'])
            else:
                print('\n' + "Can't use a PGM without stochprop installed." + '\n' + "Built-in celerity model options are: 'regional_hf', 'regional_lf', and 'infGEM'" + '\n')
                return
        else:
            if "auto" in loc_params['celerity_model']:
                # pre-load the infGEM model and update once integration region is identified
                infrasound._load_celerity_model('infgem')
            else:
                infrasound._load_celerity_model(loc_params["celerity_model"])
            pgm = None
    else:
        loc_params["celerity_model"] = None
        pgm = None

    return loc_params, tm_lims, pgm


def set_ev_char_config(cnfg_file):
    char_cnfg = cnfg.ConfigParser()
    if os.path.isfile(cnfg_file):
        char_cnfg.read(cnfg_file)
    else:
        print('\nInvalid configuration file (' + cnfg_file + ' not found).\nDefault parameters will be used.\n')

    char_params = {}
    char_params['det_mask'] = set_param(char_cnfg, 'YIELD', 'det_mask', None, 'string')
    char_params['loc_index'] = set_param(char_cnfg, 'YIELD', 'loc_index', None, 'str')
    char_params['tlm_label'] = set_param(char_cnfg, 'YIELD', 'tlm_label', None, 'str')
    char_params['freq_min'] = set_param(char_cnfg, 'YIELD', 'freq_min', None, 'float')
    char_params['freq_max'] = set_param(char_cnfg, 'YIELD', 'freq_max', None, 'float')
    char_params['yld_min'] = set_param(char_cnfg, 'YIELD', 'yld_min', None, 'float')
    char_params['yld_max'] = set_param(char_cnfg, 'YIELD', 'yld_max', None, 'float')
    char_params['resolution'] = set_param(char_cnfg, 'YIELD', 'resolution', None, 'int')
    char_params['ref_rng'] = set_param(char_cnfg, 'YIELD', 'ref_rng', None, 'float')
    char_params['amb_press'] = set_param(char_cnfg, 'YIELD', 'amb_press', None, 'float')
    char_params['amb_temp'] = set_param(char_cnfg, 'YIELD', 'amb_temp', None, 'float')
    char_params['grnd_burst'] = set_param(char_cnfg, 'YIELD', 'grnd_burst', None, 'bool')
    char_params['exp_type'] = set_param(char_cnfg, 'YIELD', 'exp_type', None, 'str')

    return char_params
