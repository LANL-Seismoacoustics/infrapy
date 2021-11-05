#!/usr/bin/env python

import imp 

import configparser as cnfg

from click.utils import echo

# Set up default configuation
defaults = cnfg.ConfigParser()
defaults.read(imp.find_module('infrapy')[1] + '/resources/default.config')

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
                    return user_config[section][param]
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
                        return defaults[section][param]
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
                    return defaults[section][param]
        except:
            return None

