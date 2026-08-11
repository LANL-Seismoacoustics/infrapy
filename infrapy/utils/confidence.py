# confidence.py
#
# Method to compute confidence intervals for 1D function f(x)

import numpy as np

from scipy.integrate import simpson


def find_confidence(func, lims, conf_aim, resol=1e3):
    if conf_aim > 1.0:
        print("WARNING - find_confidence cannot use conf > 1.0")
        return lims, 1.0, 0.0

    resol = int(resol)
    x_vals = np.linspace(lims[0], lims[1], resol)
    f_vals = func(x_vals)

    norm = simpson(f_vals, x_vals)

    f_max = np.max(f_vals)
    thresh_vals = np.linspace(0.0, f_max, resol)

    for n in range(resol):
        temp = f_vals.copy()
        temp[temp < thresh_vals[n]] = 0.0

        conf = simpson(temp, x_vals) / norm

        if conf < conf_aim:
            thresh = thresh_vals[n]
            bnds = [np.min(x_vals[f_vals > thresh]), np.max(x_vals[f_vals > thresh])]
            break
        
    return bnds, conf, thresh
