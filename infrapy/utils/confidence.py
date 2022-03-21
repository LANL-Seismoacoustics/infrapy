# confidence.py
#
# Method to compute confidence intervals for 1D function f(x)

import numpy as np

from scipy.integrate import quad


def find_confidence(func, lims, conf_aim):
    if conf_aim > 1.0:
        print("WARNING - find_confidence cannot use conf > 1.0")
        return lims

    def conf_func(x, thresh):
        val = func(x)
        if val >= thresh:
            return val
        else:
            return 0.0

    resol = int(1e3)
    x_vals = np.linspace(lims[0], lims[1], resol)
    f_vals = func(x_vals)

    f_max = max(f_vals)
    thresh_vals = np.linspace(0.0, f_max, resol)
    
    norm = quad(func, lims[0], lims[1], limit=100, epsrel=1.0e-3)[0]
    
    conf_prev=1.0
    bnds=[]
    for n in range(resol):
        conf = quad(conf_func, lims[0], lims[1], (thresh_vals[n],), limit=100, epsrel=1.0e-3)[0] / norm
        
        if conf < conf_aim < conf_prev:
            thresh = thresh_vals[n - 1] - (thresh_vals[n-1] - thresh_vals[n]) / (conf_prev - conf) * (conf_aim - conf)
            conf = quad(conf_func, lims[0], lims[1], (thresh,), limit=100, epsrel=1.0e-3)[0] / norm
            
            for n in range(resol - 1):
                if (f_vals[n] - thresh) * (f_vals[n+1] - thresh) < 0.0:
                    bnds.append(x_vals[n])
            break
        
        conf_prev=conf

    return bnds, conf, thresh
