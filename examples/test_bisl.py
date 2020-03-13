#!/usr/bin/env python

# test_bisl.py
#
# Author    Philip Blom (pblom@lanl.gov)

import numpy as np

from infrapy.location import bisl
from infrapy.propagation import likelihoods as lklhds
from infrapy.propagation import infrasound as infsnd
from infrapy.location import plot_loc_res


# ######################### #
#       Define Inputs       #
# ######################### #

# Define ground_truth if known (41.131, -112.896 for UTTR; Test includes show in June 2004)
grnd_trth = [41.131, -112.896, np.datetime64('2004-06-02T17:23:04.0')]

# Define localization parameters
bm_width = 12.5
rad_min, rad_max = 50.0, 500.0
rng_max = np.pi / 2.0 * 6370.0
resolution = int(np.sqrt(1e5))

'''
# Define the list of detections (output from association)
# detection format: (lat, lon, arrival time, back az, F stat, elements)
# arrival time format: datetime.datetime(year, month, day, hour, minute, second)
det1 = lklhds.InfrasoundDetection(42.7668, -109.5939, np.datetime64('2004-06-02T17:42:14.0'), -125.6, 75.0, 4)
det2 = lklhds.InfrasoundDetection(38.4296, -118.3036, np.datetime64('2004-06-02T17:50:38.0'),   56.6, 75.0, 4)
det3 = lklhds.InfrasoundDetection(48.2641, -117.1257, np.datetime64('2004-06-02T18:09:14.0'),  157.5, 75.0, 4)
det_list = [det1, det2, det3]
'''

# Load detection list from flat file
#det_list = lklhds.file2dets("data/detection_set2.dat")

# Load detection list from json file
det_list = lklhds.json_to_detection_list('data/detection_set2.json')


# Plot detections
plot_loc_res.plot_dets(det_list, display=True,save_fig=None)

# ########################## #
#          Run BISL          #
#       in Verbose Mode      #
# ########################## #

# Run analysis without priors
result,pdf = bisl.run(det_list,
                    bm_width=bm_width,
                    rad_min=rad_min, 
                    rad_max=rad_max, 
                    rng_max=rng_max, 
                    resol=resolution,angle=[-180,180])

summary = bisl.summarize(result)

print('-' * 75)
print('BISL Summary\n')
print(summary)
print('\n' + '-'*75 + '\n')

# ###############################
# plot the results from bisl ####

plot_loc_res.plot_loc(det_list, result,pdf,confidence=[95,99],display=True,save_fig=None,grnd_trth =grnd_trth)
plot_loc_res.plot_time(result,display=True, save_fig=None,grnd_trth =grnd_trth)
# Define priors, load from file, and display
model = infsnd.PathGeometryModel()
model.load("../infrapy/propagation/priors/UTTR_models/UTTR_06_1800UTC.pgm")
#model.display()

# Re-run analysis with priors
result,pdf = bisl.run(det_list, 
                    bm_width=bm_width,
                    rad_min=rad_min, 
                    rad_max=rad_max, 
                    rng_max=rng_max, 
                    resol=resolution,angle=[-180,180],
                    path_geo_model=model)


summary = bisl.summarize(result)

print('-' * 75)
print('BISL Summary\n')
print(summary)
print('\n' + '-'*75 + '\n')

plot_loc_res.plot_loc(det_list, result,pdf,confidence=[95,99],display=True,save_fig=None,grnd_trth =grnd_trth)
plot_loc_res.plot_time(result,display=True, save_fig=None,grnd_trth =grnd_trth)
