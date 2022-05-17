#!/usr/bin/env python -W ignore::DeprecationWarning

# test_assoc.py
#
# Author    Philip Blom (pblom@lanl.gov)

import numpy as np

import pathos.multiprocessing as mp
from multiprocessing import cpu_count

from infrapy.association import hjl
from infrapy.utils import data_io

from multiprocess import Pool


if __name__ == '__main__':
    #########################
    ### Define parameters ###
    #########################

    # Read in detections from file
    det_list = data_io.json_to_detection_list('data/example1.dets.json')

    # define joint-likelihood calculation parameters
    width = 10.0
    rng_max = 3000.0

    # define clustering parameters
    dist_max = 10.0
    clustering_threshold = 5.0
    trimming_thresh = 3.0
    
    pl = Pool(cpu_count() - 1)
    ######################
    #### Run analysis ####
    ######################
    labels, dists = hjl.run(det_list, clustering_threshold, dist_max=dist_max, bm_width=width, rng_max=rng_max, trimming_thresh=trimming_thresh, pool=pl,show_result=True)
    
    # Summarize clusters
    clusters, qualities = hjl.summarize_clusters(labels, dists)
    for n in range(len(clusters)):
        print("Cluster:", clusters[n], '\t', "Cluster Quality:", 10.0**(-qualities[n]))
    
    pl.close()
    pl.terminate()






