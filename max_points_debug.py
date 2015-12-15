import copy
import logging
import os
import sys
import time
import warnings

import yaml
import cPickle as pickle
import pprint

from hexrd import config
from ge_processor.ge_pre_processor import *

import numpy as np
from scipy import ndimage as ndimage

from skimage.feature import peak_local_max
from skimage.morphology import watershed

from joblib import Parallel, delayed  
import multiprocessing

def get_local_maxima(blob, roi):
	markers = np.zeros_like(roi)
	print 'Blob label:', blob.blob_label, ', blob size:', blob.blob_size

	max_points = peak_local_max(roi, min_distance=min_peak_separation, threshold_rel=0.2, 
				    exclude_border=False, indices=False)
	max_points[roi < 20] = 0
	max_points = np.nonzero(max_points)

	for max_x, max_y, max_z, max_id in zip(max_points[0], max_points[1], max_points[2], 
					       range(len(max_points[0]))):
		print '\t', max_x, max_y, max_z, roi[max_x][max_y][max_z]
		markers[max_x][max_y][max_z] = max_id+1

	return max_points
#--- END

min_peak_separation = 5
#
blobs = pickle.load(open('ge_blobs.cpl', 'rb'))
rois = []

for blob in blobs:
	roi = pickle.load(open('roi' + str(blob.blob_label) + '.cpl', 'rb'))
	rois.append(roi)

num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=num_cores, verbose=5)(delayed(get_local_maxima)(blob, roi) for blob, roi in zip(blobs, rois)) 
