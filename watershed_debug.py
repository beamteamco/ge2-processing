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

blobs = pickle.load(open('ge_blobs.cpl', 'rb'))
max_points_array = pickle.load(open('ge_spot_max_points.cpl', 'rb'))
roi = pickle.load(open('roi10.cpl', 'rb'))
#roi = roi[10, ...]
distance = ndimage.distance_transform_edt(roi > 0)
write_image('roi_slice.png', roi[10, ...], vmin=0)
#local_max = peak_local_max(distance, min_distance=10, threshold_rel=0.5, exclude_border=False, indices=False, labels=(roi > 0))
local_max = np.zeros_like(roi)
for x, y, z, i in zip(max_points_array[3][0], max_points_array[3][1], max_points_array[3][2], range(len(max_points_array[3][0]))):
	print x, y, z
	local_max[x][y][z] = i+1
print 'Local max num:', np.count_nonzero(local_max)

#markers = ndimage.label(local_max)[0]
labels = watershed(-distance, local_max, mask=(roi > 0))
#write_image('markers_slice.png', markers[10, ...], vmin=0)
write_image('labels_slice.png', labels[10, ...], vmin=0)
write_image('distance_slice.png', distance[10, ...], vmin=0)

print 'Markers', np.unique(markers)
