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

for blob_id, blob, max_points in zip(range(len(blobs)), blobs, max_points_array):
	
        label_num = blob.blob_label
	if label_num != 10:
	    continue
        roi = pickle.load(open('roi'+str(label_num)+'.cpl', 'rb'))
	markers = np.zeros_like(roi)
	print 'Blob id:', blob_id, ', blob label:', blob.blob_label, ', blob size:', blob.blob_size

	for max_x, max_y, max_z, max_id in zip(max_points[0], max_points[1], max_points[2], range(len(max_points[0]))):
		print '\t', max_x, max_y, max_z, roi[max_x][max_y][max_z]
		markers[max_x][max_y][max_z] = max_id+1

        slice_x, slice_y, slice_z = blob.slice_x, blob.slice_y, blob.slice_z
#	local_max = peak_local_max(roi, min_distance=10, threshold_rel=0.1, exclude_border=True, indices=False)
#	markers = ndimage.label(local_max)[0]
	labels = watershed(-roi, markers, mask=(roi>0))
	print 'Markers', label_num, np.unique(markers)

	#distance = ndimage.distance_transform_edt(roi > 0)
#	write_image('roi_slice.png', roi[10, ...], vmin=0)
#	local_max = peak_local_max(distance, min_distance=10, threshold_rel=0.5, exclude_border=False, indices=False, labels=(roi > 0))
	local_max = np.zeros_like(roi)
	for x, y, z, i in zip(max_points_array[3][0], max_points_array[3][1], max_points_array[3][2], range(len(max_points_array[3][0]))):
        #	print x, y, z
        	local_max[x][y][z] = i+1
	print 'Local max num:', np.count_nonzero(local_max)

#markers = ndimage.label(local_max)[0]
	#labels = watershed(-distance, local_max, mask=(roi > 0))
#write_image('markers_slice.png', markers[10, ...], vmin=0)
#write_image('labels_slice.png', labels[10, ...], vmin=0)
#write_image('distance_slice.png', distance[10, ...], vmin=0)

	#print 'Markers', np.unique(markers)

	directory = 'stack_' + str(label_num)
	if not os.path.exists(directory):
    	    os.makedirs(directory)
	    
	for i in range(np.shape(labels)[0]):
	    	write_image(directory + '/roi_' + str(i) + '.png', roi[i, :, :], vmin=0, vmax=np.amax(roi))
		write_image(directory + '/layer_' + str(i) + '_watershed.png', labels[i, :, :], vmin=0)
