import copy
import logging
import os
import sys
import time
import warnings
# Pickle saves a data structure to a binary file
try:
   import cPickle as pickle
except:
   import pickle
import yaml
# Numpy and Scipy for the good stuff (array operations, image processing)
import numpy as np
import scipy.ndimage as ndimage
from scipy.ndimage.filters import gaussian_filter
# Matplotlib plots like Matlab
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    from matplotlib import pyplot as plt
# hexrd helpers to read a config file and load GE2 data
from hexrd import config
from hexrd.coreutil import initialize_experiment
# Image analysis that is out of league for scipy
from skimage.morphology import watershed
from skimage.feature import peak_local_max
# Clustering stuff
from sklearn.cluster import DBSCAN
# Parallelization for spped
from joblib import Parallel, delayed
import multiprocessing
from random import randint

# Helper to save a 2D array as an image (thanks Branden Kappes @ mines.edu)
def write_image(filename, arr, pts=None, minsize=None, **kwds):
    '''
    Write a 2D array to a PNG image file. Optionally superimpose
    points that are described in terms of their x, y coordinates.
    '''
    xsize, ysize = 1024, 1024

    fig = plt.figure(figsize=(20, 20), dpi=120)
    # plot the image
    ax = fig.gca()
    kwds['interpolation'] = kwds.get('interpolation', 'none')
    image = ax.imshow(arr, **kwds)
    ax.set_xlabel(r'X', fontsize='large')
    ax.set_ylabel(r'Y', fontsize='large')
    cbar = plt.colorbar(image)
    # plot any points
    if pts is not None:
        pts = np.asarray(pts)
        ax.plot(pts[:,1], pts[:,0], 'go', markersize=6)
        # resize (since adding points often adds padding)
        ax.set_xlim(0, 2048)
        ax.set_ylim(0, 2048)
    fig.savefig(filename, bbox_inches='tight', pad_inches=1./3.)
    fig.clf()
    plt.close()
#--
def write_ge2(filename, arr, nbytes_header=8192, pixel_type=np.uint16):
    '''
    Write a 3D array to a GE2 file.
    '''
    fid = open(filename, 'wb')
    fid.seek(nbytes_header)
    fid.write(arr.astype(pixel_type))
    fid.close()
#--
def find_blobs_mp(ge_data, int_scale_factor, min_size, min_peak_separation, cfg):
   '''
   Multiprocessing worker subroutine for blob and local maxima finding.
   First finds blobs using a connected component type algorithm. Then
   finds the local maxima (spots). Then calculates the area corresponding to 
   each of the spots using the watershed algorithm.
   '''
   # Once again, remove noise
   ge_mask = ge_data > (int_scale_factor * cfg.fit_grains.threshold)
   # Dilate the mask to merge any small spots
   ge_mask = ndimage.binary_dilation(ge_mask, ndimage.generate_binary_structure(3, 3))
   # Label each connected component (spot)
   label_ge, number_of_labels = ndimage.label(ge_mask)
   ge_labeled = np.amax(label_ge, 1)
   label_rand = str(randint(0, 1000000))
   # Get size of each blob by counting number of pixels with the same label
   blob_sizes = ndimage.sum(ge_mask, label_ge, range(1, number_of_labels + 1))
   blob_labels = np.unique(label_ge)
   # Loop over all detected regions and filter based on size/aspect ratio
   # Save filtered blob information in an array 'blobs' of GEBlob object
   blobs = []
   blob_centroids = []
   max_points_global = []
   for label_num, blob_size in zip(blob_labels, blob_sizes):
        # Label 0 is for the whole image background. Do not want.
        if label_num == 0:
            continue
        # Total pixels in blob < min_size? Move on
        if blob_size < min_size:
            continue
        # Get the minimal region of interest (ROI) for the blob
        slice_x, slice_y, slice_z = ndimage.find_objects(label_ge==label_num)[0]
        bbox = [slice_x.stop -  slice_x.start, slice_y.stop -  slice_y.start, slice_z.stop -  slice_z.start]
        # Is any of roi bounding box dim < min? Move on
        if(min(bbox) < (min_size ** (1.0/3.0))):
            #print 'Finished one spot but its small'
            continue
        #print 'Yay finished one spot'
        #print 'Finished processing a blob with bbox', bbox
        blob_centroids.append([(slice_x.stop + slice_x.start)/2.0, (slice_y.stop + slice_y.start)/2.0, (slice_z.stop + slice_z.start)/2.0])
        #
        # Now run local maxima finding and then watershed
        roi = ge_data[slice_x, slice_y, slice_z]
        roi_original = ge_data[slice_x, slice_y, slice_z]

        markers = np.zeros_like(roi)
        # Find local maxima
        max_points = peak_local_max(roi, min_distance=(min_peak_separation),
                                    threshold_rel=0.05, exclude_border=False, indices=False)
        max_points[roi < int_scale_factor*cfg.fit_grains.threshold] = 0
        max_points = np.nonzero(max_points)
        #
        for max_x, max_y, max_z, max_id in zip(max_points[0], max_points[1], max_points[2],
                                               range(len(max_points[0]))):
           markers[max_x][max_y][max_z] = max_id+1
           # Store global x, y, z for the maxima and the intensity
           max_points_global.append([max_x + slice_x.start, max_y + slice_y.start, max_z + slice_z.start, roi[max_x][max_y][max_z]])
        # Run watershed now
        labels = watershed(-roi, markers, mask=(roi>0.1*np.amax(roi)))
        # Done watershed
        # This looks like a legit spot. Add to blobs array
        blobs.append(GEBlob(slice_x, slice_y, slice_z, label_num, blob_size, max_points))

   ge_labeled = np.amax(label_ge, 0)

   return {'blobs': blobs, 'label_ge': label_ge, 'blob_centroids': blob_centroids, 'local_maxima': max_points_global}
#--

# A blob is a set of pixels in an image that are connected to each other
class GEBlob:
    def __init__(self, slice_x, slice_y, slice_z,
                 blob_label, blob_size, max_points):
       # slice_* defines the bounding box of a blob in the GE2 data
       self.slice_x    = slice_x
       self.slice_y    = slice_y
       self.slice_z    = slice_z
       # Blob id
       self.blob_label = blob_label
       # Number of pixels
       self.blob_size  = blob_size
       # Local maxima in the blob ([[x], [y], [z]])
       # xyz for the local maxima are w.r.t to the blob
       # To get global xyz for the maxima, add slice_x.start
       # etc. to the xyz
       self.max_points = max_points
#--
# An object for all the GE2 pre-processing routines
class GEPreProcessor:
    '''
        Pre-processing on GE files to smooth data,
        detect local maxima etc

        cfg = configuration dictionary read from a heXRD config file
        logger = a logger object to log progress/comments
    	gauss_sigma = std deviation for Gaussian smoothing on GE data
    	ge_data = GE raw data
    	ge_smooth_data = GE smoothed data with Gaussian
        ge_labeled_data = GE data with connected components labeled
        number_of_labels = Number of connected components in GE
        min_blob_size = min size of connected objects
        int_scale_factor = GE intensity is scaled by this.
        min_peak_separation = min_peak_separation
        blobs = An array of GEBlob (information about blobs)
    '''

    def __init__(self, cfg, logger):
    	self.cfg                 = cfg                      # An open hexrd config file object
        self.logger              = logger                   # An open logger object
    	self.ge_data             = []                       # GE2 image data
    	self.ge_smooth_data      = []                       # Above + smoothed using Gauss
        self.ge_labeled_data     = []                       # Above + connected components labeled
        self.number_of_labels    = 0                        # number of labels = number of blobs in the data
        self.min_blob_size       = cfg.get('pre_processing')['min_blob_size'] # Blobs smaller than this are removed (user input)
        self.int_scale_factor    = 1                        # A scale factor for intensity (auto-calculated)
        self.min_peak_separation = cfg.get('pre_processing')['min_peak_separation'] # Minimum separation in the local maxima (user input)
        self.blobs               = []                       # An array of blob objects
	self.max_points          = []                       # An array of local maxima coordinates in the blobs
        self.omega_start         = []                       # In the parallelized frame data, start omega number for each portion

        return
    #--
    def load_data(self):
        '''
            Read the config file and load appropriate GE2
            frames.
        '''
        cfg          = self.cfg
        logger       = self.logger
        omega_start  = self.omega_start

        # process the data
        pd, reader, detector = initialize_experiment(cfg)
        n_frames = reader.getNFrames()
        logger.info("Reading %d frames of data, storing values > %.1f",
                    n_frames, cfg.get('pre_processing')['ge_reader_threshold'])
        # Loop over all frames and save them in a 3D array
        frame_list = []
        for i in range(n_frames):
            frame = reader.read()
            #omega = reader.getFrameOmega()
            frame_list.append(frame)
        # Turn the frame array into a Numpy array
        frame_list = np.array(frame_list)
        # Remove low intensity noise
        frame_list[frame_list < cfg.get('pre_processing')['ge_reader_threshold']] = 0
        # Scale the intensity to 16k
        int_scale_factor = float(2**14)/float(np.amax(frame_list))
        frame_list = frame_list*int_scale_factor
        if cfg.get('pre_processing')['print_diag_images']:
        	# Flatten along omega and write the frame array to an image
        	write_image('slice.png', np.amax(frame_list, axis=0), vmin=0)
        # Split the frame array into chunks for multiprocessing
        num_cores = cfg.multiprocessing
        frame_list_split = np.array_split(frame_list, num_cores, axis=0)
        ge_data_ang_red = ()
        omega_start.append(0)
        for array_piece in frame_list_split:
           ge_data_ang_red = ge_data_ang_red + (array_piece,)
           omega_start.append(np.shape(array_piece)[0])

        omega_start.pop()
        omega_start = np.cumsum(omega_start)
        logger.info("Finished reading frames")

        self.ge_data          = frame_list
        self.int_scale_factor = int_scale_factor
	self.ge_data_ang_red  = ge_data_ang_red
        self.omega_start      = omega_start
        self.input_data_shape = np.shape(frame_list)

        return frame_list
    #--
    def find_blobs(self):
        '''
            Find connected componnets (spots)
        '''
        min_size            = self.min_blob_size
        ge_data_smooth      = self.ge_smooth_data
        cfg                 = self.cfg
        logger              = self.logger
        int_scale_factor    = self.int_scale_factor
        ge_data_ang_red     = self.ge_data_ang_red
        min_peak_separation = self.min_peak_separation
        omega_start         = self.omega_start

        num_cores = cfg.multiprocessing
        logger.info("Starting spot finding with %d cores", num_cores)
        blobs_mp_output = Parallel(n_jobs=num_cores, verbose=5, max_nbytes=1e6)(delayed(find_blobs_mp)(ge_data, int_scale_factor, min_size, min_peak_separation, cfg) for ge_data in ge_data_ang_red)

        logger.info("Finished multiprocessing spot finding algorithm")

        blobs = []
        label_ge = []
        blob_centroids_oxy = []
        local_maxima_oxy = []
        local_maxima_oxyi = []
        for blobs_mp_output_i, omega_start_i in zip(blobs_mp_output, omega_start):
           for blob_i in blobs_mp_output_i['blobs']:
              blobs.append(blob_i)
           #
           for maxima_o, maxima_x, maxima_y, max_intensity in blobs_mp_output_i['local_maxima']:
              if max_intensity > (0.05*float(2**14)):
                 local_maxima_oxyi.append([maxima_o + omega_start_i, maxima_x, maxima_y, max_intensity])
                 local_maxima_oxy.append([maxima_o + omega_start_i, maxima_x, maxima_y])
           #
           label_ge.append(blobs_mp_output_i['label_ge'])
           

        logger.info("Clustering spots")
        # Cluster the local minima
        local_maxima_oxy = np.array(local_maxima_oxy)
        local_maxima_oxyi = np.array(local_maxima_oxyi)
        db = DBSCAN(eps=2.5, min_samples=1).fit(local_maxima_oxy)
        local_maxima_labels = db.labels_
        #
        o_sum = np.bincount(local_maxima_labels, weights=local_maxima_oxyi[:, 0])
        x_sum = np.bincount(local_maxima_labels, weights=local_maxima_oxyi[:, 1])
        y_sum = np.bincount(local_maxima_labels, weights=local_maxima_oxyi[:, 2])
        i_sum = np.bincount(local_maxima_labels, weights=local_maxima_oxyi[:, 3])
        l_sum = np.bincount(local_maxima_labels)
        #
        local_maxima_oxyi_clustered = zip(np.divide(o_sum, l_sum), np.divide(x_sum, l_sum), np.divide(y_sum, l_sum), np.divide(i_sum, l_sum))
        local_maxima_xy = zip(np.divide(x_sum, l_sum), np.divide(y_sum, l_sum))

        blob_centroids = []
        for blob in blobs:
           if blob.blob_size > cfg.get('pre_processing')['min_blob_size']:
              blob_centroids.append([(blob.slice_y.start + blob.slice_y.stop)/2.0, (blob.slice_z.start + blob.slice_z.stop)/2.0])
              blob_centroids_oxy.append([(blob.slice_x.start + blob.slice_x.stop)/2.0, (blob.slice_y.start + blob.slice_y.stop)/2.0, (blob.slice_z.start + blob.slice_z.stop)/2.0])

        self.blobs                       = blobs
        self.blob_centroids_oxy          = blob_centroids_oxy
        self.local_maxima_oxyi           = local_maxima_oxyi
        self.local_maxima_oxyi_clustered = local_maxima_oxyi_clustered
        self.local_maxima_clusters       = db

        logger.info("Finished clustering")
        logger.info("Found %d connected components", len(blobs))
        logger.info("Found %d spots", np.shape(local_maxima_oxyi)[0])

	if cfg.get('pre_processing')['print_diag_images']:
        	logger.info("Writing diagnostic images")
        	# Superimpose the centroids of blobs and the local maxima on the original data 
        	# and write to an image file.
       		write_image('slice_blob_centroids.png', np.amax(self.ge_data, axis=0), pts=blob_centroids)
        	write_image('slice_local_maxima.png', np.amax(self.ge_data, axis=0), pts=local_maxima_xy)

        if(cfg.get('pre_processing')['print_spots_info']):
           # Because I want to pretty print columns
           template = "{0:>12}{1:>12}{2:>12}{3:>12}"
           print template.format("Omega", "X", "Y", "Intensity")

           template = "{0:12.2f}{1:12.2f}{2:12.2f}{3:12.3f}"
           for o, x, y, i in local_maxima_oxyi_clustered:
              print template.format(o, x, y, i)

        if cfg.get('pre_processing')['print_ge']:
           logger.info("Writing GE2 files")
           # Synthesize a GE2 file based on the IDed spots
           frames_synth = np.zeros(self.input_data_shape)
           for o, x, y, i in local_maxima_oxyi_clustered:
              frames_synth[int(round(o)), int(round(x)), int(round(y))] = i

           frames_synth = ndimage.morphology.grey_dilation(frames_synth, size=2)
           write_ge2('synth_spots.ge2', frames_synth)
        else:
           logger.info("Skipped writing GE2 files")

        return label_ge
    #--

    def find_local_maxima(self):
        '''
            Find local maxima in each of the detected blobs
        '''
        ge_labeled_data     = self.ge_labeled_data
        ge_data_ang_red     = self.ge_data_ang_red
        cfg                 = self.cfg
        logger              = self.logger
        int_scale_factor    = self.int_scale_factor
        number_of_labels    = self.number_of_labels
        min_peak_separation = self.min_peak_separation
        blobs               = self.blobs

        logger.info("Detecting local maxima")

        num_cores = cfg.multiprocessing
        spot_max_points = Parallel(n_jobs=num_cores, verbose=5)(delayed(get_local_maxima)(blob, ge_data, min_peak_separation, cfg, int_scale_factor) for blob, ge_data in zip(blobs, ge_data_ang_red))

        for pt, layer_num in zip(spot_max_points, range(num_cores)):
           for x, y, z in pt:
              print layer_num, x, y, z

	self.max_points = spot_max_points

        return blobs
