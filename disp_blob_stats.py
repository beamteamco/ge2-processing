import sys

import cPickle as pickle
import pprint

blobs = pickle.load(open('ge_blobs.cpl', 'rb'))

for blob in blobs:
	print blob.blob_label, blob.blob_size

