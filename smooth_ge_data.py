import copy
import logging
import os
import sys
import time
import warnings

import yaml

from hexrd import config
from ge_processor.ge_pre_processor import *

if __name__ == '__main__':
    # Read args
    if len(sys.argv) < 2:
        print 'USAGE: python load_data.py config.yml'
        sys.exit(1)
    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print 'USAGE: python load_data.py config.yml'
        sys.exit(1)
    else:
        cfg_file = sys.argv[1]

    log_level = logging.DEBUG
    logger = logging.getLogger('hexrd')
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    cf = logging.Formatter('%(asctime)s - %(message)s', '%y-%m-%d %H:%M:%S')
    ch.setFormatter(cf)
    logger.addHandler(ch)
    # load the configuration settings
    cfg = config.open(cfg_file)
    # cfg is a list. We only need the first cfg data.
    cfg = cfg[0]

    gepp = GEPreProcessor(cfg=cfg, logger=logger, gauss_sigma=3, min_blob_size=729)

    logger.info('=== begin image-smoothing ===')

    gepp.load_data()
    gepp.smooth_data()
    gepp.find_blobs()
    gepp.find_local_maxima()

    # ge_data = load_data(cfg)
    # smooth_data(cfg, ge_data)
