"""Microscope Processing
Jonathan Mortlock 25/06/20

"""

import os
import sys
import numpy as np
from collections import OrderedDict
import time
from scipy.signal import find_peaks
from scipy.stats import norm
from skimage.filters import threshold_minimum
from astropy.stats import binom_conf_interval
from imageanalysis.analysis import Analysis, BOOL
import logging
logger = logging.getLogger(__name__)



class image_handler(Analysis):
    """Class to handle image metadata
    (NB Slightly different implementation to Tweezers)
    Inherits the types and stats dictionaries, and reset_arrays, 
    load, and save methods from Analysis.    
    """
    def __init__(self):
        super().__init__()
        self.types = OrderedDict([('File ID', int), # number ID of image
            ('Atom Number', int), # Total occupied
            ('Lattice 00',float),
            ('Lattice 01',float),
            ('Lattice 10',float),
            ('Lattice 11',float),
            ('Offset 0',float),
            ('Offset 1',float),
            ('RL iterations',int),
            ('Threshold',float),
            ('Include', BOOL)])# whether to include in further analysis
        self.stats = OrderedDict([(key, []) for key in self.types.keys()]) 
        
        self.delim = ' '                # delimieter to use when opening image files
        self.bias = 697                 # bias offset from EMCCD
        self.peak_indexes = [0,0]       # indexes of peaks in histogram
        self.peak_heights = [0,0]       # heights of peaks in histogram
        self.peak_widths  = [0,0]       # widths of peaks in histogram
        self.peak_centre  = [0,0]       # peak position in counts in histogram
        self.fidelity     = 0           # fidelity of detecting atom presence
        self.err_fidelity = 0           # error in fidelity
        # self.mask      = np.zeros((1,1))# normalised mask to apply to image for ROI
        self.xc        = 1              # ROI centre x position 
        self.yc        = 1              # ROI centre y position
        self.roi_size  = 1              # ROI length in pixels. default 1 takes top left pixel
        self.pic_width = 512            # number of pixels along horizontal axis of an image
        self.pic_height= 512            # number of pixels along vertical axis of an image
        self.thresh    = 1              # initial threshold for atom detection
        self.fid       = 0              # file ID number for the next image
        self.ind       = 0              # number of images processed
        self.im_vals   = np.array([])   # the data from the last image is accessible to an image_handler instance
        self.bin_array = []        # if bins for the histogram are supplied, plotting can be faster
