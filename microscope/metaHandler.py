"""
a class to collect and calculate microscope metadata

"""
import numpy as np
from collections import OrderedDict
from astropy.stats import binom_conf_interval
from analysis import Analysis, BOOL
import fitCurve as fc