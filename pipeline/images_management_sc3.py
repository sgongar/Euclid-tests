#!/usr/bin/python
# -*- coding: utf-8 -*-

"""


Todo:
    * Improve log messages
    * Improve usability
"""


from astropy.io import fits
from numpy import arange
from collections import Counter
from decimal import Decimal
from math import hypot
from multiprocessing import cpu_count
from os import listdir, remove, rename, makedirs, path, mkdir
from platform import platform

from ConfigParser import ConfigParser
import numpy as np
from pandas import Series
import statsmodels.api as sm

from errors import BadSettings
from logging import getLogger, config


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def extract_images(fits_file):
    """

    :param fits_file:
    :return:
    """
    images_idxs = arange(1, 109, 3)
    hdu_list = fits.open(fits_file)

    for order, idx in enumerate(images_idxs):
        fits.writeto('{}_{}.fits'.format(fits_file[:-5], order),
                     hdu_list[idx].data, header=hdu_list[idx].header)


def extract_flags(fits_file):
    """

    :param fits_file:
    :return:
    """
    images_idxs = arange(3, 109, 3)
    hdu_list = fits.open(fits_file)

    for order, idx in enumerate(images_idxs):
        fits.writeto('{}_f{}.fits'.format(fits_file[:-5], order),
                     hdu_list[idx].data, header=hdu_list[idx].header)