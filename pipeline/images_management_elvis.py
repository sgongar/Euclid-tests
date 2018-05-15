#!/usr/bin/python
# -*- coding: utf-8 -*-

"""


Todo:
    * Improve log messages
    * Improve usability
"""

from astropy.io import fits
import numpy as np

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def extract_quadrants(fits_file):
    """

    :param fits_file:
    :return:
    """
    images_idxs = np.arange(1, 144, 4)
    hdu_list = fits.open(fits_file)

    for order, idx in enumerate(images_idxs):
        for quadrant in range(1, 5, 1):
            fits.writeto('{}_{}_{}.fits'.format(fits_file[:-5],
                                                quadrant, order),
                         hdu_list[idx].data, header=hdu_list[idx].header)
