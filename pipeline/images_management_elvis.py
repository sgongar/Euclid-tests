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


def get_position(order):
    """

    :param order:
    :return:
    """
    order_d = {0: [1, 6], 1: [2, 6], 2: [3, 6], 3: [4, 6],
               4: [5, 6], 5: [6, 6], 6: [1, 5], 7: [2, 5],
               8: [3, 5], 9: [4, 5], 10: [5, 5], 11: [6, 5],
               12: [1, 4], 13: [2, 4], 14: [3, 4], 15: [4, 4],
               16: [4, 5], 17: [4, 6], 18: [1, 3], 19: [2, 3],
               20: [3, 3], 21: [4, 3], 22: [5, 3], 23: [6, 3],
               24: [1, 2], 25: [2, 2], 26: [3, 2], 27: [4, 2],
               28: [5, 2], 29: [6, 2], 30: [1, 1], 31: [2, 1],
               32: [3, 1], 33: [4, 1], 34: [5, 1], 35: [6, 1]}

    coords = 'x{}_y{}'.format(order_d[order][0], order_d[order][1])

    return coords


def extract_quadrants(fits_file):
    """

    :param fits_file:
    :return:
    """
    images_idxs = np.arange(1, 144, 4)
    hdu_list = fits.open(fits_file)

    for order, idx in enumerate(images_idxs):
        for quadrant in range(1, 5, 1):
            coords = get_position(order)
            name = 'CCD_{}_q{}_d{}.fits'.format(coords, quadrant,
                                                fits_file[-6:-5])
            fits.writeto(name, hdu_list[idx].data, header=hdu_list[idx].header)
