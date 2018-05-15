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
    order_d = {0: [1, 6], 1: [2, 6], 2: [3, 6], 3: [4, 6], 4: [5, 6], 5: [6, 6],
               6: [1, 5], 7: [2, 5], 8: [3, 5], 9: [4, 5], 10: [5, 5],
               11: [6, 5], 12: [1, 4], 13: [2, 4], 14: [3, 4], 15: [4, 4],
               16: [5, 4], 17: [6, 4], 18: [1, 3], 19: [2, 3], 20: [3, 3],
               21: [4, 3], 22: [5, 3], 23: [6, 3], 24: [1, 2], 25: [2, 2],
               26: [3, 2], 27: [4, 2], 28: [5, 2], 29: [6, 2], 30: [1, 1],
               31: [2, 1], 32: [3, 1], 33: [4, 1], 34: [5, 1], 35: [6, 1]}

    coords = 'x{}_y{}'.format(order_d[order][0], order_d[order][1])

    return coords


def create_ccds(proc, fits_dir, fpa_dir, fpa_file):
    """

    :param fits_dir:
    :param fpa_dir:
    :param fpa_file:
    :return:
    """
    quadrants_d = {}
    dither = fpa_file[-6:-5]

    images_idxs = np.arange(1, 144, 4)
    hdu_list = fits.open('{}/{}'.format(fpa_dir, fpa_file))

    for order, idx in enumerate(images_idxs):
        coords = get_position(order)
        quadrants_l = []

        for quadrant in range(1, 5, 1):
            """
            print('order {} - quadrant {} - dither {}'.format(order, quadrant,
                                                              dither))
            name = 'CCD_{}_q{}_d{}.fits'.format(fits_dir, coords, quadrant,
                                                fits_file[-6:-5])
            """
            quadrants_l.append(hdu_list[idx])
            # quadrants_l.append(hdu_list[idx].data)

        quadrant_name = 'CCD_{}_d{}'.format(coords, dither)
        quadrants_d[quadrant_name] = quadrants_l

    for key_ in quadrants_d.keys():
        for quadrant_ in range(0, len(quadrants_d[key_]), 1):
            create_ccd(quadrants_d[key_], key_)


def create_ccd(quadrants, key_):
    """

    :param quadrants:
    :return:
    """
    print(key_)

    # hdus = fits.open(args.inputFile)
    prex = quadrants[0].header['PRESCANX']
    ovrx = quadrants[0].header['OVRSCANX']
    try:
        ovry = quadrants[0].header['OVRSCANY']
    except:
        ovry = 0

    # Code for injected lines, when they will be implemented in ELViS
    for i in range(0, len(quadrants), 1):
        img = np.zeros([4132, 4096])
        # Quadrant E
        if i == 0:
            img[:2066, :2048] = quadrants[i].data[:, prex:-ovrx]
            #
            # Correct reference pixel coordinate
            hdr = quadrants[i].header
            hdr.set('CRPIX1', hdr['CRPIX1'] - prex)
        # Quadrant F
        if i == 1:
            img[:2066, 2048:] = quadrants[i].data[:, ovrx:-prex]
        # Quadrant G
        if i == 2:
            img[2066:, :2048] = quadrants[i].data[:, prex:-ovrx]
        # Quadrant H
        if i == 3:
            img[2066:, 2048:] = quadrants[i].data[:, ovrx:-prex]
            #
            # Save to FITS file
            outputfits = fits.HDUList()
            outputfits.append(fits.PrimaryHDU())
            outputfits.append(fits.ImageHDU(data=img, header=hdr))
            outputfits.writeto('{}.fits'.format(key_), overwrite=True)

    return True
