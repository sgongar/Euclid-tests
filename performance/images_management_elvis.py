#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release. Split from ssos_catalog_creation.py

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    *

*GNU Terry Pratchett*

"""
from math import cos, sin
from multiprocessing import Process
from sys import stdout

from astropy.io import fits
from astropy.table import Table
from numpy import median
from pandas import concat, DataFrame, read_csv

from misc import extract_settings_elvis, get_fits, get_fits_limits

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_borders():
    """

    :return:
    """
    prfs_d = extract_settings_elvis()
    borders_d = {}
    for dither_ in range(1, 5, 1):
        borders_d[dither_] = {}
        fits_list = get_fits(dither_)
        for fits_ in fits_list:
            fits_loc = '{}/{}'.format(prfs_d['fits_dir'], fits_)
            borders_d[dither_][fits_] = get_fits_limits(fits_loc)

    return borders_d
