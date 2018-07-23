# !/usr/bin/python
# -*- coding: utf-8 -*-

""" Gets
   - factors from magnitude bins

Versions:
- 0.1

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""

from astropy.io import fits
from astropy.table import Table
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from numpy import array, float64
from pandas import concat, DataFrame, read_csv, Series
from pyds9 import DS9
import statsmodels.api as sm


from misc import extract_settings_elvis, setting_logger


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def main():
    catalog_loc = '/home/sgongora/Dev/Euclid-tests/performance/stats/total.csv'

    for pm_ in prfs_dict['pms']:
        print(pm_)


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    main()
