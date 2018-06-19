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
from astropy.io import fits
from astropy.table import Table
from pandas import concat, Series

from misc import extract_settings_elvis, get_cats_elvis_d

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def create_regions():
    """

    :return:
    """
    for dither_ in range(1, 5, 1):
        cat_list = get_cats_elvis_d(dither_)
        for cat_ in cat_list:
            catalog = fits.open('{}/{}'.format(prfs_d['fits_dir'], cat_))
            catalog_data = Table(catalog[2].data).to_pandas()

            alpha_list = Series(catalog_data['ALPHA_J2000'].tolist(),
                                name='ALPHA_J2000')
            delta_list = Series(catalog_data['DELTA_J2000'].tolist(),
                                name='DELTA_J2000')
    
            positions_table = concat([alpha_list, delta_list], axis=1)
            positions_table.to_csv('{}.reg'.format(cat_[:-4]),
                                   index=False, header=False, sep=" ")


if __name__ == "__main__":
    prfs_d = extract_settings_elvis()

    create_regions()
