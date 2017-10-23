#!/usr/bin/python
# -*- coding: utf-8 -*-

import astropy

from astropy.io import fits
from astropy.table import Table

catalog_n = input('Catalog number: ')
i_alpha = input('ALPHA_J2000: ')
i_delta = input('DELTA_J2000: ')

o_file = fits.open('/home/sgongora/Documents/CarpetaCompartida/full_10_1.2_2.5_0.64_20-21_1.cat')
o_cat = Table(o_file[2].data).to_pandas()


tolerance = 0.0002

o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
o_df = o_df[o_df['ALPHA_J2000'] + tolerance > i_alpha]
o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - tolerance]
o_df = o_df[o_df['DELTA_J2000'] + tolerance > i_delta]
o_df = o_df[i_delta > o_df['DELTA_J2000'] - tolerance]

print o_df
