# !/usr/bin/python
# -*- coding: utf-8 -*-

""" Extracts SSOs regions from CCDs or filtered catalogs

Versions:
- 0.1

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""

from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pandas import concat, DataFrame, read_csv

from misc import extract_settings_elvis, get_cats_elvis_d
from sys import argv

from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import array, isnan, nan, nanmean
from pandas import concat, DataFrame, read_csv, Series

from regions import Create_regions


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_dither(catalog_n):
    """

    :param catalog_n:
    :return: dither_n
    """
    dither_n = 0

    if catalog_n <= 36:
        dither_n = 1
    elif 36 < catalog_n <= 72:
        dither_n = 2
    elif 72 < catalog_n <= 108:
        dither_n = 3
    elif 108 < catalog_n <= 144:
        dither_n = 4

    return dither_n


def check_source(o_df, i_alpha, i_delta, keys):
    """

    :param o_df:
    :param i_alpha:
    :param i_delta:
    :param keys:
    :return:
    """
    prfs_d = extract_settings_elvis()

    o_df = o_df[o_df[keys[0]] + prfs_d['tolerance']*2 > i_alpha]
    o_df = o_df[i_alpha > o_df[keys[0]] - prfs_d['tolerance']*2]
    o_df = o_df[o_df[keys[1]] + prfs_d['tolerance']*2 > i_delta]
    o_df = o_df[i_delta > o_df[keys[1]] - prfs_d['tolerance']*2]

    return o_df


class FactorsScampPerformance:

    def __init__(self):
        """ Me da los valores de salida de todos las estrellas  y galaxias
        presentes en filt 3 obtenidos o no

        """
        self.filter_p_number = 3  # First one with enough data for statistics
        self.prfs_d = extract_settings_elvis()

        ccds = True
        filtered = False
        scamp = False

        input_df = read_csv('cats/cat_clean_ssos.csv', index_col=0)
        filt_cat = self.gets_filtered_catalog()  # Gets data from filtered

        if ccds:
            cats_d = self.extract_cats()
            self.extract_stats_ccds(cats_d, input_df, filt_cat)
        elif filtered:
            self.extract_stats_filt(filt_cat, input_df)
        elif scamp:
            pass
            # self.extract_stats_scamp(input_df)
        else:
            pass

    def extract_cats(self):
        """

        :return:
        """
        cats_d = {}
        for dither_ in range(1, 5, 1):
            cat_list = get_cats_elvis_d(dither_)
            cats_d[dither_] = {}
            for cat_ in cat_list:
                cat = fits.open('{}/{}'.format(self.prfs_d['fits_dir'], cat_))
                cat_df = Table(cat[2].data).to_pandas()
                cats_d[dither_][cat_] = cat_df

        return cats_d

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt__{}.csv'.format(self.filter_p_number)
        # filter_n = 'filt__full_1.csv'
        filter_o_n = '{}/{}'.format(self.prfs_d['filtered'], filter_n)

        print('Opens filtered catalogue {}'.format(filter_o_n))
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def extract_stats_ccds(self, cats_d, input_df, filt_cat):
        """

        :param cats_d:
        :param input_df:
        :param filt_cat:
        :return:
        """
        # Unique sources (?)
        unique_sources = list(set(input_df['SOURCE'].tolist()))

        test_dict = {'PM': [], 'A_IMAGE': [], 'B_IMAGE': [], 'CLASS_STAR': []}

        print('total {}'.format(len(unique_sources)))
        for idx_source_, source_ in enumerate(unique_sources):
            print('idx {}'.format(idx_source_))
            source_df = input_df[input_df['SOURCE'].isin([source_])]

            # Loops over CCD catalogues
            for i, row in enumerate(source_df.itertuples(), 1):
                dither_df = source_df[source_df['DITHER'].isin([row.DITHER])]
                print('dither {}'.format(row.DITHER))
                i_alpha = float(dither_df['RA'].iloc[0])
                i_delta = float(dither_df['DEC'].iloc[0])

                pm = ''
                a_image = ''
                b_image = ''
                class_star = ''
                test = False
                for cat_ in cats_d[row.DITHER]:
                    out_df = check_source(cats_d[row.DITHER][cat_],
                                          i_alpha, i_delta,
                                          keys=['ALPHA_J2000', 'DELTA_J2000'])

                    if out_df.empty is not True:
                        alpha = float(out_df['ALPHA_J2000'].iloc[0])
                        delta = float(out_df['DELTA_J2000'].iloc[0])
                        pm_df = check_source(filt_cat, alpha, delta,
                                             keys=['ALPHA_J2000',
                                                   'DELTA_J2000'])
                        print(pm_df['PM'].size)
                        a_image = out_df['A_IMAGE'].iloc[0]
                        b_image = out_df['B_IMAGE'].iloc[0]
                        class_star = out_df['CLASS_STAR'].iloc[0]
                        test = True

                if test:
                    test_dict['PM'].append(pm)
                    test_dict['A_IMAGE'].append(a_image)
                    test_dict['B_IMAGE'].append(b_image)
                    test_dict['CLASS_STAR'].append(class_star)

        pm_list = Series(test_dict['PM'], name='PM')
        a_image_list = Series(test_dict['A_IMAGE'], name='A_IMAGE')
        b_image_list = Series(test_dict['B_IMAGE'], name='B_IMAGE')
        class_star_list = Series(test_dict['CLASS_STAR'], name='CLASS_STAR')

        positions_table = concat([pm_list, a_image_list,
                                  b_image_list, class_star_list], axis=1)
        positions_table.to_csv('catalogue.csv')

    def extract_stats_filt(self, filt_cat, input_df):
        """

        :param filt_cat:
        :param input_df:
        :return:
        """
        # Unique sources (?)
        unique_sources = list(set(input_df['SOURCE'].tolist()))

        ok = 0
        no = 0
        filt_cat_size = filt_cat['SOURCE_NUMBER'].size
        total_input_size = input_df['SOURCE'].size
        print('sources in input catalogue: {}'.format(total_input_size))
        print('sources in filtered catalogue {}'.format(filt_cat_size))

        for source_ in unique_sources:
            source_df = input_df[input_df['SOURCE'].isin([source_])]
            print(source_df.columns)
            for i, row in enumerate(source_df.itertuples(), 1):
                dither_df = source_df[source_df['DITHER'].isin([row.DITHER])]
                i_alpha = float(dither_df['RA'].iloc[0])
                i_delta = float(dither_df['DEC'].iloc[0])
                out_df = check_source(filt_cat, i_alpha, i_delta,
                                      keys=['ALPHA_J2000', 'DELTA_J2000'])

                if out_df.empty:
                    no += 1
                else:
                    ok += 1

        print('detected: {}'.format(ok))
        print('not detected: {}'.format(no))


if __name__ == "__main__":
    FactorsScampPerformance()
