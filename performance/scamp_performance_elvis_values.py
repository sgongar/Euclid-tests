# !/usr/bin/python
# -*- coding: utf-8 -*-

""" Gets
    - 'median_a_image'
    - 'median_erra_image'
    - 'median_b_image'
    - 'median_errb_image'
    - 'median_class_star'
    - 'ellipticity'
    - 'median_mag_iso'
    - 'median_magerr_iso'
    - 'median_flux_iso'
    from scamp's output. Saves them to different csv files.

Versions:
- 0.1

Todo:
    * Improve log messages
    * Get out check_source

*GNU Terry Pratchett*
"""
from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pandas import concat, DataFrame, read_csv

from misc import extract_settings_elvis, setting_logger

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


def create_output_dicts():
    # prfs_d = extract_settings_elvis()

    # Creates a dictionary
    ssos_d = {'median_a_image': [], 'median_erra_image': [],
              'median_b_image': [], 'median_errb_image': [],
              'median_class_star': [],
              'ellipticity': [], 'median_mag_iso': [],
              'median_magerr_iso': [], 'median_flux_iso': []}
    # for idx in range(0, len(prfs_d['pms']), 1):
    #     ssos_d['median_mag_iso'].append([])
    #     ssos_d['median_magerr_iso'].append([])
    #     ssos_d['median_a_image'].append([])
    #     ssos_d['median_erra_image'].append([])
    #     ssos_d['median_b_image'].append([])
    #     ssos_d['median_errb_image'].append([])
    #     ssos_d['median_class_star'].append([])
    #     ssos_d['ellipticity'].append([])
    #     ssos_d['output_pm'].append([])
    #     ssos_d['median_flux_iso'].append([])

    stars_d = {'median_a_image': [], 'median_erra_image': [],
               'median_b_image': [], 'median_errb_image': [],
               'median_class_star': [],
               'ellipticity': [], 'median_mag_iso': [],
               'median_magerr_iso': [], 'median_flux_iso': []}
    # for idx in range(0, len(prfs_d['pms']), 1):
    #     stars_d['median_mag_iso'].append([])
    #     stars_d['median_magerr_iso'].append([])
    #     stars_d['median_a_image'].append([])
    #     stars_d['median_erra_image'].append([])
    #     stars_d['median_b_image'].append([])
    #     stars_d['median_errb_image'].append([])
    #     stars_d['median_class_star'].append([])
    #     stars_d['ellipticity'].append([])
    #     stars_d['output_pm'].append([])
    #     stars_d['median_flux_iso'].append([])

    galaxies_d = {'median_a_image': [], 'median_erra_image': [],
                  'median_b_image': [], 'median_errb_image': [],
                  'median_class_star': [],
                  'ellipticity': [], 'median_mag_iso': [],
                  'median_magerr_iso': [], 'median_flux_iso': []}
    # for idx in range(0, len(prfs_d['pms']), 1):
    #     galaxies_d['median_mag_iso'].append([])
    #     galaxies_d['median_magerr_iso'].append([])
    #     galaxies_d['median_a_image'].append([])
    #     galaxies_d['median_erra_image'].append([])
    #     galaxies_d['median_b_image'].append([])
    #     galaxies_d['median_errb_image'].append([])
    #     galaxies_d['median_class_star'].append([])
    #     galaxies_d['ellipticity'].append([])
    #     galaxies_d['output_pm'].append([])
    #     galaxies_d['median_flux_iso'].append([])

    lost_d = {'median_a_image': [], 'median_erra_image': [],
              'median_b_image': [], 'median_errb_image': [],
              'median_class_star': [],
              'ellipticity': [], 'median_mag_iso': [],
              'median_magerr_iso': [], 'median_flux_iso': []}

    output_d = {'ssos': ssos_d, 'stars': stars_d, 'galaxies': galaxies_d,
                'lost': lost_d}

    return output_d


def check_source(o_df, i_alpha, i_delta):
    """

    :param o_df:
    :param i_alpha:
    :param i_delta:
    :return:
    """
    prfs_d = extract_settings_elvis()

    o_df = o_df[o_df['ALPHA_J2000'] + prfs_d['tolerance'] > i_alpha]
    o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - prfs_d['tolerance']]
    o_df = o_df[o_df['DELTA_J2000'] + prfs_d['tolerance'] > i_delta]
    o_df = o_df[i_delta > o_df['DELTA_J2000'] - prfs_d['tolerance']]

    return o_df


class TotalScampPerformance:

    def __init__(self):
        """ Me da los valores de salida de todos las estrellas  y galaxias
        presentes en filt 3 obtenidos o no

        """
        self.filter_p_number = 3  # First one with enough data for statistics
        self.prfs_d = extract_settings_elvis()
        self.data_d = create_output_dicts()

        self.save = True

        logger_name = 'scamp_performance'  # Set as desired
        self.logger = setting_logger(self.prfs_d, logger_name)

        filt_cat = self.gets_filtered_catalog()  # Gets data from filtered
        input_df = self.gets_data()  # Gets data from catalogs
        self.splits_data(filt_cat, input_df)  # Splits due type

        print(self.data_d.keys())
        print(type(self.data_d['ssos']))

        # Saves data
        for key_ in ['ssos', 'stars', 'galaxies']:
            out_df = DataFrame(self.data_d[key_])
            out_df.to_csv('{}.csv'.format(key_))

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt__{}.csv'.format(self.filter_p_number)
        filter_o_n = '{}/{}'.format(self.prfs_d['filtered'], filter_n)

        print('Opens filtered catalogue {}'.format(filter_o_n))
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def gets_data(self):
        """ Creates an input dictionary. Each key contains SSOs' information
        for each dither.

        :return: input_dict
        """
        # # For now we only have data for dither 1
        # input_df = {1: {}, 2: {}, 3: {}, 4: {}}
        #
        # columns = ['ALPHA_J2000', 'DELTA_J2000', 'pm', 'pa',
        #            'mag', 'mag', 'mag', 'mag']
        # input_cat = read_csv('{}/SSO_Cat.txt'.format(self.prfs_d['references']),
        #                      delim_whitespace=True, header=None,
        #                      names=columns)
        # print('Opens SSO catalogue {}'.format('{}/SSO_Cat.txt'.format(self.prfs_d['references'])))
        # input_cat = input_cat[['ALPHA_J2000', 'DELTA_J2000', 'pm', 'mag']]
        #
        # # Merge output
        # reject_l = {1: [], 2: [], 3: [], 4:[]}
        #
        # for d in range(1, 2, 1):
        #     cat_name = '{}/coadd_{}.cat'.format(self.prfs_d['references'], d)
        #     print('Opens reference catalogue {}'.format(cat_name))
        #     coadd_cat = fits.open(cat_name)
        #     cat_data = Table(coadd_cat[2].data).to_pandas()
        #
        #     for idx_source, source_ in enumerate(cat_data['NUMBER']):
        #         source_df = cat_data[cat_data['NUMBER'].isin([source_])]
        #         alpha = source_df['ALPHA_J2000'].iloc[0]
        #         delta = source_df['DELTA_J2000'].iloc[0]
        #
        #         test_df = check_source(input_cat, alpha, delta)
        #
        #         if test_df.empty is not True:
        #             reject_l[d].append(source_)
        #
        #     print('Removing SSOs from catalog')
        #     input_df[d]['SSOs'] = cat_data[cat_data['NUMBER'].isin(reject_l)]
        #     fixed_data = cat_data[~cat_data['NUMBER'].isin(reject_l)]
        #     input_df[d]['stars'] = fixed_data[fixed_data['CLASS_STAR'] > 0.95]
        #     input_df[d]['galaxies'] = fixed_data[fixed_data['CLASS_STAR'] < 0.95]

        return input_df

    def splits_data(self, filt_cat, input_df):
        """ Splits filtered catalogue due object type.
        This method creates a dictionary with three different keys.
        Each one (SSOs, stars, galaxies) it is a Dataframe with all valuable
        data.

        :param filt_cat:
        :param input_df:
        :return:
        """
        # out_d = {}

        # Unique sources (?)
        unique_sources = list(set(filt_cat['SOURCE_NUMBER'].tolist()))

        print('Creating new catalogues from filtered catalogue due type')
        for source_ in unique_sources:
            source_df = filt_cat[filt_cat['SOURCE_NUMBER'].isin([source_])]

            for i, row in enumerate(source_df.itertuples(), 1):
                dither_n = get_dither(int(row.CATALOG_NUMBER))
                if dither_n == 1:
                    # Checks object type
                    alpha = source_df['ALPHA_J2000'].iloc[0]
                    delta = source_df['DELTA_J2000'].iloc[0]
                    test_sso = check_source(input_df[dither_n]['SSOs'],
                                            alpha, delta)
                    test_star = check_source(input_df[dither_n]['stars'],
                                             alpha, delta)
                    test_galaxy = check_source(input_df[dither_n]['galaxies'],
                                               alpha, delta)
                    if test_sso.empty is not True:
                        key_ = 'ssos'
                    elif test_star.empty is not True:
                        key_ = 'stars'
                    elif test_galaxy.empty is not True:
                        key_ = 'galaxies'
                    else:
                        key_ = 'lost'

                    a_image = row.MEDIAN_A_IMAGE
                    self.data_d[key_]['median_a_image'].append(a_image)
                    erra_image = row.MEDIAN_ERRA_IMAGE
                    self.data_d[key_]['median_erra_image'].append(erra_image)
                    b_image = row.MEDIAN_B_IMAGE
                    self.data_d[key_]['median_b_image'].append(b_image)
                    errb_image = row.MEDIAN_ERRB_IMAGE
                    self.data_d[key_]['median_errb_image'].append(errb_image)
                    class_star = row.MEDIAN_CLASS_STAR
                    self.data_d[key_]['median_class_star'].append(class_star)
                    ellipticity = row.ELLIPTICITY
                    self.data_d[key_]['ellipticity'].append(ellipticity)
                    mag_iso = row.MEDIAN_MAG_ISO
                    self.data_d[key_]['median_mag_iso'].append(mag_iso)
                    magerr_iso = row.MEDIAN_MAGERR_ISO
                    self.data_d[key_]['median_magerr_iso'].append(magerr_iso)
                    flux_iso = row.MEDIAN_FLUX_ISO
                    self.data_d[key_]['median_flux_iso'].append(flux_iso)

        print(self.data_d)


class PlotTotalScampPerformance:

    def __init__(self):
        """
        Read data from files
        Creates a page per data
        """
        self.dpi = 100
        self.keys_to_plot = ['median_a_image', 'median_b_image',
                             'median_class_image', 'ellipticity',
                             'median_mag_iso', 'median_flux_iso']
        self.data_d = {}

        self.read_data()
        self.plot()

    def read_data(self):
        for key_ in self.keys_to_plot:
            self.data_d[key_] = read_csv('{}.csv'.format(key_), index_col=0)

    def plot(self):
        fig = plt.figure(figsize=(16.53, 11.69), dpi=self.dpi)
        ax = fig.add_subplot(1, 1, 1)

        # mag_iso vs. median_b_image
        # SSOs first

        # Stars second

        # Galaxies third
        pass


if __name__ == "__main__":
    TotalScampPerformance()
