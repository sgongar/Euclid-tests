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


def check_source(o_df, i_alpha, i_delta, keys):
    """

    :param o_df:
    :param i_alpha:
    :param i_delta:
    :param keys:
    :return:
    """
    prfs_d = extract_settings_elvis()

    o_df = o_df[o_df[keys[0]] + prfs_d['tolerance'] > i_alpha]
    o_df = o_df[i_alpha > o_df[keys[0]] - prfs_d['tolerance']]
    o_df = o_df[o_df[keys[1]] + prfs_d['tolerance'] > i_delta]
    o_df = o_df[i_delta > o_df[keys[1]] - prfs_d['tolerance']]

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

        """
        # Saves data
        for key_ in ['ssos', 'stars', 'galaxies']:
            out_df = DataFrame(self.data_d[key_])
            out_df.to_csv('{}.csv'.format(key_))
        """

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
        input_df = {1: {}, 2: {}, 3: {}, 4: {}}

        for key_ in input_df.keys():
            ssos_cat = 'cat_clean_ssos_{}.csv'.format(key_)
            input_df[key_]['SSOs'] = read_csv(ssos_cat, index_col=0)
            stars_cat = 'stars.csv'
            input_df[key_]['stars'] = read_csv(stars_cat, index_col=0)
            galaxies_cat = 'galaxies.csv'
            input_df[key_]['galaxies'] = read_csv(galaxies_cat, index_col=0)

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
        # Unique sources (?)
        unique_sources = list(set(filt_cat['SOURCE_NUMBER'].tolist()))

        print('Creating new catalogues from filtered catalogue due type')
        print('Total sources: {}'.format(filt_cat['SOURCE_NUMBER'].size))
        for source_ in unique_sources:
            source_df = filt_cat[filt_cat['SOURCE_NUMBER'].isin([source_])]

            for i, row in enumerate(source_df.itertuples(), 1):
                dither_n = get_dither(int(row.CATALOG_NUMBER))
                # Checks object type
                alpha = source_df['ALPHA_J2000'].iloc[0]
                delta = source_df['DELTA_J2000'].iloc[0]
                keys = ['ALPHA_J2000', 'DELTA_J2000']
                test_sso = check_source(input_df[dither_n]['SSOs'],
                                        alpha, delta, keys)
                keys = ['X_WORLD', 'Y_WORLD']
                test_star = check_source(input_df[dither_n]['stars'],
                                         alpha, delta, keys)

                keys = ['X_WORLD', 'Y_WORLD']
                test_galaxy = check_source(input_df[dither_n]['galaxies'],
                                           alpha, delta, keys)
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

        for type_key in self.data_d.keys():
            data_df = DataFrame(self.data_d[type_key])
            data_df.to_csv('cat_{}.csv'.format(type_key))


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

        objects = ['cat_ssos', 'cat_stars', 'cat_galaxies']
        for object_ in objects:
            self.read_data(object_)

        self.plot()

    def read_data(self, object_):
        """

        :param object_:
        :return:
        """
        self.data_d[object_] = read_csv('{}.csv'.format(object_), index_col=0)

    def plot(self):
        """

        :return:
        """
        for key_ in self.data_d.keys():
            pdf_name = '{}.pdf'.format(key_)
            data_d = self.data_d[key_]

            with PdfPages(pdf_name) as pdf:
                # MAG_ISO Galaxies
                fig_1 = plt.figure(figsize=(16.53, 11.69), dpi=self.dpi)
                ax_1 = fig_1.add_subplot(1, 1, 1)
                ax_1.set_title('MAG_ISO - Galaxies')

                ax_1.scatter(data_d['median_b_image'],
                             data_d['median_a_image'],
                             label='a_image', c='b')

                ax_1.legend(loc=4)
                ax_1.grid(True)

                pdf.savefig()  # saves figure
                plt.clf()  # clear current figure
                plt.close(fig_1)  # removes figure


if __name__ == "__main__":
    TotalScampPerformance()
    # PlotTotalScampPerformance()
