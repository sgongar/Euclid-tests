#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

La idea es que compare los objetos extraidos en total con los que tengo.

Todo:
    * Improve log messages

"""

from sys import argv

from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pandas import concat, read_csv

from misc import extract_settings_sc3

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"

"""
'2-000000-0000000__20170630T011642.0Z_00.00_0',
'2-000000-0000000__20170630T011642.0Z_00.00_10',
'2-000000-0000000__20170630T011642.0Z_00.00_11',
'2-000000-0000000__20170630T011642.0Z_00.00_12',
'2-000000-0000000__20170630T011642.0Z_00.00_13',
'2-000000-0000000__20170630T011642.0Z_00.00_14',
'2-000000-0000000__20170630T011642.0Z_00.00_15',
'2-000000-0000000__20170630T011642.0Z_00.00_16',
'2-000000-0000000__20170630T011642.0Z_00.00_17',
'2-000000-0000000__20170630T011642.0Z_00.00_18',
'2-000000-0000000__20170630T011642.0Z_00.00_19',
'2-000000-0000000__20170630T011642.0Z_00.00_1',
'2-000000-0000000__20170630T011642.0Z_00.00_20',
'2-000000-0000000__20170630T011642.0Z_00.00_21',
'2-000000-0000000__20170630T011642.0Z_00.00_22',
'2-000000-0000000__20170630T011642.0Z_00.00_23',
'2-000000-0000000__20170630T011642.0Z_00.00_24',
'2-000000-0000000__20170630T011642.0Z_00.00_25',
'2-000000-0000000__20170630T011642.0Z_00.00_26',
'2-000000-0000000__20170630T011642.0Z_00.00_27',
'2-000000-0000000__20170630T011642.0Z_00.00_28',
'2-000000-0000000__20170630T011642.0Z_00.00_29',
'2-000000-0000000__20170630T011642.0Z_00.00_2',
'2-000000-0000000__20170630T011642.0Z_00.00_30',
'2-000000-0000000__20170630T011642.0Z_00.00_31',
'2-000000-0000000__20170630T011642.0Z_00.00_32',
'2-000000-0000000__20170630T011642.0Z_00.00_33',
'2-000000-0000000__20170630T011642.0Z_00.00_34',
'2-000000-0000000__20170630T011642.0Z_00.00_35',
'2-000000-0000000__20170630T011642.0Z_00.00_3',
'2-000000-0000000__20170630T011642.0Z_00.00_4',
'2-000000-0000000__20170630T011642.0Z_00.00_5',
'2-000000-0000000__20170630T011642.0Z_00.00_6',
'2-000000-0000000__20170630T011642.0Z_00.00_7',
'2-000000-0000000__20170630T011642.0Z_00.00_8',
'2-000000-0000000__20170630T011642.0Z_00.00_9',
'3-000000-0000000__20170630T011848.6Z_00.00_0',
'3-000000-0000000__20170630T011848.6Z_00.00_10',
'3-000000-0000000__20170630T011848.6Z_00.00_11',
'3-000000-0000000__20170630T011848.6Z_00.00_12',
'3-000000-0000000__20170630T011848.6Z_00.00_13',
'3-000000-0000000__20170630T011848.6Z_00.00_14',
'3-000000-0000000__20170630T011848.6Z_00.00_15',
'3-000000-0000000__20170630T011848.6Z_00.00_16',
'3-000000-0000000__20170630T011848.6Z_00.00_17',
'3-000000-0000000__20170630T011848.6Z_00.00_18',
'3-000000-0000000__20170630T011848.6Z_00.00_19',
'3-000000-0000000__20170630T011848.6Z_00.00_1',
'3-000000-0000000__20170630T011848.6Z_00.00_20',
'3-000000-0000000__20170630T011848.6Z_00.00_21',
'3-000000-0000000__20170630T011848.6Z_00.00_22',
'3-000000-0000000__20170630T011848.6Z_00.00_23',
'3-000000-0000000__20170630T011848.6Z_00.00_24',
'3-000000-0000000__20170630T011848.6Z_00.00_25',
'3-000000-0000000__20170630T011848.6Z_00.00_26',
'3-000000-0000000__20170630T011848.6Z_00.00_27',
'3-000000-0000000__20170630T011848.6Z_00.00_28',
'3-000000-0000000__20170630T011848.6Z_00.00_29',
'3-000000-0000000__20170630T011848.6Z_00.00_2',
'3-000000-0000000__20170630T011848.6Z_00.00_30',
'3-000000-0000000__20170630T011848.6Z_00.00_31',
'3-000000-0000000__20170630T011848.6Z_00.00_32',
'3-000000-0000000__20170630T011848.6Z_00.00_33',
'3-000000-0000000__20170630T011848.6Z_00.00_34',
'3-000000-0000000__20170630T011848.6Z_00.00_35',
'3-000000-0000000__20170630T011848.6Z_00.00_3',
'3-000000-0000000__20170630T011848.6Z_00.00_4',
'3-000000-0000000__20170630T011848.6Z_00.00_5',
'3-000000-0000000__20170630T011848.6Z_00.00_6',
'3-000000-0000000__20170630T011848.6Z_00.00_7',
'3-000000-0000000__20170630T011848.6Z_00.00_8',
'3-000000-0000000__20170630T011848.6Z_00.00_9',
'4-000000-0000000__20170630T012050.1Z_00.00_0',
'4-000000-0000000__20170630T012050.1Z_00.00_10',
'4-000000-0000000__20170630T012050.1Z_00.00_11',
'4-000000-0000000__20170630T012050.1Z_00.00_12',
'4-000000-0000000__20170630T012050.1Z_00.00_13',
'4-000000-0000000__20170630T012050.1Z_00.00_14',
'4-000000-0000000__20170630T012050.1Z_00.00_15',
'4-000000-0000000__20170630T012050.1Z_00.00_16',
'4-000000-0000000__20170630T012050.1Z_00.00_17',
'4-000000-0000000__20170630T012050.1Z_00.00_18',
'4-000000-0000000__20170630T012050.1Z_00.00_19',
'4-000000-0000000__20170630T012050.1Z_00.00_1',
'4-000000-0000000__20170630T012050.1Z_00.00_20',
'4-000000-0000000__20170630T012050.1Z_00.00_21',
'4-000000-0000000__20170630T012050.1Z_00.00_22',
'4-000000-0000000__20170630T012050.1Z_00.00_23',
'4-000000-0000000__20170630T012050.1Z_00.00_24',
'4-000000-0000000__20170630T012050.1Z_00.00_25',
'4-000000-0000000__20170630T012050.1Z_00.00_26',
'4-000000-0000000__20170630T012050.1Z_00.00_27',
'4-000000-0000000__20170630T012050.1Z_00.00_28',
'4-000000-0000000__20170630T012050.1Z_00.00_29',
'4-000000-0000000__20170630T012050.1Z_00.00_2',
'4-000000-0000000__20170630T012050.1Z_00.00_30',
'4-000000-0000000__20170630T012050.1Z_00.00_31',
'4-000000-0000000__20170630T012050.1Z_00.00_32',
'4-000000-0000000__20170630T012050.1Z_00.00_33',
'4-000000-0000000__20170630T012050.1Z_00.00_34',
'4-000000-0000000__20170630T012050.1Z_00.00_35',
'4-000000-0000000__20170630T012050.1Z_00.00_3',
'4-000000-0000000__20170630T012050.1Z_00.00_4',
'4-000000-0000000__20170630T012050.1Z_00.00_5',
'4-000000-0000000__20170630T012050.1Z_00.00_6',
'4-000000-0000000__20170630T012050.1Z_00.00_7',
'4-000000-0000000__20170630T012050.1Z_00.00_8',
'4-000000-0000000__20170630T012050.1Z_00.00_9']
"""


def get_cat(cat_n):
    """  TODO get scamp organization in order to reduce the lines number

    :param cat_n:
    :return: cat_file
    """
    cat_part_1 = 'EUC_VIS_SWL-DET-001-000000-0000000' \
                 '__20170630T011437.3Z_00.00_'

    cats_final = []
    cat_indexes = ['0', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                   '19', '1', '20', '21', '22', '23', '24', '25', '26', '27',
                   '28', '29', '2', '30', '31', '32', '33', '34', '35', '3',
                   '4', '5', '6', '7', '8', '9']

    cats_final.append(' ')
    for cat_ in cat_indexes:
        cats_final.append('{}{}.cat'.format(cat_part_1, cat_))

    cat_file = cats_final[cat_n]

    return cat_file


def check_source(i_alpha, i_delta, e_df):
    """

    :param i_alpha:
    :param i_delta:
    :param e_df:
    :return: e_df
    """
    prfs_d = extract_settings_sc3()

    e_df = e_df[e_df['ALPHA_J2000'] + prfs_d['tolerance'] > i_alpha]
    e_df = e_df[i_alpha > e_df['ALPHA_J2000'] - prfs_d['tolerance']]
    e_df = e_df[e_df['DELTA_J2000'] + prfs_d['tolerance'] > i_delta]
    e_df = e_df[i_delta > e_df['DELTA_J2000'] - prfs_d['tolerance']]

    return e_df


def load_sextractor_cats():
    """

    :return: cat_d
    """
    prfs_d = extract_settings_sc3()

    cats_number = 36  # todo hardcoded!
    cat_d = {}
    for cat_n in range(1, cats_number + 1, 1):
        cat_file = get_cat(cat_n)
        cat_data = fits.open('{}/{}'.format(prfs_d['fits_dir'], cat_file))

        ccd_df = Table(cat_data[2].data)
        # print('CCD catalog {} to Pandas'.format(cat_n))
        cat_d[cat_n] = ccd_df.to_pandas()

    return cat_d


def merge_cats(cat_d):
    """

    :param cat_d:
    :return:
    """
    cat_list = []

    for idx, cat_ in enumerate(cat_d.keys()):
        cat_list.append(cat_d[cat_])

    full_cat = concat(cat_list)
    full_cat = full_cat.reset_index(drop=True)

    return full_cat


def main():
    """

    :return:
    """
    input_catalog = read_csv('1523634221737L.csv')
    input_catalog = input_catalog[input_catalog['dither'].isin([1])]
    # Load all catalogs
    cat_d = load_sextractor_cats()
    # Get boundaries for all catalogs
    full_cat = merge_cats(cat_d)

    return input_catalog, full_cat


def extracted(input_catalog, full_cat):
    """

    :param input_catalog:
    :param full_cat:
    :return:
    """
    # Look for sources
    total_ones = input_catalog['rightascension'].size
    right_ones = 0
    false_ones = 0

    for i, row in enumerate(input_catalog.itertuples(), 1):
        # print('source number: {} - total number: {}'.format(i, total_ones))
        i_alpha = row.rightascension
        i_delta = row.declination

        e_df = check_source(i_alpha, i_delta, full_cat)

        if e_df.empty is not True:
            right_ones += 1
        else:
            false_ones += 1

    print('total ones {}'.format(total_ones))
    print('right ones {}'.format(right_ones))
    print('false ones {}'.format(false_ones))


def ab_size(input_catalog, full_cat, pdf):
    """

    :param input_catalog:
    :param full_cat:
    :return:
    """
    # Look for sources
    total_ones = input_catalog['rightascension'].size
    a_image_g_l = []
    b_image_g_l = []
    mag_iso_g_l = []
    mag_aper_g_l = []
    a_image_s_l = []
    b_image_s_l = []
    mag_iso_s_l = []
    mag_aper_s_l = []

    for i, row in enumerate(input_catalog.itertuples(), 1):
        print('source number: {} - total number: {}'.format(i, total_ones))
        i_alpha = row.rightascension
        i_delta = row.declination
        i_class = str(row.starflag)

        e_df = check_source(i_alpha, i_delta, full_cat)

        if e_df.empty is not True:
            if i_class == '0':
                a_image_g_l.append(e_df['A_IMAGE'].iloc[0])
                b_image_g_l.append(e_df['B_IMAGE'].iloc[0])
                mag_iso_g_l.append(e_df['MAG_ISO'].iloc[0])
                mag_aper_g_l.append(e_df['MAG_APER'].iloc[0])
            elif i_class == '1':
                a_image_s_l.append(e_df['A_IMAGE'].iloc[0])
                b_image_s_l.append(e_df['B_IMAGE'].iloc[0])
                mag_iso_s_l.append(e_df['MAG_ISO'].iloc[0])
                mag_aper_s_l.append(e_df['MAG_APER'].iloc[0])
        else:
            pass

    if pdf:
        # PDF parameters
        plot_size = [16.53, 11.69]
        plot_dpi = 100
        pdf_name = 'movement.pdf'

        with PdfPages(pdf_name) as pdf:
            # MAG_ISO Galaxies
            fig_1 = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax_1 = fig_1.add_subplot(1, 1, 1)
            ax_1.set_title('MAG_ISO - Galaxies')

            ax_1.scatter(mag_iso_g_l, a_image_g_l, label='a_image', c='b')
            ax_1.scatter(mag_iso_g_l, b_image_g_l, label='b_image', c='g')

            ax_1.set_xlim(10, 26)
            ax_1.set_ylim(0, 40)

            ax_1.legend(loc=4)
            ax_1.grid(True)

            pdf.savefig()  # saves figure
            plt.clf()  # clear current figure
            plt.close(fig_1)  # removes figure

            # MAG_ISO Stars
            fig_2 = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax_2 = fig_2.add_subplot(1, 1, 1)
            ax_2.set_title('MAG_ISO - Stars')

            ax_2.scatter(mag_iso_s_l, a_image_s_l, label='a_image', c='b')
            ax_2.scatter(mag_iso_s_l, b_image_s_l, label='b_image', c='g')

            ax_2.set_xlim(10, 26)
            ax_2.set_ylim(0, 40)

            ax_2.legend(loc=4)
            ax_2.grid(True)

            pdf.savefig()  # saves figure
            plt.clf()  # clear current figure
            plt.close(fig_2)  # removes figure

            # MAG_APER
            fig_3 = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax_3 = fig_3.add_subplot(1, 1, 1)
            ax_3.set_title('MAG_APER - Galaxies')

            ax_3.scatter(mag_aper_g_l, a_image_g_l, label='a_image', c='b')
            ax_3.scatter(mag_aper_g_l, b_image_g_l, label='b_image', c='g')

            ax_3.set_xlim(10, 26)
            ax_3.set_ylim(0, 40)

            ax_3.legend(loc=4)
            ax_3.grid(True)

            pdf.savefig()  # saves figure
            plt.clf()  # clear current figure
            plt.close(fig_3)  # removes figure

            # MAG_APER
            fig_4 = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax_4 = fig_4.add_subplot(1, 1, 1)
            ax_4.set_title('TEST')

            ax_4.scatter(mag_aper_s_l, a_image_s_l, label='a_image', c='b')
            ax_4.scatter(mag_aper_s_l, b_image_s_l, label='b_image', c='g')

            ax_4.set_xlim(10, 26)
            ax_4.set_ylim(0, 40)

            ax_4.legend(loc=4)
            ax_4.grid(True)

            pdf.savefig()  # saves figure
            plt.clf()  # clear current figure
            plt.close(fig_3)  # removes figure


if __name__ == "__main__":
    input_catalog, full_cat = main()

    if argv[1] == '-extracted':
        extracted(input_catalog, full_cat)
    elif argv[1] == '-ab_size':
        ab_size(input_catalog, full_cat, pdf=True)

