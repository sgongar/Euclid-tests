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

            # B_IMAGE / A_IMAGE Galaxies
            fig_3 = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax_3 = fig_3.add_subplot(1, 1, 1)
            ax_3.set_title('B/A_IMAGE - Galaxies')

            ax_3.scatter(b_image_g_l, a_image_g_l, label='b_image/a_image',
                         c='b')

            # ax_3.set_xlim(10, 26)
            # ax_3.set_ylim(0, 40)

            ax_3.legend(loc=4)
            ax_3.grid(True)

            pdf.savefig()  # saves figure
            plt.clf()  # clear current figure
            plt.close(fig_3)  # removes figure

            # MAG_APER
            fig_4 = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax_4 = fig_4.add_subplot(1, 1, 1)
            ax_4.set_title('MAG_APER - Galaxies')

            ax_4.scatter(mag_aper_g_l, a_image_g_l, label='a_image', c='b')
            ax_4.scatter(mag_aper_g_l, b_image_g_l, label='b_image', c='g')

            ax_4.set_xlim(10, 26)
            ax_4.set_ylim(0, 40)

            ax_4.legend(loc=4)
            ax_4.grid(True)

            pdf.savefig()  # saves figure
            plt.clf()  # clear current figure
            plt.close(fig_4)  # removes figure

            # MAG_APER
            fig_5 = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax_5 = fig_5.add_subplot(1, 1, 1)
            ax_5.set_title('MAG_APER - Stars')

            ax_5.scatter(mag_aper_s_l, a_image_s_l, label='a_image', c='b')
            ax_5.scatter(mag_aper_s_l, b_image_s_l, label='b_image', c='g')

            ax_5.set_xlim(10, 26)
            ax_5.set_ylim(0, 40)

            ax_5.legend(loc=4)
            ax_5.grid(True)

            pdf.savefig()  # saves figure
            plt.clf()  # clear current figure
            plt.close(fig_5)  # removes figure

            # B_IMAGE / A_IMAGE Stars
            fig_6 = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax_6 = fig_6.add_subplot(1, 1, 1)
            ax_6.set_title('B/A_IMAGE - Galaxies')

            ax_6.scatter(b_image_s_l, a_image_s_l, label='b_image/a_image',
                         c='b')

            # ax_3.set_xlim(10, 26)
            # ax_3.set_ylim(0, 40)

            ax_6.legend(loc=4)
            ax_6.grid(True)

            pdf.savefig()  # saves figure
            plt.clf()  # clear current figure
            plt.close(fig_6)  # removes figure


if __name__ == "__main__":
    input_catalog, full_cat = main()

    if argv[1] == '-extracted':
        extracted(input_catalog, full_cat)
    elif argv[1] == '-ab_size':
        ab_size(input_catalog, full_cat, pdf=True)

