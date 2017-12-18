#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.1

Todo:
    * Improve log messages

"""
from os import makedirs, path

from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import arange
from pandas import concat

from misc import extract_settings
from regions import Create_regions


def search(o_cat, i_alpha, i_delta):
    """

    :param o_cat:
    :param i_alpha:
    :param i_delta:
    :return:
    """
    tolerance = 0.0001

    o_df = o_cat[o_cat['alpha_j2000'] + tolerance > i_alpha]
    o_df = o_df[i_alpha > o_df['alpha_j2000'] - tolerance]
    o_df = o_df[o_df['delta_j2000'] + tolerance > i_delta]
    o_df = o_df[i_delta > o_df['delta_j2000'] - tolerance]

    return o_df


class SextractorPerformance:

    def __init__(self):
        """

        """
        self.prfs_d = extract_settings()

    def error(self, logger, mag, sex_cf):
        """

        :param logger:
        :param mag:
        :param sex_cf:
        :return:
        """

        # Input sources
        input_ssos_d = {}
        for d in range(1, 5, 1):
            cat_name = '{}/Cat_20-21_d{}'.format(self.prfs_d['input_ref'], d)
            input_ssos_d[d] = '{}.dat'.format(cat_name)
        input_ssos_d = Create_regions(input_ssos_d).check_luca(True, True)

        # Creates a DataFrame from an input dictionary
        input_ssos_l = []
        for key_ in input_ssos_d.keys():
            input_ssos_l.append(input_ssos_d[key_])

        i_ssos_df = concat(input_ssos_l, axis=0)
        i_ssos_df = i_ssos_df.reset_index(drop=True)

        # Stars
        input_stars_d = {}
        for d in range(1, 5, 1):
            cat_name = '{}/Cat_20-21_d{}'.format(self.prfs_d['input_ref'], d)
            input_stars_d[d] = '{}.dat'.format(cat_name)
        input_stars_d = Create_regions(input_stars_d).check_stars(True, True)

        # Creates a DataFrame from an input dictionary
        input_stars_l = []
        for key_ in input_stars_d.keys():
            input_stars_l.append(input_stars_d[key_])

        i_stars_df = concat(input_stars_l, axis=0)
        i_stars_df = i_stars_df.reset_index(drop=True)

        # Galaxies
        input_galaxies_d = {}
        for d in range(1, 5, 1):
            cat_name = '{}/Cat_20-21_d{}'.format(self.prfs_d['input_ref'], d)
            input_galaxies_d[d] = '{}.dat'.format(cat_name)
        input_galaxies_d = Create_regions(input_galaxies_d).check_galaxies(True,
                                                                           True)

        # Creates a DataFrame from an input dictionary
        input_galaxies_l = []
        for key_ in input_galaxies_d.keys():
            input_galaxies_l.append(input_galaxies_d[key_])

        i_galaxies_df = concat(input_galaxies_l, axis=0)
        i_galaxies_df = i_galaxies_df.reset_index(drop=True)

        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {}'.format(sex_cf))

        opts = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1],
                [1, 2], [2, 0], [2, 1], [2, 2]]
        for opt in opts:
            for dither in range(1, 5, 1):
                # Lists definition
                mags_stars = []
                errors_a_stars = []
                errors_b_stars = []
                mags_galaxies = []
                errors_a_galaxies = []
                errors_b_galaxies = []
                mags_ssos = []
                errors_a_ssos = []
                errors_b_ssos = []

                cat_n = 'mag_{}_CCD_x{}_y{}_d{}.cat'.format(mag, opt[0],
                                                            opt[1], dither)
                cat_o_n = '{}/{}/{}'.format(self.prfs_d['fits_dir'], sex_cf,
                                            cat_n)

                hdu_list = fits.open(cat_o_n)
                o_cat = Table(hdu_list[2].data).to_pandas()

                i_stars_d_df = i_stars_df[
                    i_stars_df['dither_values'].isin([dither])]
                i_galaxies_d_df = i_galaxies_df[
                    i_galaxies_df['dither_values'].isin([dither])]
                i_ssos_d_df = i_ssos_df[
                    i_ssos_df['dither_values'].isin([dither])]
                for idx, row in o_cat.iterrows():
                    tmp_count = 0

                    i_alpha = row['ALPHA_J2000']
                    i_delta = row['DELTA_J2000']
                    o_df = search(i_stars_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        mags_stars.append(row['MAG_AUTO'])
                        errors_a_stars.append(row['ERRA_IMAGE'])
                        errors_b_stars.append(row['ERRB_IMAGE'])
                        tmp_count += 1
                    o_df = search(i_galaxies_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        mags_galaxies.append(row['MAG_AUTO'])
                        errors_a_galaxies.append(row['ERRA_IMAGE'])
                        errors_b_galaxies.append(row['ERRB_IMAGE'])
                        tmp_count += 1
                    o_df = search(i_ssos_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        mags_ssos.append(row['MAG_AUTO'])
                        errors_a_ssos.append(row['ERRA_IMAGE'])
                        errors_b_ssos.append(row['ERRB_IMAGE'])
                        tmp_count += 1

                sextractor_folder = '{}/{}'.format(self.prfs_d['plots_dir'],
                                                   sex_cf)
                if not path.exists(sextractor_folder):
                    makedirs(sextractor_folder)

                with PdfPages('{}/{}_log.pdf'.format(sextractor_folder,
                                                     cat_n[:-4])) as pdf:
                    fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                    ax = fig.add_subplot(1, 1, 1)

                    ax.semilogy(mags_stars, errors_a_stars, 'bs', ms=1)
                    ax.semilogy(mags_galaxies, errors_a_galaxies, 'gs', ms=1)
                    ax.semilogy(mags_ssos, errors_a_ssos, 'rs', ms=1)

                    ax.set_xlabel('mag')
                    ax.set_ylabel('error a')
                    ax.set_xlim([14, 28])

                    # x-scale
                    x_major_ticks = arange(14, 28, 0.5)
                    x_minor_ticks = arange(14, 28, 0.1)
                    ax.set_xticks(x_major_ticks, minor=False)
                    ax.set_xticks(x_minor_ticks, minor=True)

                    ax.grid(b=True, which='major', ls='-', lw=2)
                    ax.grid(b=True, which='minor', ls='--', lw=1)
                    pdf.savefig()  # saves current figure
                    plt.clf()  # clear current figure

                    # B parameters
                    fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                    ax = fig.add_subplot(1, 1, 1)

                    ax.semilogy(mags_stars, errors_b_stars, 'bs', ms=1)
                    ax.semilogy(mags_galaxies, errors_b_galaxies, 'gs', ms=1)
                    ax.semilogy(mags_ssos, errors_b_ssos, 'rs', ms=1)

                    ax.set_xlabel('mag')
                    ax.set_ylabel('error b')
                    ax.set_xlim([14, 28])

                    # x-scale
                    x_major_ticks = arange(14, 28, 0.5)
                    x_minor_ticks = arange(14, 28, 0.1)
                    ax.set_xticks(x_major_ticks, minor=False)
                    ax.set_xticks(x_minor_ticks, minor=True)

                    ax.grid(b=True, which='major', ls='-', lw=2)
                    ax.grid(b=True, which='minor', ls='--', lw=1)

                    pdf.savefig()  # saves current figure
                    plt.clf()  # clear current figure
