#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.2 Now supports confidence intervals
- 0.3 Sextractor performance's added
- 0.4 Proper motion performance's added
Todo:
    * Improve log messages

"""
from os import makedirs, path

from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import arange, array
import numpy as np
from pandas import concat, read_csv, DataFrame

from misc import all_same, speeds_range
from regions import Create_regions


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.4"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class StatsPerformance:
    def __init__(self):
        """

        """
        pass

    def error(self, logger, prfs_d, mag, sex_cf):
        """

        :param logger:
        :param prfs_d:
        :param mag:
        :param sex_cf:
        :return:
        """
        mags_stars = []
        errors_a_stars = []
        errors_b_stars = []
        mags_galaxies = []
        errors_a_galaxies = []
        errors_b_galaxies = []
        mags_ssos = []
        errors_a_ssos = []
        errors_b_ssos = []

        # Input sources
        input_ssos_d = {}
        for d in range(1, 5, 1):
            input_ssos_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                            d)
        input_ssos_d = Create_regions(input_ssos_d, prfs_d).check_luca(True, True)

        # Creates a DataFrame from an input dictionary
        input_ssos_l = []
        for key_ in input_ssos_d.keys():
            input_ssos_l.append(input_ssos_d[key_])

        i_ssos_df = concat(input_ssos_l, axis=0)
        i_ssos_df = i_ssos_df.reset_index(drop=True)

        # Stars
        input_stars_d = {}
        for d in range(1, 5, 1):
            input_stars_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                             d)
        input_stars_d = Create_regions(input_stars_d, prfs_d).check_stars(True, True)

        # Creates a DataFrame from an input dictionary
        input_stars_l = []
        for key_ in input_stars_d.keys():
            input_stars_l.append(input_stars_d[key_])

        i_stars_df = concat(input_stars_l, axis=0)
        i_stars_df = i_stars_df.reset_index(drop=True)

        # Galaxies
        input_galaxies_d = {}
        for d in range(1, 5, 1):
            input_galaxies_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                                d)
        input_galaxies_d = Create_regions(input_galaxies_d,
                                          prfs_d).check_galaxies(True, True)

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
                errors_a_stars = []
                errors_b_stars = []
                errors_a_galaxies = []
                errors_b_galaxies = []
                errors_a_ssos = []
                errors_b_ssos = []

                cat_n = 'mag_{}_CCD_x{}_y{}_d{}.cat'.format(mag, opt[0], opt[1],
                                                            dither)
                cat_o_n = '{}/{}/{}'.format(prfs_d['fits_dir'], sex_cf, cat_n)

                hdu_list = fits.open(cat_o_n)
                o_cat = Table(hdu_list[2].data).to_pandas()

                i_stars_d_df = i_stars_df[i_stars_df['dither_values'].isin([dither])]
                i_galaxies_d_df = i_galaxies_df[i_galaxies_df['dither_values'].isin([dither])]
                i_ssos_d_df = i_ssos_df[i_ssos_df['dither_values'].isin([dither])]
                for idx, row in o_cat.iterrows():
                    i_alpha = row['ALPHA_J2000']
                    i_delta = row['DELTA_J2000']
                    o_df = self.search(i_stars_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        errors_a_stars.append(row['ERRA_IMAGE'])
                        errors_b_stars.append(row['ERRB_IMAGE'])
                    o_df = self.search(i_galaxies_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        errors_a_galaxies.append(row['ERRA_IMAGE'])
                        errors_b_galaxies.append(row['ERRB_IMAGE'])
                    o_df = self.search(i_ssos_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        errors_a_ssos.append(row['ERRA_IMAGE'])
                        errors_b_ssos.append(row['ERRB_IMAGE'])

        """
        mean_a_stars = array(errors_a_stars).mean()
        mean_b_stars = array(errors_b_stars).mean()

        mean_a_galaxies = array(errors_a_galaxies).mean()
        mean_b_galaxies = array(errors_b_galaxies).mean()

        mean_a_ssos = array(errors_a_ssos).mean()
        mean_b_ssos = array(errors_b_ssos).mean()

        stats_d = {'m_a_stars': mean_a_stars, 'm_b_stars': mean_b_stars,
                   'm_a_gals': mean_a_galaxies, 'm_b_gals': mean_b_galaxies,
                   'm_a_ssos': mean_a_ssos, 'm_b_ssos': mean_b_ssos}
        """

        std_a_stars = array(errors_a_stars).std()
        std_b_stars = array(errors_b_stars).std()

        std_a_galaxies = array(errors_a_galaxies).std()
        std_b_galaxies = array(errors_b_galaxies).std()

        std_a_ssos = array(errors_a_ssos).std()
        std_b_ssos = array(errors_b_ssos).std()

        stats_d = {'conf': sex_cf, 's_a_stars': std_a_stars, 's_b_stars': std_b_stars,
                   's_a_gals': std_a_galaxies, 's_b_gals': std_b_galaxies,
                   's_a_ssos': std_a_ssos, 's_b_ssos': std_b_ssos}



        return stats_d

    def search(self, o_cat, i_alpha, i_delta):
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


class PMPerformance:

    def __init__(self):
        """

        """
        pass

    def check(self, logger, prfs_d, mag, scmp_cf, sex_cf, confidence_):
        """

        :param logger:
        :param prfs_d:
        :param mag:
        :param scmp_cf:
        :param sex_cf:
        :param confidence_:
        :return:
        """
        input_pm = []
        output_pm = []

        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {} and {}'.format(scmp_cf,
                                                                 sex_cf))
        input_d = {}
        for d in range(1, 5, 1):
            input_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                       d)
        input_d = Create_regions(input_d, prfs_d).check_luca(True, True)

        # Creates a DataFrame from an input dictionary
        input_l = []
        for key_ in input_d.keys():
            input_l.append(input_d[key_])

        i_df = concat(input_l, axis=0)
        # Look for < 3 coincidences
        i_df = concat(g for _, g in i_df.groupby('source')
                      if len(g) >= 3)
        i_df = i_df.reset_index(drop=True)

        # Open particular file!
        filt_n = 'filt_{}_{}_4.csv'.format(scmp_cf, mag)
        filter_o_n = '{}/{}/{}/{}'.format(prfs_d['filter_dir'],
                                          sex_cf, scmp_cf, filt_n)

        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        # Gets unique sources from input data
        unique_sources = list(set(i_df['source'].tolist()))

        # Loops over input data
        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            cat = i_df[i_df['source'].isin([source_])]
            # cat has, at least, three or more rows, one for each dither

            # Creates lists for each source
            boolean_l = []
            tmp_catalog = []
            tmp_source = []
            # tmp lists
            tmp_input_pm = []
            tmp_output_pm = []

            # Iterate over each detection of each source
            for i, row in enumerate(cat.itertuples(), 1):
                catalog_n = row.catalog
                i_pm = row.pm_values
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # Checks if there is a source closed to input one
                # o_cat contains data from output (filtered) catalog
                o_df = self.check_source(catalog_n, o_cat,
                                         i_alpha, i_delta)

                # If there is one saves data from input data
                if o_df.empty is not True and o_df['PM'].size == 1:
                    pm_mask = self.pm_filter(o_df, i_pm, prfs_d, confidence_)
                    if pm_mask:
                        if o_df['SOURCE_NUMBER'].size != 1:
                            boolean_l.append('False')
                        else:
                            boolean_l.append('True')
                        tmp_catalog.append(catalog_n)
                        tmp_source.append(o_df['SOURCE_NUMBER'].iloc[0])
                        tmp_input_pm.append(i_pm)
                        o_pm = o_df['PM'].iloc[0]
                        tmp_output_pm.append(o_pm)
                else:
                    boolean_l.append('False')

            # check
            if len(tmp_input_pm) >= 3:
                input_pm.append(i_pm)
                output_pm.append(o_pm)

        if not self.plot_comp(prfs_d, sex_cf, scmp_cf, input_pm, output_pm, confidence_):
            raise Exception

        stats_d = {}
        return stats_d

    def pm_filter(self, o_df, pm, prfs_d, confidence_):
        """

        :param o_df:
        :param pm:
        :param prfs_d:
        :param confidence_:
        :return:
        """
        pm_ranges = speeds_range(prfs_d, confidence_)
        pm_range = pm_ranges[pm]

        if pm_range[0] < float(o_df['PM']) < pm_range[1]:
            return True
        else:
            return False

    def check_source(self, catalog_n, o_cat, i_alpha, i_delta):
        """

        :param catalog_n:
        :param o_cat:
        :param i_alpha:
        :param i_delta:
        :return:
        """
        tolerance = 0.001

        o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
        o_df = o_df[o_df['ALPHA_J2000'] + tolerance > i_alpha]
        o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - tolerance]
        o_df = o_df[o_df['DELTA_J2000'] + tolerance > i_delta]
        o_df = o_df[i_delta > o_df['DELTA_J2000'] - tolerance]

        return o_df

    def plot_comp(self, prfs_d, sex_cf, scmp_cf,
                  input_pm, output_pm, confidence_):
        """

        :param prfs_d:
        :param sex_cf:
        :param scmp_cf:
        :param input_pm:
        :param output_pm:
        :param confidence_:
        :return:
        """
        sextractor_folder = '{}/pm/{}'.format(prfs_d['plots_dir'], sex_cf)
        if not path.exists(sextractor_folder):
            makedirs(sextractor_folder)

        with PdfPages('{}/pm/{}/{}_{}_cmp.pdf'.format(prfs_d['plots_dir'], sex_cf,
                                                      scmp_cf, confidence_)) as pdf:
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)

            input_pm_cutted = filter(lambda a: a < 4, input_pm)
            output_pm_cutted = filter(lambda a: a < 4, output_pm)
            ax.plot(input_pm_cutted, output_pm_cutted, 'bs', markersize=1)
            ax.plot([0, 4], [0, 4])

            ax.set_xlabel('pm input [0 - 4]')
            ax.set_ylabel('pm output [0 - 4]')

            # x-scale
            x_major_ticks = arange(0, 4, 0.5)
            x_minor_ticks = arange(0, 4, 0.1)
            ax.set_xticks(x_major_ticks, minor=False)
            ax.set_xticks(x_minor_ticks, minor=True)

            y_major_ticks = arange(0, 4, 0.5)
            y_minor_ticks = arange(0, 4, 0.1)
            ax.set_yticks(y_major_ticks, minor=False)
            ax.set_yticks(y_minor_ticks, minor=True)

            ax.grid(b=True, which='major', linestyle='-', linewidth=2)
            ax.grid(b=True, which='minor', linestyle='--', linewidth=1)
            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure

            # [8-10]
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)

            input_pm_cutted = filter(lambda a: a > 4, input_pm)
            output_pm_cutted = filter(lambda a: a > 4, output_pm)
            ax.plot(input_pm_cutted, output_pm_cutted, 'bs', markersize=1)
            ax.plot([9.5, 10.5], [9.5, 10.5])

            ax.set_xlabel('pm input [9.5 - 10.5]')
            ax.set_ylabel('pm output [9.5 - 10.5]')

            # x-scale
            x_major_ticks = arange(9.5, 10.5, 0.5)
            x_minor_ticks = arange(9.5, 10.5, 0.05)
            ax.set_xticks(x_major_ticks, minor=False)
            ax.set_xticks(x_minor_ticks, minor=True)

            y_major_ticks = arange(9.5, 10.5, 0.5)
            y_minor_ticks = arange(9.5, 10.5, 0.05)
            ax.set_yticks(y_major_ticks, minor=False)
            ax.set_yticks(y_minor_ticks, minor=True)

            ax.grid(b=True, which='major', linestyle='-', linewidth=2)
            ax.grid(b=True, which='minor', linestyle='--', linewidth=1)
            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure


        return True


class SextractorPerformance:

    def __init__(self):
        """

        """
        pass

    def error(self, logger, prfs_d, mag, sex_cf):
        """

        :param logger:
        :param prfs_d:
        :param mag:
        :param sex_cf:
        :return:
        """

        # Input sources
        input_ssos_d = {}
        for d in range(1, 5, 1):
            input_ssos_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                            d)
        input_ssos_d = Create_regions(input_ssos_d, prfs_d).check_luca(True, True)

        # Creates a DataFrame from an input dictionary
        input_ssos_l = []
        for key_ in input_ssos_d.keys():
            input_ssos_l.append(input_ssos_d[key_])

        i_ssos_df = concat(input_ssos_l, axis=0)
        i_ssos_df = i_ssos_df.reset_index(drop=True)

        # Stars
        input_stars_d = {}
        for d in range(1, 5, 1):
            input_stars_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                             d)
        input_stars_d = Create_regions(input_stars_d, prfs_d).check_stars(True, True)

        # Creates a DataFrame from an input dictionary
        input_stars_l = []
        for key_ in input_stars_d.keys():
            input_stars_l.append(input_stars_d[key_])

        i_stars_df = concat(input_stars_l, axis=0)
        i_stars_df = i_stars_df.reset_index(drop=True)

        # Galaxies
        input_galaxies_d = {}
        for d in range(1, 5, 1):
            input_galaxies_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                                d)
        input_galaxies_d = Create_regions(input_galaxies_d,
                                          prfs_d).check_galaxies(True, True)

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

                cat_n = 'mag_{}_CCD_x{}_y{}_d{}.cat'.format(mag, opt[0], opt[1],
                                                            dither)
                cat_o_n = '{}/{}/{}'.format(prfs_d['fits_dir'], sex_cf, cat_n)

                hdu_list = fits.open(cat_o_n)
                o_cat = Table(hdu_list[2].data).to_pandas()

                i_stars_d_df = i_stars_df[i_stars_df['dither_values'].isin([dither])]
                i_galaxies_d_df = i_galaxies_df[i_galaxies_df['dither_values'].isin([dither])]
                i_ssos_d_df = i_ssos_df[i_ssos_df['dither_values'].isin([dither])]
                for idx, row in o_cat.iterrows():
                    tmp_count = 0

                    i_alpha = row['ALPHA_J2000']
                    i_delta = row['DELTA_J2000']
                    o_df = self.search(i_stars_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        mags_stars.append(row['MAG_AUTO'])
                        errors_a_stars.append(row['ERRA_IMAGE'])
                        errors_b_stars.append(row['ERRB_IMAGE'])
                        tmp_count += 1
                    o_df = self.search(i_galaxies_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        mags_galaxies.append(row['MAG_AUTO'])
                        errors_a_galaxies.append(row['ERRA_IMAGE'])
                        errors_b_galaxies.append(row['ERRB_IMAGE'])
                        tmp_count += 1
                    o_df = self.search(i_ssos_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        mags_ssos.append(row['MAG_AUTO'])
                        errors_a_ssos.append(row['ERRA_IMAGE'])
                        errors_b_ssos.append(row['ERRB_IMAGE'])
                        tmp_count += 1

                sextractor_folder = '{}/{}'.format(prfs_d['plots_dir'], sex_cf)
                if not path.exists(sextractor_folder):
                    makedirs(sextractor_folder)

                with PdfPages('{}/{}_log.pdf'.format(sextractor_folder,
                                                     cat_n[:-4])) as pdf:
                    fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                    ax = fig.add_subplot(1, 1, 1)

                    ax.semilogy(mags_stars, errors_a_stars, 'bs', markersize=1)
                    ax.semilogy(mags_galaxies, errors_a_galaxies, 'gs', markersize=1)
                    ax.semilogy(mags_ssos, errors_a_ssos, 'rs', markersize=1)

                    ax.set_xlabel('mag')
                    ax.set_ylabel('error a')
                    ax.set_xlim([14, 28])
            
                    # x-scale
                    x_major_ticks = arange(14, 28, 0.5)
                    x_minor_ticks = arange(14, 28, 0.1)
                    ax.set_xticks(x_major_ticks, minor=False)
                    ax.set_xticks(x_minor_ticks, minor=True)

                    ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                    ax.grid(b=True, which='minor', linestyle='--', linewidth=1)
                    pdf.savefig()  # saves current figure
                    plt.clf()  # clear current figure

                    # B parameters
                    fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                    ax = fig.add_subplot(1, 1, 1)

                    ax.semilogy(mags_stars, errors_b_stars, 'bs', markersize=1)
                    ax.semilogy(mags_galaxies, errors_b_galaxies, 'gs', markersize=1)
                    ax.semilogy(mags_ssos, errors_b_ssos, 'rs', markersize=1)

                    ax.set_xlabel('mag')
                    ax.set_ylabel('error b')
                    ax.set_xlim([14, 28])

                    # x-scale
                    x_major_ticks = arange(14, 28, 0.5)
                    x_minor_ticks = arange(14, 28, 0.1)
                    ax.set_xticks(x_major_ticks, minor=False)
                    ax.set_xticks(x_minor_ticks, minor=True)

                    ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                    ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                    pdf.savefig()  # saves current figure
                    plt.clf()  # clear current figure

    def search(self, o_cat, i_alpha, i_delta):
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


class ScampPerformance:

    def __init__(self):
        """

        """
        pass

    def check(self, logger, prfs_d, mag, scmp_cf, sex_cf,
              idx_file, confidence_):
        """

        :param logger:
        :param prfs_d:
        :param mag:
        :param scmp_cf:
        :param sex_cf:
        :param idx_file:
        :param confidence_:
        :return:
        """
        # For now any file is saved
        save = False

        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {} and {}'.format(scmp_cf,
                                                                 sex_cf))
        input_d = {}
        for d in range(1, 5, 1):
            input_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                       d)
        input_d = Create_regions(input_d, prfs_d).check_luca(True, True)

        # Creates a DataFrame from an input dictionary
        input_l = []
        for key_ in input_d.keys():
            input_l.append(input_d[key_])

        i_df = concat(input_l, axis=0)
        # Look for < 3 coincidences
        i_df = concat(g for _, g in i_df.groupby('source')
                      if len(g) >= 3)
        i_df = i_df.reset_index(drop=True)

        # if save is true saves a csv file populated by input sources
        # and a regions file of them.
        if save:
            i_df.to_csv('input_sources.csv')

            alpha_df = i_df['alpha_j2000']
            delta_df = i_df['delta_j2000']

            df = concat([alpha_df, delta_df], axis=1)
            df.to_csv('input_sources.reg')

        # Open particular file!
        filt_n = 'filt_{}_{}_3.csv'.format(scmp_cf, mag)
        filter_o_n = '{}/{}/{}/{}'.format(prfs_d['filter_dir'],
                                          sex_cf, scmp_cf, filt_n)

        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)
        stats_d, out_d = self.create_dict(scmp_cf, sex_cf, confidence_)

        # Gets unique sources from input data
        unique_sources = list(set(i_df['source'].tolist()))

        # Loops over input data
        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            cat = i_df[i_df['source'].isin([source_])]
            # Creates lists for each source
            boolean_l = []
            tmp_catalog = []
            tmp_source = []
            tmp_pm = []
            tmp_alpha = []
            tmp_delta = []
            # Creates a flag for right detections
            # Initial value will set to False
            flag_detection = False

            # Iterate over each detection of each source
            for i, row in enumerate(cat.itertuples(), 1):
                source_ = row.source
                ccd_ = row.CCD
                dither_ = row.dither_values
                catalog_n = row.catalog
                pm = row.pm_values
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # Checks if there is a source closed to input one
                # o_cat contains data from output (filtered) catalog
                o_df = self.check_source(catalog_n, o_cat,
                                         i_alpha, i_delta)

                # If there is one saves data from input data
                if o_df.empty is not True and o_df['PM'].size == 1:
                    pm_mask = self.pm_filter(o_df, pm, prfs_d, confidence_)
                    if pm_mask:
                        if o_df['SOURCE_NUMBER'].size != 1:
                            boolean_l.append('False')
                        else:
                            boolean_l.append('True')
                        tmp_catalog.append(catalog_n)
                        tmp_source.append(o_df['SOURCE_NUMBER'].iloc[0])
                        tmp_pm.append(pm)
                        tmp_alpha.append(i_alpha)
                        tmp_delta.append(i_delta)
                else:
                    # Create a big region file with all faults
                    out_d['source'].append(source_)
                    out_d['CCD'].append(ccd_)
                    out_d['dither'].append(dither_)
                    out_d['catalog'].append(catalog_n)
                    out_d['PM'].append(pm)
                    out_d['alpha_j2000'].append(i_alpha)
                    out_d['delta_j2000'].append(i_delta)

                    boolean_l.append('False')

            # Total number
            idx = stats_d['PM'].index(pm)
            stats_d['total'][idx] += 1

            if len(tmp_source) is not 0:
                flag_detection, sources_number = all_same(tmp_source)
                if len(list(set(tmp_source))) > 1:
                    print tmp_source

#            if len(list(set(boolean_l))) == 1 and list(set(boolean_l))[0] == True:
            if flag_detection and sources_number >= 3:
                # print tmp_source
                idx = stats_d['PM'].index(pm)
                stats_d['right'][idx] += 1
                # print "True", boolean_l
                # print "catalog", tmp_catalog
                # print "source", tmp_source
                # print "pm", tmp_pm
                # print "tmp_alpha", tmp_alpha
                # print "tmp_delta", tmp_delta
            else:
                # self.create_regions(tmp_alpha, tmp_delta, source_)
                # self.show_regions()
                # print "False", boolean_l
                pass

        out_df = DataFrame(out_d)
        # Saves output to easily readable files
        out_df.to_csv('errors_{}_{}.csv'.format(idx_file, confidence_))
        out_df.to_csv('errors_{}_{}.reg'.format(idx_file, confidence_),
                      index=False, header=False, sep=" ")

        return stats_d

    def pm_filter(self, o_df, pm, prfs_d, confidence_):
        """

        :param o_df:
        :param pm:
        :param prfs_d:
        :param confidence_:
        :return:
        """
        pm_ranges = speeds_range(prfs_d, confidence_)
        pm_range = pm_ranges[pm]

        if pm_range[0] < float(o_df['PM']) < pm_range[1]:
            return True
        else:
            return False

    def check_source(self, catalog_n, o_cat, i_alpha, i_delta):
        """

        :param catalog_n:
        :param o_cat:
        :param i_alpha:
        :param i_delta:
        :return:
        """
        tolerance = 0.001

        o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
        o_df = o_df[o_df['ALPHA_J2000'] + tolerance > i_alpha]
        o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - tolerance]
        o_df = o_df[o_df['DELTA_J2000'] + tolerance > i_delta]
        o_df = o_df[i_delta > o_df['DELTA_J2000'] - tolerance]

        return o_df

    def create_regions(self, i_alpha_l, i_delta_l, source_):
        """

        :param i_alpha_l:
        :param i_delta_l:
        :param source_:
        :return:
        """
        i_dict = {'i_alpha': i_alpha_l, 'i_delta': i_delta_l}
        i_df = DataFrame(i_dict)
        i_df.to_csv('{}.csv'.format(source_), index=False,
                    header=False, sep=" ")

    def create_dict(self, scmp_cf, sex_cf, confidence_):
        """

        :param scmp_cf:
        :param sex_cf:
        :param confidence_:
        :return:
        """
        stats_keys = ['total', 'right', 'false',
                      'f_dr', 'f_pur', 'f_com']

        stats_d = {}
        stats_d['PM'] = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3,
                         1, 3, 10, 30, 100, 300]

        scamp_parameters = scmp_cf.split('_')
        sex_parameters = sex_cf.split('_')

        stats_d['crossid'] = []
        stats_d['pixscale'] = []
        stats_d['posangle'] = []
        stats_d['position'] = []
        stats_d['deblending'] = []
        stats_d['threshold'] = []
        stats_d['mincount'] = []
        stats_d['area'] = []
        stats_d['confidence'] = []

        for value_ in range(len(stats_d['PM'])):
            stats_d['crossid'].append(scamp_parameters[0])
            stats_d['pixscale'].append(scamp_parameters[1])
            stats_d['posangle'].append(scamp_parameters[2])
            stats_d['position'].append(scamp_parameters[3])
            stats_d['deblending'].append(sex_parameters[0])
            stats_d['threshold'].append(sex_parameters[1])
            stats_d['mincount'].append(sex_parameters[3])
            stats_d['area'].append(sex_parameters[4])
            # Confidence
            stats_d['confidence'].append(confidence_)

        for key_ in stats_keys:
            stats_d[key_] = []
            for value_ in range(len(stats_d['PM'])):
                stats_d[key_].append(0)

        # out dictionary
        out_keys = ['alpha_j2000', 'delta_j2000',
                    'catalog', 'PM', 'source', 'CCD', 'dither']
        out_d = {}

        for key_ in out_keys:
            out_d[key_] = []

        return (stats_d, out_d)
