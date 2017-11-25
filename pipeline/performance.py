#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.2 Now supports confidence intervals
- 0.3 Sextractor performance's added
- 0.4 Proper motion performance's added
- 0.5 Proper motion now creates an stats file

Todo:
    * Improve log messages

"""
from datetime import datetime, timedelta
from decimal import Decimal
from math import modf
from os import makedirs, path

from astropy import units as u
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from numpy import arange, array, median
from pandas import concat, read_csv, DataFrame
import statsmodels.api as sm

from misc import all_same, speeds_range
from regions import Create_regions


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0."
__maintainer__ = "Samuel Góngora García"
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
        errors_a_stars = []
        errors_b_stars = []
        errors_a_galaxies = []
        errors_b_galaxies = []
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
            i_cat_n = '/Cat_20-21_d'
            input_stars_d[d] = '{}{}{}.dat'.format(prfs_d['input_ref'],
                                                   i_cat_n, d)
        input_stars_d = Create_regions(input_stars_d,
                                       prfs_d).check_stars(True, True)

        # Creates a DataFrame from an input dictionary
        input_stars_l = []
        for key_ in input_stars_d.keys():
            input_stars_l.append(input_stars_d[key_])

        i_stars_df = concat(input_stars_l, axis=0)
        i_stars_df = i_stars_df.reset_index(drop=True)

        # Galaxies
        input_galaxies_d = {}
        for d in range(1, 5, 1):
            i_cat_n = '/Cat_20-21_d'
            input_galaxies_d[d] = '{}{}{}.dat'.format(prfs_d['input_ref'],
                                                      i_cat_n, d)
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

        mean_a_stars = array(errors_a_stars).mean()
        mean_b_stars = array(errors_b_stars).mean()

        mean_a_galaxies = array(errors_a_galaxies).mean()
        mean_b_galaxies = array(errors_b_galaxies).mean()

        mean_a_ssos = array(errors_a_ssos).mean()
        mean_b_ssos = array(errors_b_ssos).mean()

        std_a_stars = array(errors_a_stars).std()
        std_b_stars = array(errors_b_stars).std()

        std_a_galaxies = array(errors_a_galaxies).std()
        std_b_galaxies = array(errors_b_galaxies).std()

        std_a_ssos = array(errors_a_ssos).std()
        std_b_ssos = array(errors_b_ssos).std()

        stats_d = {'conf': sex_cf,
                   'm_a_stars': mean_a_stars, 'm_b_stars': mean_b_stars,
                   'm_a_gals': mean_a_galaxies, 'm_b_gals': mean_b_galaxies,
                   'm_a_ssos': mean_a_ssos, 'm_b_ssos': mean_b_ssos,
                   's_a_stars': std_a_stars, 's_b_stars': std_b_stars,
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
        self.bypassfilter = True
        self.plot = False

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

        i_df.to_csv('input.csv')

        # Open particular file!
        filt_n = 'filt_{}_{}_3.csv'.format(scmp_cf, mag)
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
                    pm_mask = self.pm_filter(o_df, i_pm, prfs_d,
                                             confidence_, self.bypassfilter)
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

        pm_d = {}
        for key_ in list(set(input_pm)):
            pm_d[key_] = []

        for idx, o_pm in enumerate(output_pm):
            i_pm = input_pm[idx]
            pm_d[i_pm].append(o_pm)

        stats_d = {'pm': [], 'mean': [], 'median': [], 'std': [],
                   'max': [], 'min': [], 'detected': []}

        # mean and std determination
        for key_ in pm_d.keys():
            stats_d['pm'].append(key_)
            stats_d['mean'].append(array(pm_d[key_]).mean())
            stats_d['median'].append(median(array(pm_d[key_])))
            stats_d['std'].append(array(pm_d[key_]).std())
            stats_d['max'].append(max(pm_d[key_]))
            stats_d['min'].append(min(pm_d[key_]))
            stats_d['detected'].append(len(pm_d[key_]))
            # print(key_, len(pm_d[key_]))
            # print(array(pm_d[key_]).mean())
            # print(median(array(pm_d[key_])))
            # print(array(pm_d[key_]).std())
            # print(max(pm_d[key_]))
            # print(min(pm_d[key_]))
            # print(" ")

        if self.plot:
            if not self.plot_comp(prfs_d, sex_cf, scmp_cf,
                                  input_pm, output_pm, confidence_):
                raise Exception

        return stats_d

    def pm_filter(self, o_df, pm, prfs_d, confidence_, bypassfilter):
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
        elif bypassfilter:
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
        tolerance = 0.0001389  # 0.5 arcseconds

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
        self.bypassfilter = True
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
        # For now any file is saved
        save = False
        self.prfs_d = prfs_d

        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {} and {}'.format(scmp_cf,
                                                                 sex_cf))
        input_d = {}
        for d in range(1, 5, 1):
            input_d[d] = '{}/Cat_20-21_d{}.dat'.format(prfs_d['input_ref'],
                                                       d)
        input_d = Create_regions(input_d, self.prfs_d).check_luca(True, True)

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
        filter_o_n = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
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
            tmp_epoch = []
            tmp_pm = []
            tmp_i_alpha = []
            tmp_i_delta = []
            tmp_o_alpha = []
            tmp_o_delta = []
            tmp_a_error = []
            tmp_b_error = []
            # Creates a flag for right detections
            # Initial value will set to False
            flag_detection = False

            # Iterate over each detection of each source
            for i, row in enumerate(cat.itertuples(), 1):
                # source_ = row.source
                # ccd_ = row.CCD
                # dither_ = row.dither_values
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
                    pm_mask = self.pm_filter(o_df, pm, confidence_)
                    if pm_mask:
                        if o_df['SOURCE_NUMBER'].size != 1:
                            boolean_l.append('False')
                        else:
                            boolean_l.append('True')
                        tmp_catalog.append(catalog_n)
                        source = int(o_df['SOURCE_NUMBER'].iloc[0])
                        tmp_source.append(source)
                        epoch = float(o_df['EPOCH'].iloc[0])
                        tmp_epoch.append(epoch)
                        tmp_pm.append(pm)
                        tmp_i_alpha.append(i_alpha)
                        tmp_i_delta.append(i_delta)
                        o_alpha = float(o_df['ALPHA_J2000'].iloc[0])
                        tmp_o_alpha.append(o_alpha)
                        o_delta = float(o_df['DELTA_J2000'].iloc[0])
                        tmp_o_delta.append(o_delta)
                        errora_world = float(o_df['ERRA_WORLD'].iloc[0])
                        tmp_a_error.append(errora_world)
                        errorb_world = float(o_df['ERRB_WORLD'].iloc[0])
                        tmp_b_error.append(errorb_world)
                else:
                    boolean_l.append('False')

            # Total number
            idx = stats_d['PM'].index(pm)
            stats_d['total'][idx] += 1

            order_mask = self.check_cat_order(tmp_catalog)

            if not order_mask:
                raise Exception

            if len(tmp_source) is not 0:
                flag_detection, sources_number = all_same(tmp_source)

            if flag_detection and sources_number >= 3:
                # print('tmp_source {}'.format(tmp_source))  # testing purposes
                idx = stats_d['PM'].index(pm)
                stats_d['right'][idx] += 1
                fitted_d = self.confidence(tmp_source[0], scmp_cf, sex_cf)

                pm = float(tmp_pm[0])
                if pm == 0.1:
                    self.plot(tmp_source[0], tmp_epoch,
                              tmp_i_alpha, tmp_i_delta,
                              tmp_o_alpha, tmp_o_delta,
                              tmp_a_error, tmp_b_error,
                              pm, fitted_d,
                              scmp_cf, sex_cf)
            else:
                pass

        return stats_d

    def pm_filter(self, o_df, pm, confidence_):
        """

        :param o_df:
        :param pm:
        :param confidence_:
        :return:
        """
        pm_ranges = speeds_range(self.prfs_d, confidence_)
        pm_range = pm_ranges[pm]

        if pm_range[0] < float(o_df['PM']) < pm_range[1]:
            return True
        elif self.bypassfilter:
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
        tolerance = 0.0001

        o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
        o_df = o_df[o_df['ALPHA_J2000'] + tolerance > i_alpha]
        o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - tolerance]
        o_df = o_df[o_df['DELTA_J2000'] + tolerance > i_delta]
        o_df = o_df[i_delta > o_df['DELTA_J2000'] - tolerance]

        return o_df

    def check_cat_order(self, cat_list):
        """

        :param cat_list:
        :return:
        """
        tmp_order = []

        for idx_cat in range(0, len(cat_list), 1):
            cat = cat_list[idx_cat]
            dither = cat - (cat / 4) * 4
            if dither == 0:
                dither = 4
            tmp_order.append(dither)

        if sorted(tmp_order) == tmp_order:
            return True
        else:
            return False

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
        stats_keys = ['total', 'right', 'false', 'f_dr', 'f_pur', 'f_com']

        stats_d = {'PM': [0.001, 0.003, 0.01, 0.03, 0.1, 0.3,
                          1, 3, 10, 30]}

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

        return stats_d, out_d

    def plot(self, source, tmp_epoch, tmp_i_alpha, tmp_i_delta,
             tmp_o_alpha, tmp_o_delta, tmp_a_error, tmp_b_error,
             pm, fitted_d, scmp_cf, sex_cf):
        """

        :param source:
        :param tmp_epoch:
        :param tmp_i_alpha:
        :param tmp_i_delta:
        :param tmp_o_alpha:
        :param tmp_o_delta:
        :param tmp_a_error:
        :param tmp_b_error:
        :param pm:
        :param fitted_d:
        :return:
        """
        # Set True to plot input and output data
        both = False

        output_path = '{}/mv/{}/{}/{}'.format(self.prfs_d['plots_dir'],
                                              scmp_cf, sex_cf, pm)
        if not path.exists(output_path):
            makedirs(output_path)

        if both:
            self.plot_both(tmp_epoch, tmp_i_alpha, tmp_i_delta,
                           tmp_o_alpha, tmp_o_delta, output_path,
                           source, tmp_a_error, tmp_b_error, pm, fitted_d)
        else:
            from plots import PlotConfidence
            mode = 'o'
            tmp_d = {'alpha': tmp_o_alpha, 'delta': tmp_o_delta,
                     'error_a': tmp_a_error, 'error_b': tmp_b_error,
                     'epoch': tmp_epoch}
            test = PlotConfidence(output_path, source, pm, mode,
                                  fitted_d, tmp_d)
            """
            self.plot_input(tmp_epoch, tmp_i_alpha, tmp_i_delta, output_path,
                            source, tmp_a_error, tmp_b_error, pm, fitted_d)
            self.plot_output(tmp_epoch, tmp_o_alpha, tmp_o_delta, output_path,
                             source, tmp_a_error, tmp_b_error, pm, fitted_d)
            """

    def plot_both(self, tmp_epoch, tmp_i_alpha, tmp_i_delta, tmp_o_alpha,
                  tmp_o_delta, output_path, source, tmp_a_error, tmp_b_error,
                  pm, fitted_d):
        """
        # TODO Change variables name
        1 - For y-axis, alpha values
        1.1 - Creates a timestamp object from alpha angle
        1.2 - Creates a datetime object fom timestamp
        1.3 - Meanwhile creates a list populate only with seconds information
        2 - For x-axis, epoch values

        :param tmp_epoch:
        :param tmp_i_alpha:
        :param tmp_i_delta:
        :param tmp_o_alpha:
        :param tmp_o_delta:
        :param output_path:
        :param source:
        :param tmp_a_error:
        :param tmp_b_error:
        :param pm:
        :param fitted_d:
        :return:
        """
        with PdfPages('{}/{}_{}_b.pdf'.format(output_path,
                                              source, pm)) as pdf:
            # ALPHA PARAMETERS
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)
            ax.set_title('catalog - chi_squared: {}'.format(fitted_d['ra']))

            i_alpha_tmpstmp = []
            i_alpha_seconds = []
            tmp_hour = []
            tmp_minute = []
            for alpha_ in tmp_i_alpha:
                a = Angle(alpha_, u.degree)
                hms = a.hms
                hour = int(hms[0])
                tmp_hour.append(hour)
                minute = int(hms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(hms[2]))
                i_alpha_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
                i_alpha_seconds.append('{}'.format(second))

            i_datetimes = [datetime.strptime(t, "%H:%M:%S.%f") for t in i_alpha_tmpstmp]

            o_alpha_tmpstmp = []
            o_alpha_seconds = []
            tmp_hour = []
            tmp_minute = []
            for alpha_ in tmp_o_alpha:
                a = Angle(alpha_, u.degree)
                hms = a.hms
                hour = int(hms[0])
                tmp_hour.append(hour)
                minute = int(hms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(hms[2]))
                o_alpha_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
                o_alpha_seconds.append('{}'.format(second))

            o_datetimes = [datetime.strptime(t, "%H:%M:%S.%f") for t in o_alpha_tmpstmp]

            # X-AXIS Shared between input and output catalogs
            epoch_seconds = []
            for epoch_ in tmp_epoch:
                frac, whole = modf(epoch_)
                seconds_ = frac * 365.25 * 24 * 60
                epoch_seconds.append(seconds_)

            x_difference = float(max(epoch_seconds)) - float(min(epoch_seconds))
            x_major_step = (x_difference / 3)
            x_minor_step = (x_difference / 6)

            # X-SCALE
            x_major_ticks = [epoch_seconds[0] - x_major_step,
                             epoch_seconds[0],
                             epoch_seconds[0] + x_major_step,
                             epoch_seconds[0] + x_major_step * 2,
                             epoch_seconds[0] + x_major_step * 3,
                             epoch_seconds[0] + x_major_step * 4]

            x_minor_ticks = [epoch_seconds[0] - x_minor_step * 2,
                             epoch_seconds[0] - x_minor_step,
                             epoch_seconds[0],
                             epoch_seconds[0] + x_minor_step,
                             epoch_seconds[0] + x_minor_step * 2,
                             epoch_seconds[0] + x_minor_step * 3,
                             epoch_seconds[0] + x_minor_step * 4,
                             epoch_seconds[0] + x_minor_step * 5,
                             epoch_seconds[0] + x_minor_step * 6,
                             epoch_seconds[0] + x_minor_step * 7,
                             epoch_seconds[0] + x_minor_step * 8]

            ax.set_xlabel('EPOCH')

            # Y-AXIS
            myformat = mdates.DateFormatter('%S.%f')
            ax.yaxis.set_major_formatter(myformat)

            # format elements in alpha_seconds list to floats
            i_alpha_seconds = [float(i) for i in i_alpha_seconds]
            i_alpha_seconds = [float("{0:.6f}".format(i)) for i in
                               i_alpha_seconds]

            # format elements in alpha_seconds list to floats
            o_alpha_seconds = [float(i) for i in o_alpha_seconds]
            o_alpha_seconds = [float("{0:.6f}".format(i)) for i in
                               o_alpha_seconds]

            # Gets maximum and minium value of alpha and delta (in seconds)
            max_alpha_i = max(i_alpha_seconds)
            max_alpha_o = max(o_alpha_seconds)
            max_alpha = max([max_alpha_i, max_alpha_o])
            min_alpha_i = min(i_alpha_seconds)
            min_alpha_o = min(o_alpha_seconds)
            min_alpha = min([min_alpha_i, min_alpha_o])

            # Gets the major step between ticks thought the difference
            # between the maximum and the minimum value of alpha
            difference = float(max_alpha) - float(min_alpha)
            major_step = (difference / 4)
            major_step_t = self.significant_l(major_step)

            # Rounds the float value of difference
            first_digit = '%.2E' % Decimal(major_step - major_step_t)

            if first_digit[0] == '-':
                first_digit = int(first_digit[1])
            else:
                first_digit = int(first_digit[0])

            # Redondea al alza el ultimo digio
            if first_digit > 5:
                last_digit = str(major_step_t)[-1]
                last_digit = int(last_digit)
                last_digit += 1
                major_step_t = list(str(major_step_t))
                major_step_t[-1] = last_digit
                major_step_t = [str(i) for i in major_step_t]
                major_step_t = ''.join(major_step_t)
                major_step_t = float(major_step_t)
            else:
                pass  # do nothing

            major_step = major_step_t
            minor_step = (major_step / 4)

            # Gets maximum decimal position of major step
            step_decimals = int(str(major_step)[::-1].find('.'))

            # Major step list starts two times before and end two times after
            # known values
            major_steps = arange(round(min_alpha,
                                       step_decimals) - major_step * 2,
                                 round(max_alpha,
                                       step_decimals) + major_step * 2,
                                 major_step)
            # Minor step list starts eight times before and end eight times
            # after know values
            minor_steps = arange(round(min_alpha,
                                       step_decimals) - minor_step * 8,
                                 round(max_alpha,
                                       step_decimals) + minor_step * 8,
                                 minor_step)

            # Formats major steps for be readable by matplotlib axis
            i_alpha_major_steps = []
            for idx_step, step_ in enumerate(major_steps):
                # Check if all hours/minutes are same
                if len(list(set(tmp_hour))) != 1:
                    raise Exception
                if len(list(set(tmp_minute))) != 1:
                    raise Exception

                # Gets hour and minute value
                hour = tmp_hour[0]
                minute = tmp_minute[0]
                i_alpha_major_steps.append('{}:{}:{}'.format(hour, minute, step_))

            # Creates a list of datatime objects
            # Sometimes, seconds are greater than 59, this values should be
            # filter in a better way
            try:
                i_alpha_steps_d = [datetime.strptime(t, "%H:%M:%S.%f") for t in i_alpha_major_steps]

                # Formats minor steps
                i_alpha_minor_steps = []
                for idx_step, step_ in enumerate(minor_steps):
                    # Check if all hours/minutes are same
                    if len(list(set(tmp_hour))) != 1:
                        raise Exception
                    if len(list(set(tmp_minute))) != 1:
                        raise Exception

                    # Gets hour and minute value
                    hour = tmp_hour[0]
                    minute = tmp_minute[0]
                    i_alpha_minor_steps.append('{}:{}:{}'.format(hour, minute, step_))

                # Creates a list of datetime objects
                i_alpha_minor_steps_d = [datetime.strptime(t, "%H:%M:%S.%f") for t in i_alpha_minor_steps]

                # x-ticks assignation
                ax.set_xticks(x_major_ticks, minor=False)
                ax.set_xticks(x_minor_ticks, minor=True)
                # x-ticks labels

                # y-ticks assignation
                ax.set_yticks(i_alpha_steps_d, minor=False)
                ax.set_yticks(i_alpha_minor_steps_d, minor=True)
                # y-ticks labels
                empty_string_labels = [''] * len(major_steps)
                ax.set_yticklabels(empty_string_labels, minor=False)
                ax.set_yticklabels(minor_steps, minor=True)

                # Formats grids
                ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    hour = datetime_.hour
                    minute = datetime_.minute
                    second = datetime_.second
                    msecond = datetime_.microsecond
                    alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
                                                              second, msecond)
                    error = float(
                        "{0:.6f}".format(tmp_a_error[idx_datetime_] * 3600))
                    error_str = '   error  {}"'.format(error)
                    # Annotate position and error associated
                    ax.annotate('{}\n{}'.format(alpha_str, error_str),
                                xy=(epoch_seconds[idx_datetime_], datetime_),
                                textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    ax.errorbar(epoch_seconds[idx_datetime_], datetime_,
                                yerr=timedelta(0, tmp_a_error[
                                    idx_datetime_] * 3600),
                                fmt='o', ecolor='g', capthick=2, elinewidth=4)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    hour = datetime_.hour
                    minute = datetime_.minute
                    second = datetime_.second
                    msecond = datetime_.microsecond
                    alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
                                                              second, msecond)
                    error = float(
                        "{0:.6f}".format(tmp_a_error[idx_datetime_] * 3600))
                    error_str = '   error  {}"'.format(error)
                    # Annotate position and error associated
                    ax.annotate('{}\n{}'.format(alpha_str, error_str),
                                xy=(epoch_seconds[idx_datetime_], datetime_),
                                textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    ax.errorbar(epoch_seconds[idx_datetime_], datetime_,
                                yerr=timedelta(0, tmp_a_error[
                                    idx_datetime_] * 3600),
                                fmt='o', ecolor='g', capthick=2, elinewidth=4)

                # Label creation
                y_label_ra = 'Right ascension (")\n'
                y_label_major_step = 'major step size {}"\n'.format(major_step)
                y_label_minor_step = 'minor step size {}"'.format(minor_step)
                ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                              y_label_minor_step))
                # Plots data
                ax.plot(epoch_seconds, o_datetimes, 'bs', markersize=6,
                        label='extracted position')
                ax.plot(epoch_seconds, i_datetimes, 'bs', markersize=6,
                        label='extracted position')
                # In order to avoid scientific notation plot should be redrawn
                ax = plt.gca()
                ax.get_xaxis().get_major_formatter().set_useOffset(False)

                # Plots a legend
                plt.legend(loc=0, ncol=2, borderaxespad=0.)
                plt.draw()

                # Saves the current figure in pdf file
                pdf.savefig()  # saves current figure
                plt.clf()  # clear current figure

                #
                #
                # DELTA PARAMETERS
                #
                #
                fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                ax = fig.add_subplot(1, 1, 1)
                ax.set_title('catalog - chi_squared: {}'.format(fitted_d['dec']))

                o_delta_tmpstmp = []
                o_delta_seconds = []
                for delta_ in tmp_i_delta:
                    d = Angle(delta_, u.degree)
                    dms = d.dms
                    o_delta_tmpstmp.append(
                        '{}:{}.{}'.format(int(dms[0]), int(dms[1]),
                                          float("{0:.6f}".format(dms[2]))))
                    o_delta_seconds.append(
                        '{}'.format(float("{0:.6f}".format(dms[2]))))

                ax.set_xlabel('EPOCH')

                # Y-SCALE
                # format elements in alpha_seconds list to floats
                o_delta_seconds = [float(i) for i in o_delta_seconds]
                o_delta_seconds = [float("{0:.6f}".format(i)) for i in
                                   o_delta_seconds]

                max_delta = max(o_delta_seconds)
                min_delta = min(o_delta_seconds)

                delta_difference = float(max_delta) - float(min_delta)
                delta_major_step = (delta_difference / 4)
                delta_major_step_t = self.significant_l(delta_major_step)

                # Rounds the float value of difference
                first_digit = '%.2E' % Decimal(delta_major_step - delta_major_step_t)

                if first_digit[0] == '-':
                    first_digit = int(first_digit[1])
                else:
                    first_digit = int(first_digit[0])

                # Redondea al alza el ultimo digio
                if first_digit > 5:
                    last_digit = str(delta_major_step_t)[-1]
                    last_digit = int(last_digit)
                    last_digit += 1
                    delta_major_step_t = list(str(delta_major_step_t))
                    delta_major_step_t[-1] = last_digit
                    delta_major_step_t = [str(i) for i in delta_major_step_t]
                    delta_major_step_t = ''.join(delta_major_step_t)
                    delta_major_step_t = float(delta_major_step_t)
                else:
                    pass  # do nothing

                delta_major_step = delta_major_step_t
                delta_minor_step = (delta_major_step / 4)

                # Gets maximum decimal position of major step
                step_decimals = int(str(delta_major_step)[::-1].find('.'))

                # Major step list starts two times before and end two times
                # after known values
                major_steps = arange(round(min_delta,
                                           step_decimals) - major_step * 2,
                                     round(max_delta,
                                           step_decimals) + major_step * 2,
                                     delta_major_step)
                # Minor step list starts eight times before and end eight times
                # after know values
                minor_steps = arange(round(min_delta,
                                           step_decimals) - minor_step * 8,
                                     round(max_delta,
                                           step_decimals) + minor_step * 8,
                                     delta_minor_step)

                # x-ticks assignation
                ax.set_xticks(x_major_ticks, minor=False)
                ax.set_xticks(x_minor_ticks, minor=True)
                # x-ticks labels
                # TODO
                # y-ticks assignation
                # rounds float to six decimal position
                # y_minor_ticks = [float("{0:.6f}".format(i)) for i in
                #                  y_minor_ticks]
                ax.set_yticks(major_steps, minor=False)
                ax.set_yticks(minor_steps, minor=True)
                # y-ticks labels
                empty_string_labels = [''] * len(major_steps)
                ax.set_yticklabels(empty_string_labels, minor=False)
                # Converts floats into strings
                # y_minor_ticks = [str(i) for i in y_minor_ticks]
                ax.set_yticklabels(minor_steps, minor=True)

                # Formats grids
                ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                # Input annotations
                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    # variables for position and error representation
                    delta_position = o_delta_tmpstmp[idx_datetime_]
                    delta_str = '   delta {}'.format(delta_position)
                    error = tmp_b_error[idx_datetime_] * 3600  # error on y
                    error_fmt = float("{0:.6f}".format(error))  # error rounded
                    error_str = '   error {}'.format(error_fmt)
                    # annotation
                    ax.annotate('{}\n{}'.format(delta_str, error_str),
                                xy=(x, y),
                                textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    error = float(tmp_b_error[idx_datetime_] * 3600)  # error y
                    ax.errorbar(x, y, yerr=error, fmt='o',
                                ecolor='g', capthick=2, elinewidth=4)

                # Output annotations
                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    # variables for position and error representation
                    delta_position = o_delta_tmpstmp[idx_datetime_]
                    delta_str = '   delta {}'.format(delta_position)
                    error = tmp_b_error[idx_datetime_] * 3600  # error on y
                    error_fmt = float("{0:.6f}".format(error))  # error rounded
                    error_str = '   error {}'.format(error_fmt)
                    # annotation
                    ax.annotate('{}\n{}'.format(delta_str, error_str),
                                xy=(x, y),
                                textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    error = float(tmp_b_error[idx_datetime_] * 3600)  # error y
                    ax.errorbar(x, y, yerr=error, fmt='o',
                                ecolor='g', capthick=2, elinewidth=4)
                # Label creation
                y_label_ra = 'Declination (")\n'
                y_label_major_step = 'major step size {}"\n'.format(delta_major_step)
                y_label_minor_step = 'minor step size {}"'.format(delta_minor_step)
                ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                              y_label_minor_step))

                ax.plot(epoch_seconds, o_delta_seconds, 'bs', markersize=6,
                        label='extracted position')
                ax.plot(epoch_seconds, o_datetimes, 'bs', markersize=6,
                        label='extracted position')

                ax = plt.gca()
                ax.get_xaxis().get_major_formatter().set_useOffset(False)
                plt.legend(loc=0, ncol=2, borderaxespad=0.)
                plt.draw()

                pdf.savefig()  # saves current figure
                plt.clf()  # clear current figure
            except ValueError:
                print('wrong value in source {}'.format(source))

    def plot_input(self, tmp_epoch, tmp_i_alpha, tmp_i_delta, output_path,
                   source, tmp_a_error, tmp_b_error, pm, fitted_d):
        """
        # TODO Change variables name
        1 - For y-axis, alpha values
        1.1 - Creates a timestamp object from alpha angle
        1.2 - Creates a datetime object fom timestamp
        1.3 - Meanwhile creates a list populate only with seconds information
        2 - For x-axis, epoch values

        :param tmp_epoch:
        :param tmp_i_alpha: list of alpha positions of sources (degrees)
        :param tmp_i_delta: list of delta positions of sources (degrees)
        :param output_path: output path for desired pm and configuration
        :param source: source number
        :param tmp_a_error:
        :param tmp_b_error:
        :param pm:
        :param fitted_d:
        :return:
        """
        with PdfPages('{}/{}_{}_i.pdf'.format(output_path,
                                              source, pm)) as pdf:
            # ALPHA PARAMETERS
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)
            ax.set_title('catalog - chi_squared: {}'.format(fitted_d['ra']))

            i_alpha_tmpstmp = []
            i_alpha_seconds = []
            tmp_hour = []
            tmp_minute = []
            for alpha_ in tmp_i_alpha:
                a = Angle(alpha_, u.degree)
                hms = a.hms
                hour = int(hms[0])
                tmp_hour.append(hour)
                minute = int(hms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(hms[2]))
                i_alpha_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
                i_alpha_seconds.append('{}'.format(second))

            i_datetimes = [datetime.strptime(t, "%H:%M:%S.%f") for t in i_alpha_tmpstmp]

            # X-AXIS Shared between input and output catalogs
            epoch_seconds = []
            for epoch_ in tmp_epoch:
                frac, whole = modf(epoch_)
                seconds_ = frac * 365.25 * 24 * 60
                epoch_seconds.append(seconds_)

            x_difference = float(max(epoch_seconds)) - float(min(epoch_seconds))
            x_major_step = (x_difference / 3)
            x_minor_step = (x_difference / 6)

            # X-SCALE
            x_major_ticks = [epoch_seconds[0] - x_major_step,
                             epoch_seconds[0],
                             epoch_seconds[0] + x_major_step,
                             epoch_seconds[0] + x_major_step * 2,
                             epoch_seconds[0] + x_major_step * 3,
                             epoch_seconds[0] + x_major_step * 4]

            x_minor_ticks = [epoch_seconds[0] - x_minor_step * 2,
                             epoch_seconds[0] - x_minor_step,
                             epoch_seconds[0],
                             epoch_seconds[0] + x_minor_step,
                             epoch_seconds[0] + x_minor_step * 2,
                             epoch_seconds[0] + x_minor_step * 3,
                             epoch_seconds[0] + x_minor_step * 4,
                             epoch_seconds[0] + x_minor_step * 5,
                             epoch_seconds[0] + x_minor_step * 6,
                             epoch_seconds[0] + x_minor_step * 7,
                             epoch_seconds[0] + x_minor_step * 8]

            ax.set_xlabel('EPOCH')

            # Y-AXIS
            myformat = mdates.DateFormatter('%S.%f')
            ax.yaxis.set_major_formatter(myformat)

            # format elements in alpha_seconds list to floats
            i_alpha_seconds = [float(i) for i in i_alpha_seconds]
            i_alpha_seconds = [float("{0:.6f}".format(i)) for i in
                               i_alpha_seconds]

            # Gets maximum and minium value of alpha (in seconds)
            max_alpha = max(i_alpha_seconds)
            min_alpha = min(i_alpha_seconds)

            # Gets the major step between ticks thought the difference
            # between the maximum and the minimum value of alpha
            difference = float(max_alpha) - float(min_alpha)
            major_step = (difference / 4)
            major_step_t = self.significant_l(major_step)

            # Rounds the float value of difference
            first_digit = '%.2E' % Decimal(major_step - major_step_t)

            if first_digit[0] == '-':
                first_digit = int(first_digit[1])
            else:
                first_digit = int(first_digit[0])

            # Redondea al alza el ultimo digio
            if first_digit > 5:
                last_digit = str(major_step_t)[-1]
                last_digit = int(last_digit)
                last_digit += 1
                major_step_t = list(str(major_step_t))
                major_step_t[-1] = last_digit
                major_step_t = [str(i) for i in major_step_t]
                major_step_t = ''.join(major_step_t)
                major_step_t = float(major_step_t)
            else:
                pass  # do nothing

            major_step = major_step_t
            minor_step = (major_step / 4)

            # Gets maximum decimal position of major step
            step_decimals = int(str(major_step)[::-1].find('.'))

            # Major step list starts two times before and end two times after
            # known values
            major_steps = arange(round(min_alpha,
                                       step_decimals) - major_step * 2,
                                 round(max_alpha,
                                       step_decimals) + major_step * 2,
                                 major_step)
            # Minor step list starts eight times before and end eight times
            # after know values
            minor_steps = arange(round(min_alpha,
                                       step_decimals) - minor_step * 8,
                                 round(max_alpha,
                                       step_decimals) + minor_step * 8,
                                 minor_step)

            # Formats major steps for be readable by matplotlib axis
            i_alpha_major_steps = []
            for idx_step, step_ in enumerate(major_steps):
                # Check if all hours/minutes are same
                if len(list(set(tmp_hour))) != 1:
                    raise Exception
                if len(list(set(tmp_minute))) != 1:
                    raise Exception

                # Gets hour and minute value
                hour = tmp_hour[0]
                minute = tmp_minute[0]
                i_alpha_major_steps.append('{}:{}:{}'.format(hour, minute, step_))

            # Creates a list of datatime objects
            # Sometimes, seconds are greater than 59, this values should be
            # filter in a better way
            try:
                i_alpha_steps_d = [datetime.strptime(t, "%H:%M:%S.%f") for t in i_alpha_major_steps]

                # Formats minor steps
                i_alpha_minor_steps = []
                for idx_step, step_ in enumerate(minor_steps):
                    # Check if all hours/minutes are same
                    if len(list(set(tmp_hour))) != 1:
                        raise Exception
                    if len(list(set(tmp_minute))) != 1:
                        raise Exception

                    # Gets hour and minute value
                    hour = tmp_hour[0]
                    minute = tmp_minute[0]
                    i_alpha_minor_steps.append('{}:{}:{}'.format(hour, minute, step_))

                # Creates a list of datetime objects
                i_alpha_minor_steps_d = [datetime.strptime(t, "%H:%M:%S.%f") for t in i_alpha_minor_steps]

                # x-ticks assignation
                ax.set_xticks(x_major_ticks, minor=False)
                ax.set_xticks(x_minor_ticks, minor=True)
                # x-ticks labels

                # y-ticks assignation
                ax.set_yticks(i_alpha_steps_d, minor=False)
                ax.set_yticks(i_alpha_minor_steps_d, minor=True)
                # y-ticks labels
                empty_string_labels = [''] * len(major_steps)
                ax.set_yticklabels(empty_string_labels, minor=False)
                ax.set_yticklabels(minor_steps, minor=True)

                # Formats grids
                ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    hour = datetime_.hour
                    minute = datetime_.minute
                    second = datetime_.second
                    msecond = datetime_.microsecond
                    alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
                                                              second, msecond)
                    error = float(
                        "{0:.6f}".format(tmp_a_error[idx_datetime_] * 3600))
                    error_str = '   error  {}"'.format(error)
                    # Annotate position and error associated
                    ax.annotate('{}\n{}'.format(alpha_str, error_str),
                                xy=(epoch_seconds[idx_datetime_], datetime_),
                                textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    ax.errorbar(epoch_seconds[idx_datetime_], datetime_,
                                yerr=timedelta(0, tmp_a_error[
                                    idx_datetime_] * 3600),
                                fmt='o', ecolor='g', capthick=2, elinewidth=4)

                # Label creation
                y_label_ra = 'Right ascension (")\n'
                y_label_major_step = 'major step size {}"\n'.format(major_step)
                y_label_minor_step = 'minor step size {}"'.format(minor_step)
                ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                              y_label_minor_step))
                # Plots data
                ax.plot(epoch_seconds, i_datetimes, 'bs', markersize=6,
                        label='extracted position')
                # In order to avoid scientific notation plot should be redrawn
                ax = plt.gca()
                ax.get_xaxis().get_major_formatter().set_useOffset(False)

                # Plots a legend
                plt.legend(loc=0, ncol=2, borderaxespad=0.)
                plt.draw()

                # Saves the current figure in pdf file
                pdf.savefig()  # saves current figure
                plt.clf()  # clear current figure

                #
                #
                # DELTA PARAMETERS
                #
                #
                fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                ax = fig.add_subplot(1, 1, 1)
                ax.set_title('catalog - chi_squared: {}'.format(fitted_d['dec']))

                o_delta_tmpstmp = []
                o_delta_seconds = []
                for delta_ in tmp_i_delta:
                    d = Angle(delta_, u.degree)
                    dms = d.dms
                    o_delta_tmpstmp.append(
                        '{}:{}.{}'.format(int(dms[0]), int(dms[1]),
                                          float("{0:.6f}".format(dms[2]))))
                    o_delta_seconds.append(
                        '{}'.format(float("{0:.6f}".format(dms[2]))))

                ax.set_xlabel('EPOCH')

                # Y-SCALE
                # format elements in alpha_seconds list to floats
                o_delta_seconds = [float(i) for i in o_delta_seconds]
                o_delta_seconds = [float("{0:.6f}".format(i)) for i in
                                   o_delta_seconds]

                max_delta = max(o_delta_seconds)
                min_delta = min(o_delta_seconds)

                delta_difference = float(max_delta) - float(min_delta)
                delta_major_step = (delta_difference / 4)
                delta_major_step_t = self.significant_l(delta_major_step)

                # Rounds the float value of difference
                first_digit = '%.2E' % Decimal(delta_major_step - delta_major_step_t)

                if first_digit[0] == '-':
                    first_digit = int(first_digit[1])
                else:
                    first_digit = int(first_digit[0])

                # Redondea al alza el ultimo digio
                if first_digit > 5:
                    last_digit = str(delta_major_step_t)[-1]
                    last_digit = int(last_digit)
                    last_digit += 1
                    delta_major_step_t = list(str(delta_major_step_t))
                    delta_major_step_t[-1] = last_digit
                    delta_major_step_t = [str(i) for i in delta_major_step_t]
                    delta_major_step_t = ''.join(delta_major_step_t)
                    delta_major_step_t = float(delta_major_step_t)
                else:
                    pass  # do nothing

                delta_major_step = delta_major_step_t
                delta_minor_step = (delta_major_step / 4)

                # Gets maximum decimal position of major step
                step_decimals = int(str(delta_major_step)[::-1].find('.'))

                # Major step list starts two times before and end two times
                # after known values
                major_steps = arange(round(min_delta,
                                           step_decimals) - major_step * 2,
                                     round(max_delta,
                                           step_decimals) + major_step * 2,
                                     delta_major_step)
                # Minor step list starts eight times before and end eight times
                # after know values
                minor_steps = arange(round(min_delta,
                                           step_decimals) - minor_step * 8,
                                     round(max_delta,
                                           step_decimals) + minor_step * 8,
                                     delta_minor_step)

                # x-ticks assignation
                ax.set_xticks(x_major_ticks, minor=False)
                ax.set_xticks(x_minor_ticks, minor=True)
                # x-ticks labels
                # TODO
                # y-ticks assignation
                # rounds float to six decimal position
                # y_minor_ticks = [float("{0:.6f}".format(i)) for i in
                #                  y_minor_ticks]
                ax.set_yticks(major_steps, minor=False)
                ax.set_yticks(minor_steps, minor=True)
                # y-ticks labels
                empty_string_labels = [''] * len(major_steps)
                ax.set_yticklabels(empty_string_labels, minor=False)
                # Converts floats into strings
                # y_minor_ticks = [str(i) for i in y_minor_ticks]
                ax.set_yticklabels(minor_steps, minor=True)

                # Formats grids
                ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    # variables for position and error representation
                    delta_position = o_delta_tmpstmp[idx_datetime_]
                    delta_str = '   delta {}'.format(delta_position)
                    error = tmp_b_error[idx_datetime_] * 3600  # error on y
                    error_fmt = float("{0:.6f}".format(error))  # error rounded
                    error_str = '   error {}'.format(error_fmt)
                    # annotation
                    ax.annotate('{}\n{}'.format(delta_str, error_str),
                                xy=(x, y),
                                textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    error = float(tmp_b_error[idx_datetime_] * 3600)  # error y
                    ax.errorbar(x, y, yerr=error, fmt='o',
                                ecolor='g', capthick=2, elinewidth=4)

                # Label creation
                y_label_ra = 'Declination (")\n'
                y_label_major_step = 'major step size {}"\n'.format(delta_major_step)
                y_label_minor_step = 'minor step size {}"'.format(delta_minor_step)
                ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                              y_label_minor_step))

                ax.plot(epoch_seconds, o_delta_seconds, 'bs', markersize=6,
                        label='extracted position')

                ax = plt.gca()
                ax.get_xaxis().get_major_formatter().set_useOffset(False)
                plt.legend(loc=0, ncol=2, borderaxespad=0.)
                plt.draw()

                pdf.savefig()  # saves current figure
                plt.clf()  # clear current figure
            except ValueError:
                print('wrong value in source {}'.format(source))

    def plot_output(self, tmp_epoch, tmp_o_alpha, tmp_o_delta, output_path,
                    source, tmp_a_error, tmp_b_error, pm, fitted_d):
        """

        1 - For y-axis, alpha values
        1.1 - Creates a timestamp object from alpha angle
        1.2 - Creates a datetime object fom timestamp
        1.3 - Meanwhile creates a list populate only with seconds information
        2 - For x-axis, epoch values

        :param tmp_epoch:
        :param tmp_o_alpha: list of alpha positions of sources (degrees)
        :param tmp_o_delta: list of delta positions of sources (degrees)
        :param output_path: output path for desired pm and configuration
        :param source: source number
        :param tmp_a_error:
        :param tmp_b_error:
        :param pm:
        :param fitted_d:
        :return:
        """
        with PdfPages('{}/{}_{}_o.pdf'.format(output_path,
                                              source, pm)) as pdf:
            # ALPHA PARAMETERS
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)
            ax.set_title('extracted - chi_squared: {}'.format(fitted_d['ra']))

            o_alpha_tmpstmp = []
            o_alpha_seconds = []
            tmp_hour = []
            tmp_minute = []
            for alpha_ in tmp_o_alpha:
                a = Angle(alpha_, u.degree)
                hms = a.hms
                hour = int(hms[0])
                tmp_hour.append(hour)
                minute = int(hms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(hms[2]))
                o_alpha_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
                o_alpha_seconds.append('{}'.format(second))

            o_datetimes = [datetime.strptime(t, "%H:%M:%S.%f") for t in o_alpha_tmpstmp]

            # X-AXIS Shared between input and output catalogs
            epoch_seconds = []
            for epoch_ in tmp_epoch:
                frac, whole = modf(epoch_)
                seconds_ = frac * 365.25 * 24 * 60
                epoch_seconds.append(seconds_)

            x_difference = float(max(epoch_seconds)) - float(min(epoch_seconds))
            x_major_step = (x_difference / 3)
            x_minor_step = (x_difference / 6)

            # X-SCALE
            x_major_ticks = [epoch_seconds[0] - x_major_step,
                             epoch_seconds[0],
                             epoch_seconds[0] + x_major_step,
                             epoch_seconds[0] + x_major_step * 2,
                             epoch_seconds[0] + x_major_step * 3,
                             epoch_seconds[0] + x_major_step * 4]

            x_minor_ticks = [epoch_seconds[0] - x_minor_step * 2,
                             epoch_seconds[0] - x_minor_step,
                             epoch_seconds[0],
                             epoch_seconds[0] + x_minor_step,
                             epoch_seconds[0] + x_minor_step * 2,
                             epoch_seconds[0] + x_minor_step * 3,
                             epoch_seconds[0] + x_minor_step * 4,
                             epoch_seconds[0] + x_minor_step * 5,
                             epoch_seconds[0] + x_minor_step * 6,
                             epoch_seconds[0] + x_minor_step * 7,
                             epoch_seconds[0] + x_minor_step * 8]

            ax.set_xlabel('EPOCH')

            # Y-AXIS
            myformat = mdates.DateFormatter('%S.%f')
            ax.yaxis.set_major_formatter(myformat)

            # format elements in alpha_seconds list to floats
            o_alpha_seconds = [float(i) for i in o_alpha_seconds]
            o_alpha_seconds = [float("{0:.6f}".format(i)) for i in
                               o_alpha_seconds]

            # Gets maximum and minium value of alpha (in seconds)
            max_alpha = max(o_alpha_seconds)
            min_alpha = min(o_alpha_seconds)

            # Gets the major step between ticks thought the difference
            # between the maximum and the minimum value of alpha
            difference = float(max_alpha) - float(min_alpha)
            major_step = (difference / 4)
            major_step_t = self.significant_l(major_step)

            # Rounds the float value of difference
            first_digit = '%.2E' % Decimal(major_step - major_step_t)

            if first_digit[0] == '-':
                first_digit = int(first_digit[1])
            else:
                first_digit = int(first_digit[0])

            # Redondea al alza el ultimo digio
            if first_digit > 5:
                last_digit = str(major_step_t)[-1]
                last_digit = int(last_digit)
                last_digit += 1
                major_step_t = list(str(major_step_t))
                major_step_t[-1] = last_digit
                major_step_t = [str(i) for i in major_step_t]
                major_step_t = ''.join(major_step_t)
                major_step_t = float(major_step_t)
            else:
                pass  # do nothing

            major_step = major_step_t
            minor_step = (major_step / 4)

            # Gets maximum decimal position of major step
            step_decimals = int(str(major_step)[::-1].find('.'))

            # Major step list starts two times before and end two times after
            # known values
            major_steps = arange(round(min_alpha,
                                       step_decimals) - major_step * 2,
                                 round(max_alpha,
                                       step_decimals) + major_step * 2,
                                 major_step)
            # Minor step list starts eight times before and end eight times
            # after know values
            minor_steps = arange(round(min_alpha,
                                       step_decimals) - minor_step * 8,
                                 round(max_alpha,
                                       step_decimals) + minor_step * 8,
                                 minor_step)

            # Formats major steps for be readable by matplotlib axis
            o_alpha_major_steps = []
            for idx_step, step_ in enumerate(major_steps):
                # Check if all hours/minutes are same
                if len(list(set(tmp_hour))) != 1:
                    raise Exception
                if len(list(set(tmp_minute))) != 1:
                    raise Exception

                # Gets hour and minute value
                hour = tmp_hour[0]
                minute = tmp_minute[0]
                o_alpha_major_steps.append('{}:{}:{}'.format(hour, minute, step_))

            # Creates a list of datatime objects
            # Sometimes, seconds are greater than 59, this values should be
            # filter in a better way
            try:
                o_alpha_steps_d = [datetime.strptime(t, "%H:%M:%S.%f") for t in o_alpha_major_steps]

                # Formats minor steps
                o_alpha_minor_steps = []
                for idx_step, step_ in enumerate(minor_steps):
                    # Check if all hours/minutes are same
                    if len(list(set(tmp_hour))) != 1:
                        raise Exception
                    if len(list(set(tmp_minute))) != 1:
                        raise Exception

                    # Gets hour and minute value
                    hour = tmp_hour[0]
                    minute = tmp_minute[0]
                    o_alpha_minor_steps.append('{}:{}:{}'.format(hour, minute, step_))

                # Creates a list of datetime objects
                o_alpha_minor_steps_d = [datetime.strptime(t, "%H:%M:%S.%f") for t in o_alpha_minor_steps]

                # x-ticks assignation
                ax.set_xticks(x_major_ticks, minor=False)
                ax.set_xticks(x_minor_ticks, minor=True)
                # x-ticks labels

                # y-ticks assignation
                ax.set_yticks(o_alpha_steps_d, minor=False)
                ax.set_yticks(o_alpha_minor_steps_d, minor=True)
                # y-ticks labels
                empty_string_labels = [''] * len(major_steps)
                ax.set_yticklabels(empty_string_labels, minor=False)
                ax.set_yticklabels(minor_steps, minor=True)

                # Formats grids
                ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    hour = datetime_.hour
                    minute = datetime_.minute
                    second = datetime_.second
                    msecond = datetime_.microsecond
                    alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
                                                              second, msecond)
                    error = float(
                        "{0:.6f}".format(tmp_a_error[idx_datetime_] * 3600))
                    error_str = '   error  {}"'.format(error)
                    # Annotate position and error associated
                    ax.annotate('{}\n{}'.format(alpha_str, error_str),
                                xy=(epoch_seconds[idx_datetime_], datetime_),
                                textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    ax.errorbar(epoch_seconds[idx_datetime_], datetime_,
                                yerr=timedelta(0, tmp_a_error[
                                    idx_datetime_] * 3600),
                                fmt='o', ecolor='g', capthick=2, elinewidth=4)

                # Label creation
                y_label_ra = 'Right ascension (")\n'
                y_label_major_step = 'major step size {}"\n'.format(major_step)
                y_label_minor_step = 'minor step size {}"'.format(minor_step)
                ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                              y_label_minor_step))
                # Plots data
                ax.plot(epoch_seconds, o_datetimes, 'bs', markersize=6,
                        label='extracted position')
                # In order to avoid scientific notation plot should be redrawn
                ax = plt.gca()
                ax.get_xaxis().get_major_formatter().set_useOffset(False)

                # Plots a legend
                plt.legend(loc=0, ncol=2, borderaxespad=0.)
                plt.draw()

                # Saves the current figure in pdf file
                pdf.savefig()  # saves current figure
                plt.clf()  # clear current figure

                #
                #
                # DELTA PARAMETERS
                #
                #
                fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
                ax = fig.add_subplot(1, 1, 1)
                ax.set_title('extracted - chi_squared: {}'.format(fitted_d['dec']))

                o_delta_tmpstmp = []
                o_delta_seconds = []
                for delta_ in tmp_o_delta:
                    d = Angle(delta_, u.degree)
                    dms = d.dms
                    o_delta_tmpstmp.append(
                        '{}:{}.{}'.format(int(dms[0]), int(dms[1]),
                                          float("{0:.6f}".format(dms[2]))))
                    o_delta_seconds.append(
                        '{}'.format(float("{0:.6f}".format(dms[2]))))

                ax.set_xlabel('EPOCH')

                # Y-SCALE
                # format elements in alpha_seconds list to floats
                o_delta_seconds = [float(i) for i in o_delta_seconds]
                o_delta_seconds = [float("{0:.6f}".format(i)) for i in
                                   o_delta_seconds]

                max_delta = max(o_delta_seconds)
                min_delta = min(o_delta_seconds)

                delta_difference = float(max_delta) - float(min_delta)
                delta_major_step = (delta_difference / 4)
                delta_major_step_t = self.significant_l(delta_major_step)

                # Rounds the float value of difference
                first_digit = '%.2E' % Decimal(delta_major_step - delta_major_step_t)

                if first_digit[0] == '-':
                    first_digit = int(first_digit[1])
                else:
                    first_digit = int(first_digit[0])

                # Redondea al alza el ultimo digio
                if first_digit > 5:
                    last_digit = str(delta_major_step_t)[-1]
                    last_digit = int(last_digit)
                    last_digit += 1
                    delta_major_step_t = list(str(delta_major_step_t))
                    delta_major_step_t[-1] = last_digit
                    delta_major_step_t = [str(i) for i in delta_major_step_t]
                    delta_major_step_t = ''.join(delta_major_step_t)
                    delta_major_step_t = float(delta_major_step_t)
                else:
                    pass  # do nothing

                delta_major_step = delta_major_step_t
                delta_minor_step = (delta_major_step / 4)

                # Gets maximum decimal position of major step
                step_decimals = int(str(delta_major_step)[::-1].find('.'))

                # Major step list starts two times before and end two times
                # after known values
                major_steps = arange(round(min_delta,
                                           step_decimals) - major_step * 2,
                                     round(max_delta,
                                           step_decimals) + major_step * 2,
                                     delta_major_step)
                # Minor step list starts eight times before and end eight times
                # after know values
                minor_steps = arange(round(min_delta,
                                           step_decimals) - minor_step * 8,
                                     round(max_delta,
                                           step_decimals) + minor_step * 8,
                                     delta_minor_step)

                # x-ticks assignation
                ax.set_xticks(x_major_ticks, minor=False)
                ax.set_xticks(x_minor_ticks, minor=True)
                # x-ticks labels
                # TODO
                # y-ticks assignation
                # rounds float to six decimal position
                # y_minor_ticks = [float("{0:.6f}".format(i)) for i in
                #                  y_minor_ticks]
                ax.set_yticks(major_steps, minor=False)
                ax.set_yticks(minor_steps, minor=True)
                # y-ticks labels
                empty_string_labels = [''] * len(major_steps)
                ax.set_yticklabels(empty_string_labels, minor=False)
                # Converts floats into strings
                # y_minor_ticks = [str(i) for i in y_minor_ticks]
                ax.set_yticklabels(minor_steps, minor=True)

                # Formats grids
                ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    # variables for position and error representation
                    delta_position = o_delta_tmpstmp[idx_datetime_]
                    delta_str = '   delta {}'.format(delta_position)
                    error = tmp_b_error[idx_datetime_] * 3600  # error on y
                    error_fmt = float("{0:.6f}".format(error))  # error rounded
                    error_str = '   error {}'.format(error_fmt)
                    # annotation
                    ax.annotate('{}\n{}'.format(delta_str, error_str),
                                xy=(x, y),
                                textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    error = float(tmp_b_error[idx_datetime_] * 3600)  # error y
                    ax.errorbar(x, y, yerr=error, fmt='o',
                                ecolor='g', capthick=2, elinewidth=4)

                # Label creation
                y_label_ra = 'Declination (")\n'
                y_label_major_step = 'major step size {}"\n'.format(delta_major_step)
                y_label_minor_step = 'minor step size {}"'.format(delta_minor_step)
                ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                              y_label_minor_step))

                ax.plot(epoch_seconds, o_delta_seconds, 'bs', markersize=6,
                        label='extracted position')

                ax = plt.gca()
                ax.get_xaxis().get_major_formatter().set_useOffset(False)
                plt.legend(loc=0, ncol=2, borderaxespad=0.)
                plt.draw()

                pdf.savefig()  # saves current figure
                plt.clf()  # clear current figure
            except ValueError:
                print('wrong value in source {}'.format(source))

    def confidence(self, source, scmp_cf, sex_cf):
        """

        :param source:
        :return:
        """
        catalogs_dir = self.prfs_d['catalogs_dir']
        configuration_dir = '/{}/{}'.format(sex_cf, scmp_cf)
        cat_name = '/full_{}_20-21_1.cat'.format(scmp_cf)
        hdu_list = fits.open('{}{}{}'.format(catalogs_dir,
                                             configuration_dir, cat_name))
        db = Table(hdu_list[2].data).to_pandas()

        ra = db.loc[db['SOURCE_NUMBER'] == source, 'ALPHA_J2000'].tolist()
        dec = db.loc[db['SOURCE_NUMBER'] == source, 'DELTA_J2000'].tolist()
        epoch = db.loc[db['SOURCE_NUMBER'] == source, 'EPOCH'].tolist()

        edims = []
        epochs = []
        dimensions = []
        fitted_d = {}

        for dimension in [ra, dec]:
            x = array(epoch)
            epochs.append(x)  # epochs list
            y = array(dimension)
            dimensions.append(y)  # dimensions list
            if dimension == ra:
                sigma = db.loc[db['SOURCE_NUMBER'] == source,
                               'ERRA_WORLD'].tolist()
                dimension_ = 'ra'
            if dimension == dec:
                sigma = db.loc[db['SOURCE_NUMBER'] == source,
                               'ERRB_WORLD'].tolist()
                dimension_ = 'dec'
            edim = array([1 / var for var in sigma])
            edims.append(edim)
            x = sm.add_constant(x)
            # model = sm.WLS(y, x, weigths=edim)
            model = sm.WLS(y, x)
            fitted = model.fit()
            fitted = fitted.rsquared
            fitted_d[dimension_] = fitted

        return fitted_d

    def significant_l(self, number):
        """

        :param number:
        :return:
        """
        l = len(str(number))   ####make sure there is enough precision
        a = ('%.' + str(l) + 'E') % Decimal(number)
        significate_d = a.split(".")[0]
        times = a.split("E")[1]
        result = int(significate_d) * (10 ** int(times))

        return result


"""
    def plot_input(self, tmp_epoch, tmp_o_alpha, tmp_o_delta, output_path,
                   source, tmp_a_error, tmp_b_error, pm, fitted_d):

        with PdfPages('{}/{}_{}_i.pdf'.format(output_path,
                                              source, pm)) as pdf:
            # ALPHA PARAMETERS
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)
            ax.set_title('chi_squared: {}'.format(fitted_d['ra']))

            o_alpha_tmpstmp = []
            o_alpha_seconds = []
            for alpha_ in tmp_o_alpha:
                a = Angle(alpha_, u.degree)
                hms = a.hms
                hour = int(hms[0])
                minute = int(hms[1])
                second = float("{0:.6f}".format(hms[2]))
                o_alpha_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
                o_alpha_seconds.append('{}'.format(second))

            o_datetimes = [datetime.strptime(t, "%H:%M:%S.%f") for t in o_alpha_tmpstmp]

            # X-AXIS Shared between input and output catalogs
            epoch_seconds = []
            for epoch_ in tmp_epoch:
                frac, whole = modf(epoch_)
                seconds_ = frac * 365.25 * 24 * 60
                epoch_seconds.append(seconds_)

            x_difference = float(max(epoch_seconds)) - float(
                min(epoch_seconds))
            x_major_step = (x_difference / 4)
            x_minor_step = (x_difference / 8)

            # X-SCALE
            x_major_ticks = [epoch_seconds[0] - x_major_step,
                             epoch_seconds[0],
                             epoch_seconds[0] + x_major_step,
                             epoch_seconds[0] + x_major_step * 2,
                             epoch_seconds[0] + x_major_step * 3,
                             epoch_seconds[0] + x_major_step * 4]

            x_minor_ticks = [epoch_seconds[0] - x_minor_step * 2,
                             epoch_seconds[0] - x_minor_step,
                             epoch_seconds[0],
                             epoch_seconds[0] + x_minor_step,
                             epoch_seconds[0] + x_minor_step * 2,
                             epoch_seconds[0] + x_minor_step * 3,
                             epoch_seconds[0] + x_minor_step * 4,
                             epoch_seconds[0] + x_minor_step * 5,
                             epoch_seconds[0] + x_minor_step * 6,
                             epoch_seconds[0] + x_minor_step * 7,
                             epoch_seconds[0] + x_minor_step * 8]

            ax.set_xlabel('EPOCH')

            # Y-AXIS
            myformat = mdates.DateFormatter('%S.%f')
            ax.yaxis.set_major_formatter(myformat)
            ax.set_ylabel('Right ascension(")')

            # format elements in alpha_seconds list to floats
            o_alpha_seconds = [float(i) for i in o_alpha_seconds]
            o_alpha_seconds = [float("{0:.6f}".format(i)) for i in
                               o_alpha_seconds]

            # gets maximum and minimum value in seconds
            # print('minimo {}'.format(min(o_alpha_seconds)))
            # min_alpha = '%.2E' % Decimal(min(o_alpha_seconds))
            # max_alpha = '%.2E' % Decimal(max(o_alpha_seconds))

            # print('minimo format {}'.format(min_alpha))
            # print('min_alpha {}'.format(min_alpha))
            # print('max_alpha {}'.format(max_alpha))

            max_alpha = max(o_alpha_seconds)
            min_alpha = min(o_alpha_seconds)

            difference = float(max_alpha) - float(min_alpha)
            major_step = (difference / 4)
            major_step_t = self.significant_l(major_step)

            first_digit = '%.2E' % Decimal(major_step - major_step_t)
            print('difference {}'.format(difference))
            print('major_step {}'.format(major_step))
            print('major_step_t {}'.format(major_step_t))
            print('first_digit {} - type {}'.format(first_digit,
                                                    type(first_digit)))
            if first_digit[0] == '-':
                first_digit = int(first_digit[1])
            else:
                first_digit = int(first_digit[0])

            # Redondea al alza el ultimo digio
            if first_digit > 5:
                last_digit = str(major_step_t)[-1]
                last_digit = int(last_digit)
                last_digit += 1
                major_step_t = list(str(major_step_t))
                major_step_t[-1] = last_digit
                major_step_t = [str(i) for i in major_step_t]
                major_step_t = ''.join(major_step_t)
                major_step_t = float(major_step_t)
            else:
                pass  # do nothing

            print('major_step {}'.format(major_step))
            print(' ')

            major_step = major_step_t
            minor_step = (difference / 16)

            min_datetime = min(o_datetimes)

            # Y-SCALE
            y_major_ticks = [min_datetime - timedelta(0, major_step),
                             min_datetime,
                             min_datetime + timedelta(0, major_step),
                             min_datetime + timedelta(0, major_step * 2),
                             min_datetime + timedelta(0, major_step * 3),
                             min_datetime + timedelta(0, major_step * 4)]

            y_minor_ticks = [min_datetime - timedelta(0, minor_step * 2),
                             min_datetime - timedelta(0, minor_step * 1),
                             min_datetime,
                             min_datetime + timedelta(0, minor_step),
                             min_datetime + timedelta(0, minor_step * 2),
                             min_datetime + timedelta(0, minor_step * 3),
                             min_datetime + timedelta(0, minor_step * 4),
                             min_datetime + timedelta(0, minor_step * 5),
                             min_datetime + timedelta(0, minor_step * 6),
                             min_datetime + timedelta(0, minor_step * 7),
                             min_datetime + timedelta(0, minor_step * 8),
                             min_datetime + timedelta(0, minor_step * 9),
                             min_datetime + timedelta(0, minor_step * 10),
                             min_datetime + timedelta(0, minor_step * 11),
                             min_datetime + timedelta(0, minor_step * 12),
                             min_datetime + timedelta(0, minor_step * 13),
                             min_datetime + timedelta(0, minor_step * 14),
                             min_datetime + timedelta(0, minor_step * 15),
                             min_datetime + timedelta(0, minor_step * 16),
                             min_datetime + timedelta(0, minor_step * 17),
                             min_datetime + timedelta(0, minor_step * 18)]

            y_minor_ticks_l = [round(min(o_alpha_seconds) - minor_step * 2, 4),
                               round(min(o_alpha_seconds) - minor_step, 4),
                               round(min(o_alpha_seconds), 4),
                               round(min(o_alpha_seconds) + minor_step, 4),
                               round(min(o_alpha_seconds) + minor_step * 2, 4),
                               round(min(o_alpha_seconds) + minor_step * 3, 4),
                               round(min(o_alpha_seconds) + minor_step * 4, 4),
                               round(min(o_alpha_seconds) + minor_step * 5, 4),
                               min(o_alpha_seconds) + minor_step * 6,
                               min(o_alpha_seconds) + minor_step * 7,
                               min(o_alpha_seconds) + minor_step * 8,
                               min(o_alpha_seconds) + minor_step * 9,
                               min(o_alpha_seconds) + minor_step * 10,
                               min(o_alpha_seconds) + minor_step * 11,
                               min(o_alpha_seconds) + minor_step * 12,
                               min(o_alpha_seconds) + minor_step * 13,
                               min(o_alpha_seconds) + minor_step * 14,
                               min(o_alpha_seconds) + minor_step * 15,
                               min(o_alpha_seconds) + minor_step * 16]

            # format elements in y_minor_tick_l for improved visualization
            y_minor_ticks_l = [float("{0:.6f}".format(i)) for i in
                               y_minor_ticks_l]

            print(arange(min(o_alpha_seconds), max(o_alpha_seconds),
                         major_step))

            # o_datetimes = [datetime.strptime(t, "%H:%M:%S.%f") for t in
            #                o_alpha_tmpstmp]

            for datetime_ in o_datetimes:
                print(datetime_)
            for tmpstmp_ in o_alpha_tmpstmp:
                print(tmpstmp_)

            ax.set_xticks(x_major_ticks, minor=False)
            ax.set_xticks(x_minor_ticks, minor=True)
            ax.set_yticks(arange(min(o_alpha_seconds),
                                 max(o_alpha_seconds),
                                 major_step), minor=False)

            ax.set_yticks(y_major_ticks, minor=False)
            ax.set_yticks(y_minor_ticks, minor=True)
            empty_string_labels = [''] * len(y_major_ticks)
            ax.set_yticklabels(empty_string_labels, minor=False)
            y_minor_ticks_l = [str(i) for i in y_minor_ticks_l]
            ax.set_yticklabels(y_minor_ticks_l, minor=True)

            ax.grid(b=True, which='major', linestyle='-', linewidth=2)
            ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

            for idx_datetime_, datetime_ in enumerate(o_datetimes):
                hour = datetime_.hour
                minute = datetime_.minute
                second = datetime_.second
                msecond = datetime_.microsecond
                alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
                                                          second, msecond)
                error = float(
                    "{0:.6f}".format(tmp_a_error[idx_datetime_] * 3600))
                error_str = '   error  {}"'.format(error)
                # Annotate position and error associated
                ax.annotate('{}\n{}'.format(alpha_str, error_str),
                            xy=(epoch_seconds[idx_datetime_], datetime_),
                            textcoords='data', fontsize=13)

            for idx_datetime_, datetime_ in enumerate(o_datetimes):
                ax.errorbar(epoch_seconds[idx_datetime_], datetime_,
                            yerr=timedelta(0, tmp_a_error[
                                idx_datetime_] * 3600),
                            fmt='o', ecolor='g', capthick=2)

            ax.plot(epoch_seconds, o_datetimes, 'bs', markersize=5,
                    label='extracted position')
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            plt.legend(loc=0, ncol=2, borderaxespad=0.)
            plt.draw()

            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure

            #
            #
            # DELTA PARAMETERS
            #
            #
            fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
            ax = fig.add_subplot(1, 1, 1)
            ax.set_title('chi_squared: {}'.format(fitted_d['dec']))

            o_delta_tmpstmp = []
            o_delta_seconds = []
            for delta_ in tmp_o_delta:
                d = Angle(delta_, u.degree)
                dms = d.dms
                o_delta_tmpstmp.append('{}:{}.{}'.format(int(dms[0]),
                                                         int(dms[1]), float("{0:.6f}".format(dms[2]))))
                o_delta_seconds.append('{}'.format(float("{0:.6f}".format(dms[2]))))

            ax.set_xlabel('EPOCH')
            ax.set_ylabel('Declination (")')

            # X-SCALE same as previous case
            ax.set_xticks(x_major_ticks, minor=False)
            ax.set_xticks(x_minor_ticks, minor=True)

            # Y-SCALE
            # format elements in alpha_seconds list to floats
            o_delta_seconds = [float(i) for i in o_delta_seconds]
            o_delta_seconds = [float("{0:.6f}".format(i)) for i in o_delta_seconds]

            max_delta = max(o_delta_seconds)
            min_delta = min(o_delta_seconds)

            delta_difference = float(max_delta) - float(min_delta)
            delta_major_step = (delta_difference / 4)
            delta_minor_step = (delta_difference / 16)

            # Y-SCALE
            y_major_ticks = [min_delta - delta_major_step,
                             min_delta,
                             min_delta + delta_major_step,
                             min_delta + delta_major_step * 2,
                             min_delta + delta_major_step * 3,
                             min_delta + delta_major_step * 4]

            y_minor_ticks = [min_delta - delta_minor_step * 2,
                             min_delta - delta_minor_step,
                             min_delta,
                             min_delta + delta_minor_step,
                             min_delta + delta_minor_step * 2,
                             min_delta + delta_minor_step * 3,
                             min_delta + delta_minor_step * 4,
                             min_delta + delta_minor_step * 5,
                             min_delta + delta_minor_step * 6,
                             min_delta + delta_minor_step * 7,
                             min_delta + delta_minor_step * 8,
                             min_delta + delta_minor_step * 9,
                             min_delta + delta_minor_step * 10,
                             min_delta + delta_minor_step * 11,
                             min_delta + delta_minor_step * 12,
                             min_delta + delta_minor_step * 13,
                             min_delta + delta_minor_step * 14,
                             min_delta + delta_minor_step * 15,
                             min_delta + delta_minor_step * 16]

            # rounds float to six decimal position
            y_minor_ticks = [float("{0:.6f}".format(i)) for i in
                             y_minor_ticks]

            # ax.set_yticks(y_major_ticks, minor=False)

            start, end = ax.get_ylim()
            # print('start {}'.format(start))
            # print('end {}'.format(end))
            # print(arange(min_delta, max_delta, 0.001))
            # print(delta_major_step)
            ax.set_yticks(arange(min_delta, max_delta,
                                 delta_major_step), minor=False)

            # ax.set_yticks(y_minor_ticks, minor=True)
            # empty_string_labels = [''] * len(y_major_ticks)
            # ax.set_yticklabels(empty_string_labels, minor=False)
            # converts floats into strings
            # y_minor_ticks = [str(i) for i in y_minor_ticks]
            # ax.set_yticklabels(y_minor_ticks, minor=True)
            ax.grid(b=True, which='major', linestyle='-', linewidth=2)
            ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

            for idx_datetime_, datetime_ in enumerate(o_datetimes):
                x = epoch_seconds[idx_datetime_]  # x position
                y = o_delta_seconds[idx_datetime_]  # y position
                # variables for position and error representation
                delta_position = o_delta_tmpstmp[idx_datetime_]
                delta_str = '   delta {}'.format(delta_position)
                error = tmp_b_error[idx_datetime_] * 3600  # error on y
                error_fmt = float("{0:.6f}".format(error))  # error rounded
                error_str = '   error {}'.format(error_fmt)
                # annotation
                ax.annotate('{}\n{}'.format(delta_str, error_str),
                            xy=(x, y), textcoords='data', fontsize=13)

            for idx_datetime_, datetime_ in enumerate(o_datetimes):
                x = epoch_seconds[idx_datetime_]  # x position
                y = o_delta_seconds[idx_datetime_]  # y position
                error = float(tmp_b_error[idx_datetime_] * 3600)  # error on y
                ax.errorbar(x, y, yerr=error, fmt='o', ecolor='g', capthick=2)

            ax.plot(epoch_seconds, o_delta_seconds, 'bs', markersize=5,
                    label='extracted position')

            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            plt.legend(loc=0, ncol=2, borderaxespad=0.)
            plt.draw()

            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure

"""