#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.2 Now supports confidence intervals
- 0.3 Sextractor performance's added
- 0.4 Proper motion performance's added
- 0.5 Proper motion now creates an stats file
- 0.6 Plotting methods out to plots.py file

Todo:
    * Improve log messages

"""
from os import makedirs, path
from subprocess import Popen

from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import arange, array, median
from pandas import concat, read_csv, DataFrame, Series
import statsmodels.api as sm

from misc import all_same, extract_settings
from misc import create_folder, speeds_range
from plots import PlotConfidence, PlotBothConfidence
from pyds9 import DS9
from regions import Create_regions


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.6"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def create_regions(alpha_l, delta_l, idx_):
    """

    :param alpha_l:
    :param delta_l:
    :param idx_:
    :return:
    """
    i_dict = {'alpha_l': alpha_l, 'delta_l': delta_l}
    i_df = DataFrame(i_dict)
    i_df.to_csv('{}.csv'.format(idx_), index=False,
                header=False, sep=" ")


def create_dict(scmp_cf, sex_cf, confidence_):
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


def cut_catalog(o_cat, margin, limits):
    """

    :param o_cat:
    :param margin:
    :param limits:
    :return:
    """
    o_df = o_cat[o_cat['ALPHA_J2000'] + margin > limits['max_alpha']]
    o_df = o_df[limits['min_alpha'] > o_df['ALPHA_J2000'] - margin]
    o_df = o_df[o_df['DELTA_J2000'] + margin > limits['max_delta']]
    o_df = o_df[limits['min_delta'] > o_df['DELTA_J2000'] - margin]

    return o_df


def check_source(catalog_n, o_cat, i_alpha, i_delta):
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


def check_star(catalog_n, i_df, o_alpha, o_delta):
    """

    :param catalog_n:
    :param i_df:0
    :param o_alpha:
    :param o_delta:
    :return:
    """
    tolerance = 0.0001

    i_df = i_df[i_df['catalog'].isin([catalog_n])]
    i_df = i_df[i_df['alpha_j2000'] + tolerance > o_alpha]
    i_df = i_df[o_alpha > i_df['alpha_j2000'] - tolerance]
    i_df = i_df[i_df['delta_j2000'] + tolerance > o_delta]
    i_df = i_df[o_delta > i_df['delta_j2000'] - tolerance]

    return i_df


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


def check_cat_order(cat_list):
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


class StatsPerformance:
    def __init__(self, prfs_d):
        """

        :param prfs_d:
        """
        self.prfs_d = prfs_d

    def error(self, logger, mag, sex_cf):
        """

        :param logger:
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
            cat_name = '{}/Cat_20-21_d{}'.format(self.prfs_d['input_ref'], d)
            input_ssos_d[d] = '{}.dat'.format(cat_name)
        input_ssos_d = Create_regions(input_ssos_d,
                                      self.prfs_d).check_luca(True, True)

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
        input_stars_d = Create_regions(input_stars_d,
                                       self.prfs_d).check_stars(True, True)

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
            input_galaxies_d[d] = '{}{}{}.dat'.format(self.prfs_d['input_ref'],
                                                      i_cat_n, d)
        input_galaxies_d = Create_regions(input_galaxies_d,
                                          self.prfs_d).check_galaxies(True,
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
                errors_a_stars = []
                errors_b_stars = []
                errors_a_galaxies = []
                errors_b_galaxies = []
                errors_a_ssos = []
                errors_b_ssos = []

                cat_n = 'mag_{}_CCD_x{}_y{}_d{}.cat'.format(mag, opt[0],
                                                            opt[1], dither)
                cat_o_n = '{}/{}/{}'.format(self.prfs_d['fits_dir'],
                                            sex_cf, cat_n)

                hdu_list = fits.open(cat_o_n)
                o_cat = Table(hdu_list[2].data).to_pandas()

                i_stars_d_df = i_stars_df[i_stars_df['dither_values'].isin([dither])]
                i_galaxies_d_df = i_galaxies_df[i_galaxies_df['dither_values'].isin([dither])]
                i_ssos_d_df = i_ssos_df[i_ssos_df['dither_values'].isin([dither])]
                for idx, row in o_cat.iterrows():
                    i_alpha = row['ALPHA_J2000']
                    i_delta = row['DELTA_J2000']
                    o_df = search(i_stars_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        errors_a_stars.append(row['ERRA_IMAGE'])
                        errors_b_stars.append(row['ERRB_IMAGE'])
                    o_df = search(i_galaxies_d_df, i_alpha, i_delta)
                    if o_df.empty is not True:
                        errors_a_galaxies.append(row['ERRA_IMAGE'])
                        errors_b_galaxies.append(row['ERRB_IMAGE'])
                    o_df = search(i_ssos_d_df, i_alpha, i_delta)
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
                o_df = check_source(catalog_n, o_cat, i_alpha, i_delta)

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
        :param bypassfilter:
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

        with PdfPages('{}/pm/{}/{}_{}_cmp.pdf'.format(prfs_d['plots_dir'],
                                                      sex_cf, scmp_cf,
                                                      confidence_)) as pdf:
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

            ax.grid(b=True, which='major', ls='-', lw=2)
            ax.grid(b=True, which='minor', ls='--', lw=1)
            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure

        return True


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
        input_ssos_d = Create_regions(input_ssos_d,
                                      self.prfs_d).check_luca(True, True)

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
        input_stars_d = Create_regions(input_stars_d,
                                       self.prfs_d).check_stars(True, True)

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
        input_galaxies_d = Create_regions(input_galaxies_d,
                                          self.prfs_d).check_galaxies(True,
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

                i_stars_d_df = i_stars_df[i_stars_df['dither_values'].isin([dither])]
                i_galaxies_d_df = i_galaxies_df[i_galaxies_df['dither_values'].isin([dither])]
                i_ssos_d_df = i_ssos_df[i_ssos_df['dither_values'].isin([dither])]
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


class ScampPerformanceSSOs:

    def __init__(self, logger, mag, scmp_cf, sex_cf, confidence_):
        """

        """
        self.bypassfilter = True
        self.prfs_d = extract_settings()

        self.mag = mag
        self.scmp_cf = scmp_cf
        self.sex_cf = sex_cf
        self.confidence_ = confidence_

        self.check(logger)
        pass

    def check(self, logger):
        """

        :param logger:
        :return:
        """
        # For now any file is saved
        save = True

        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {} and {}'.format(self.scmp_cf,
                                                                 self.sex_cf))

        input_d = {}
        for d in range(1, 5, 1):
            cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'],
                                              self.mag)
            cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
            input_d[d] = '{}.dat'.format(cat_name)
        input_d = Create_regions(input_d, self.prfs_d).check_luca(self.mag,
                                                                  False, True)

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
        filt_n = 'filt_{}_{}_4.csv'.format(self.scmp_cf, self.mag)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filt_n)

        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        stats_d, out_d = create_dict(self.scmp_cf, self.sex_cf,
                                     self.confidence_)

        # Gets unique sources from input data
        unique_sources = list(set(i_df['source'].tolist()))

        # Loops over input data
        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            cat = i_df[i_df['source'].isin([source_])]
            # Creates lists for each source
            tmp_d = {'boolean_l': [], 'catalog': [], 'source': [],
                     'epoch': [], 'i_pm': [],
                     'i_alpha': [], 'i_delta': [],
                     'pm_alpha': [], 'pm_delta': [],
                     'o_alpha': [], 'o_delta': [],
                     'error_a': [], 'error_b': []}

            i_pm = 0.0  # Raises a warning if variable pm it's not created
            # Iterate over each detection of each source
            for i, row in enumerate(cat.itertuples(), 1):
                # source_ = row.source
                # ccd_ = row.CCD
                # dither_ = row.dither_values
                catalog_n = row.catalog
                i_pm = row.pm_values
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                print('catalog_n {}'.format(catalog_n))
                # Checks if there is a source closed to input one
                # o_cat contains data from output (filtered) catalog
                o_df = check_source(catalog_n, o_cat, i_alpha, i_delta)

                # If there is one saves data from input data
                if o_df.empty is not True and o_df['PM'].size == 1:
                    pm_mask = self.pm_filter(o_df, i_pm, self.confidence_)
                    if pm_mask:
                        if o_df['SOURCE_NUMBER'].size != 1:
                            tmp_d['boolean_l'].append('False')
                        else:
                            tmp_d['boolean_l'].append('True')
                        tmp_d['catalog'].append(catalog_n)
                        source = int(o_df['SOURCE_NUMBER'].iloc[0])
                        tmp_d['source'].append(source)
                        print(o_df['EPOCH'])
                        epoch = float(o_df['EPOCH'].iloc[0])
                        tmp_d['epoch'].append(epoch)
                        tmp_d['i_pm'].append(i_pm)
                        tmp_d['i_alpha'].append(i_alpha)
                        tmp_d['i_delta'].append(i_delta)
                        pm_alpha = float(o_df['PMALPHA'])
                        tmp_d['pm_alpha'].append(pm_alpha)
                        pm_delta = float(o_df['PMDELTA'])
                        tmp_d['pm_delta'].append(pm_delta)
                        o_alpha = float(o_df['ALPHA_J2000'].iloc[0])
                        tmp_d['o_alpha'].append(o_alpha)
                        o_delta = float(o_df['DELTA_J2000'].iloc[0])
                        tmp_d['o_delta'].append(o_delta)
                        errora_world = float(o_df['ERRA_WORLD'].iloc[0])
                        tmp_d['error_a'].append(errora_world)
                        errorb_world = float(o_df['ERRB_WORLD'].iloc[0])
                        tmp_d['error_b'].append(errorb_world)
                else:
                    tmp_d['boolean_l'].append('False')

            # Total number
            idx = stats_d['PM'].index(i_pm)
            stats_d['total'][idx] += 1

            order_mask = check_cat_order(tmp_d['catalog'])

            if not order_mask:
                raise Exception

            if len(tmp_d['source']) is not 0:
                flag_detection, sources_number = all_same(tmp_d['source'])

                if flag_detection and sources_number >= 3:
                    print('epoch {}'.format(tmp_d['epoch']))
                    idx = stats_d['PM'].index(i_pm)
                    stats_d['right'][idx] += 1
                    fitted_d = self.confidence(tmp_d['source'][0],
                                               self.scmp_cf, self.sex_cf)

                    pm = float(tmp_d['i_pm'][0])

                    # Output folder creation
                    plots_dir = self.prfs_d['plots_dir']
                    self.output_path = '{}/mv/{}/{}/{}'.format(plots_dir,
                                                               self.scmp_cf,
                                                               self.sex_cf,
                                                               pm)
                    create_folder(logger, self.output_path)

                    self.plot(tmp_d, pm, fitted_d)
                    # input_fits = 'test'
                    # regions_d = self.create_regions(tmp_d['source'][0],
                    #                                 o_cat, tmp_d)
                    # self.image(input_fits)
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

    def plot(self, tmp_d, pm, fitted_d):
        """

        :param tmp_d:
        :param pm:
        :param fitted_d:
        :return:
        """
        # Set True to plot input and output data
        both = True

        if both:
            mode = 'io'
            d_ = {'i_alpha': tmp_d['i_alpha'], 'i_delta': tmp_d['i_delta'],
                  'pm_alpha': tmp_d['pm_alpha'], 'pm_delta': tmp_d['pm_delta'],
                  'o_alpha': tmp_d['o_alpha'], 'o_delta': tmp_d['o_delta'],
                  'error_a': tmp_d['error_a'], 'error_b': tmp_d['error_b'],
                  'epoch': tmp_d['epoch'], 'i_pm': tmp_d['i_pm']}
            plot = PlotBothConfidence(self.output_path, tmp_d['source'][0],
                                      pm, mode, fitted_d, d_)
            if not plot:
                raise Exception
        else:
            mode = 'o'
            dict_ = {'alpha': tmp_d['o_alpha'], 'delta': tmp_d['o_delta'],
                     'error_a': tmp_d['error_a'], 'error_b': tmp_d['error_b'],
                     'epoch': tmp_d['epoch']}
            plot = PlotConfidence(self.output_path, tmp_d['source'][0], pm,
                                  mode, fitted_d, dict_)

            if not plot:
                raise Exception

    def confidence(self, source, scmp_cf, sex_cf):
        """

        :param source:
        :param scmp_cf:
        :param sex_cf:
        :return:
        """
        catalogs_dir = self.prfs_d['catalogs_dir']
        configuration_dir = '/{}/{}/{}'.format(self.mag, sex_cf, scmp_cf)
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
                edim = array([1 / var for var in sigma])
                edims.append(edim)
                x = sm.add_constant(x)
                # model = sm.WLS(y, x, weigths=edim)
                model = sm.WLS(y, x)
                fitted = model.fit()
                fitted = fitted.rsquared
                fitted_d[dimension_] = fitted
            elif dimension == dec:
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

    def image(self, input_fits):
        """

        :return:
        """
        d = DS9()
        # Open file
        open_file = 'file {}'.format(input_fits)
        d.set(open_file)
        scale = 'scale histequ'
        d.set(scale)
        load_regions = 'regions load'
        d.set(load_regions)
        cmd_12 = '-regions -format xy '
        # load inputs regions
        cmd_13 = '{}{}.reg '.format(self.output_path)
        # cmd_14 = '-saveimage jpeg {} -exit'.format(output_image)

        return True

    def create_regions(self, source, o_cat, tmp_d):
        """

        :param source:
        :param o_cat:
        :param tmp_d:
        :return:
        """
        # En azul! Valores totales
        # Gets full list of alpha/delta coordinates
        alpha_l = tmp_d['i_alpha'] + tmp_d['o_alpha']
        delta_l = tmp_d['i_delta'] + tmp_d['o_delta']

        d_limits = {'max_alpha': float(max(alpha_l)),
                    'min_alpha': float(min(alpha_l)),
                    'max_delta': float(max(delta_l)),
                    'min_delta': float(min(delta_l))}

        margin = 0.001
        o_df = cut_catalog(o_cat, margin, d_limits)

        # o_df.to_csv('{}.reg'.format(source), index=False,
        #             header=False, sep=" ")
        # Tests reasons
        f_regions_filename = 'f_{}.reg'.format(source)
        o_df.to_csv(f_regions_filename, index=False, sep=" ")

        # En rojo! Input values from Luca's
        alpha_list = Series(tmp_d['i_alpha'], name='ALPHA_J2000')
        delta_list = Series(tmp_d['i_delta'], name='DELTA_J2000')

        i_df = concat([alpha_list, delta_list], axis=1)
        # i_df.to_csv('i_{}.reg'.format(source), index=False,
        #             header=False, sep=" ")
        i_regions_filename = 'i_{}.reg'.format(source)
        i_df.to_csv(i_regions_filename, index=False, sep=" ")

        d_regions_filename = {'f_filename': f_regions_filename,
                              'i_filename': i_regions_filename}

        return d_regions_filename


class ScampPerformanceStars:

    def __init__(self):
        """

        """
        self.bypassfilter = True
        self.prfs_d = extract_settings()

    def check(self, logger, mag, scmp_cf, sex_cf):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_cf:
        :return:
        """
        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {} and {}'.format(scmp_cf,
                                                                 sex_cf))
        input_d = {}
        for d in range(1, 5, 1):
            cat_name = '{}/Cat_20-21_d{}'.format(self.prfs_d['input_ref'], d)
            input_d[d] = '{}.dat'.format(cat_name)
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

        # Open particular file!
        filt_n = 'filt_{}_{}_3.csv'.format(scmp_cf, mag)
        filter_o_n = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                          sex_cf, scmp_cf, filt_n)

        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        # Gets unique sources from input data
        o_unique_sources = list(set(o_cat['SOURCE_NUMBER'].tolist()))

        # Creates variables
        right = 0
        not_right = 0
        total = 0
        stats_d = {'source_l': [], 'catalog_l': [], 'alpha_l': [],
                   'delta_l': [], 'pm_l': [], 'mag': [], 'error_a': [],
                   'error_b': []}

        for source_ in o_unique_sources:
            flag = False
            tmp_d = {'source': [], 'catalog': [], 'i_pm': [],
                     'i_alpha': [], 'i_delta': [], 'mag': [],
                     'error_a': [], 'error_b': []}

            o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                source = row.SOURCE_NUMBER
                catalog_n = row.CATALOG_NUMBER
                i_alpha = row.ALPHA_J2000
                i_delta = row.DELTA_J2000
                pm = row.PM
                mag = row.MAG
                error_a = row.ERRA_WORLD
                error_b = row.ERRB_WORLD

                out_df = check_star(catalog_n, i_df, i_alpha, i_delta)

                if out_df.empty:
                    # it's a star
                    right += 1
                    flag = True
                    tmp_d['source'].append(source)
                    tmp_d['catalog'].append(catalog_n)
                    tmp_d['i_alpha'].append(i_alpha)
                    tmp_d['i_delta'].append(i_delta)
                    tmp_d['i_pm'].append(pm)
                    tmp_d['mag'].append(mag)
                    tmp_d['error_a'].append(error_a)
                    tmp_d['error_b'].append(error_b)
                else:
                    # it's a SSO
                    not_right += 1
                    flag = False

                total += 1

            if flag:
                order_mask = check_cat_order(tmp_d['catalog'])
                flag_detection, sources_number = all_same(tmp_d['source'])

                if order_mask and sources_number >= 3:
                    for tmp_source_ in tmp_d['source']:
                        stats_d['source_l'].append(tmp_source_)
                    for catalog_ in tmp_d['catalog']:
                        stats_d['catalog_l'].append(catalog_)
                    for i_alpha_ in tmp_d['i_alpha']:
                        stats_d['alpha_l'].append(i_alpha_)
                    for i_delta_ in tmp_d['i_delta']:
                        stats_d['delta_l'].append(i_delta_)
                    for pm_ in tmp_d['i_pm']:
                        stats_d['pm_l'].append(pm_)
                    for mag_ in tmp_d['mag']:
                        stats_d['mag'].append(mag_)
                    for erra_ in tmp_d['error_a']:
                        stats_d['erra_world'].append(erra_)
                    for errb_ in tmp_d['error_b']:
                        stats_d['errb_world'].append(errb_)

        stats_df = DataFrame(stats_d)
        stats_df.to_csv('test_{}_{}.csv'.format(sex_cf, scmp_cf))

        return stats_d
