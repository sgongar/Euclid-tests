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

from astropy.io import fits
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import arange, array, median
from pandas import concat, read_csv, DataFrame

from misc import extract_settings
from misc import speeds_range
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
        self.mag = mag
        errors_a_stars = []
        errors_b_stars = []
        errors_a_galaxies = []
        errors_b_galaxies = []
        errors_a_ssos = []
        errors_b_ssos = []

        # Input sources
        input_ssos_d = {}
        for d in range(1, 5, 1):
            cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], mag)
            cat_name = '{}/Cat_{}_d{}'.format(cat_loc, self.mag, d)
            input_ssos_d[d] = '{}.dat'.format(cat_name)
        input_ssos_d = Create_regions(input_ssos_d).check_luca(self.mag, False, True)

        # Creates a DataFrame from an input dictionary
        input_ssos_l = []
        for key_ in input_ssos_d.keys():
            input_ssos_l.append(input_ssos_d[key_])

        i_ssos_df = concat(input_ssos_l, axis=0)
        i_ssos_df = i_ssos_df.reset_index(drop=True)

        # Input stars
        input_stars_d = {}
        for d in range(1, 5, 1):
            cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], mag)
            cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
            input_stars_d[d] = '{}.dat'.format(cat_name)
        input_stars_d = Create_regions(input_stars_d).check_luca(mag, False, True)

        # Creates a DataFrame from an input dictionary
        input_stars_l = []
        for key_ in input_stars_d.keys():
            input_stars_l.append(input_stars_d[key_])

        i_stars_df = concat(input_stars_l, axis=0)
        i_stars_df = i_stars_df.reset_index(drop=True)

        # Input stars
        input_galaxies_d = {}
        for d in range(1, 5, 1):
            cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], mag)
            cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
            input_galaxies_d[d] = '{}.dat'.format(cat_name)
        input_galaxies_d = Create_regions(input_galaxies_d).check_luca(mag, False,
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
                cat_o_n = '{}/{}/CCDs/{}/{}'.format(self.prfs_d['fits_dir'],
                                                    mag, sex_cf, cat_n)

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
        self.prfs_d = extract_settings()

        pass

    def check(self, logger, mag, scmp_cf, sex_cf, confidence_):
        """

        :param logger:
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
            input_d[d] = '{}/Cat_20-21_d{}.dat'.format(self.prfs_d['input_ref'],
                                                       d)
        input_d = Create_regions(input_d).check_luca(mag, True, True)

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
        filt_n = 'filt_{}_{}_4.csv'.format(scmp_cf, mag)
        filter_o_n = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
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
                    pm_mask = self.pm_filter(o_df, i_pm, confidence_)
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
            if not self.plot_comp(sex_cf, scmp_cf, input_pm, output_pm,
                                  confidence_):
                raise Exception

        return stats_d

    def pm_filter(self, o_df, pm, confidence_):
        """

        :param o_df:
        :param pm:
        :param prfs_d:
        :param confidence_:
        :param bypassfilter:
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

    def plot_comp(self, sex_cf, scmp_cf, input_pm, output_pm, confidence_):
        """ fixme support for multiple mags

        :param sex_cf:
        :param scmp_cf:
        :param input_pm:
        :param output_pm:
        :param confidence_:
        :return:
        """
        sextractor_folder = '{}/pm/{}'.format(self.prfs_d['plots_dir'], sex_cf)
        if not path.exists(sextractor_folder):
            makedirs(sextractor_folder)

        with PdfPages('{}/pm/{}/{}_{}_cmp.pdf'.format(self.prfs_d['plots_dir'],
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
