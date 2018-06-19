# !/usr/bin/python
# -*- coding: utf-8 -*-

""" Gets
   - factors from magnitude bins

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
from sys import argv

from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import array, isnan, nan, nanmean
from pandas import concat, DataFrame, read_csv, Series

from misc import extract_settings_luca
from regions import Create_regions


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def compute_factors(factors_d, stats_df):
    """
    N_meas: number of all detected sources(including false detections)
    N_se: number of simulated sources recovered by source extraction
    N_true: number of simulated input sources
    f_dr: detection rate f_pur: purity
    f_com: completeness

    f_dr = N_meas / N_true = (N_se + N_false) / N_true
    f_pur = N_se / N_meas = N_se / (N_se + N_false)
    f_com = f_dr * f_pur = N_se / N_true

    :param stats_d:
    :param tmp_d:
    :return:
    """
    prfs_d = extract_settings_luca()

    mags = [[14, 15], [15, 16], [16, 17], [17, 18], [18, 19],
            [19, 20], [20, 21], [21, 22], [22, 23], [23, 24],
            [24, 25], [25, 26], [26, 27], [27, 28]]
    pms = ['right-0', 'false-0', 'total-0',
           'right-0.1', 'false-0.1', 'total-0.1',
           'right-0.3', 'false-0.3', 'total-0.3',
           'right-1.0', 'false-1.0', 'total-1.0',
           'right-3.0', 'false-3.0', 'total-3.0',
           'right-10.0', 'false-10.0', 'total-10.0']
    pm_list = [0.1, 0.3, 1.0, 3.0, 10.0]
    idxs_pm = [[4, 5, 6], [7, 8, 9], [10, 11, 12],
               [13, 14, 15], [16, 17, 18]]
    for idx_mag, data in enumerate(stats_df.itertuples(), 1):
        # Loop over the different magnitudes
        if 6 < idx_mag < 13:  # Only gets values from magnitudes [20, 26]
            mag_ = data[0]  # Gets magnitude value
            for idx_pm, pm_values in enumerate(idxs_pm): # each sublist is a different pm
                pm_ = pm_list[idx_pm]
                n_se = data[pm_values[0]]
                factors_d[mag_][pm_]['n_se'] = n_se
                n_false = data[pm_values[1]]
                factors_d[mag_][pm_]['n_false'] = n_false
                n_meas = n_se + n_false
                factors_d[mag_][pm_]['n_meas'] = n_meas
                n_true = data[pm_values[2]]
                factors_d[mag_][pm_]['n_true'] = n_true

                try:
                    f_dr = float(n_meas) / float(n_true)
                    f_dr = float("{0:.2f}".format(f_dr))
                    factors_d[mag_][pm_]['f_dr'] = f_dr
                except ZeroDivisionError:
                    factors_d[mag_][pm_]['f_dr'] = nan
                try:
                    f_pur = float(n_se) / float(n_meas)
                    f_pur = float("{0:.2f}".format(f_pur))
                    factors_d[mag_][pm_]['f_pur'] = f_pur
                except ZeroDivisionError:
                    factors_d[mag_][pm_]['f_pur'] = nan
                try:
                    f_com = float(n_se) / float(n_true)
                    f_com = float("{0:.2f}".format(f_com))
                    factors_d[mag_][pm_]['f_com'] = f_com
                except ZeroDivisionError:
                    factors_d[mag_][pm_]['f_com'] = nan
        else:
            pass

    print(factors_d)

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

    o_df = o_df[o_df[keys[0]] + prfs_d['tolerance']*4 > i_alpha]
    o_df = o_df[i_alpha > o_df[keys[0]] - prfs_d['tolerance']*4]
    o_df = o_df[o_df[keys[1]] + prfs_d['tolerance']*4 > i_delta]
    o_df = o_df[i_delta > o_df[keys[1]] - prfs_d['tolerance']*4]

    return o_df


def speeds_range(confidence):
    """ given a confidence value returns a dict with speeds

    @param confidence:

    @return speeds_dict:
    """
    prfs_d = extract_settings_elvis()
    pms = prfs_d['pms']
    speeds_dict = {}
    for pm_ in pms:
        speeds_dict[pm_] = [pm_ - pm_ * confidence / 100.0,
                            pm_ + pm_ * confidence / 100.0]

    return speeds_dict


def get_norm_speed(o_pm):
    """

    :return:
    """
    speeds_d = speeds_range(confidence=50)

    pm_norm = 0
    for key_ in speeds_d.keys():
        low = speeds_d[key_][0]
        high = speeds_d[key_][1]
        if low < o_pm < high:
            pm_norm = key_

    return pm_norm


def get_norm_mag(o_mag):
    """

    :param o_mag:
    :return: mag_bin
    """
    mags = [[14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [19, 20],
            [20, 21], [21, 22], [22, 23], [23, 24], [24, 25], [25, 26],
            [26, 27], [27, 28]]

    mag_bin = ''
    for mag_ in mags:
        if mag_[0] < o_mag < mag_[1]:
            mag_bin = '{}-{}'.format(mag_[0], mag_[1])

    return mag_bin


def splits_by_mag_bin():
    """

    :return: total_d
    """
    total_d = {}
    # opens input_catalog
    cat_ssos = read_csv('cats/cat_clean_ssos.csv', index_col=0)
    test = True

    mags = [[20, 21], [21, 22], [22, 23], [23, 24],
            [24, 25], [25, 26], [26, 27], [27, 28]]
    pms = [0.1, 0.3, 1.0, 3.0, 10.0]
    pms_range = speeds_range(50)

    total_d['14-15'] = {}
    total_d['14-15'][0] = 0
    for pm_ in pms:
        total_d['14-15'][pm_] = 0

    total_d['15-16'] = {}
    total_d['15-16'][0] = 0
    for pm_ in pms:
        total_d['15-16'][pm_] = 0

    total_d['16-17'] = {}
    total_d['16-17'][0] = 0
    for pm_ in pms:
        total_d['16-17'][pm_] = 0

    total_d['17-18'] = {}
    total_d['17-18'][0] = 0
    for pm_ in pms:
        total_d['17-18'][pm_] = 0

    total_d['18-19'] = {}
    total_d['18-19'][0] = 0
    for pm_ in pms:
        total_d['18-19'][pm_] = 0

    total_d['19-20'] = {}
    total_d['19-20'][0] = 0
    for pm_ in pms:
        total_d['19-20'][pm_] = 0

    for mag_ in mags:
        mag_bin_cat = cat_ssos[cat_ssos['ABMAG'] < mag_[1]]
        mag_bin_cat = mag_bin_cat[mag_bin_cat['ABMAG'] > mag_[0]]
        total_d['{}-{}'.format(mag_[0], mag_[1])] = {}
        total_d['{}-{}'.format(mag_[0], mag_[1])][0] = 0
        for pm_ in pms:
            # gets proper motion bin
            pm_bin_cat = mag_bin_cat[mag_bin_cat['VEL'] < pms_range[pm_][1]]
            pm_bin_cat = pm_bin_cat[pm_bin_cat['VEL'] > pms_range[pm_][0]]
            try:
                pm_bin_cat = concat(g for _, g in pm_bin_cat.groupby('SOURCE')
                                    if len(g) >= 3)
                total_sources = len(list(set(pm_bin_cat['SOURCE'].tolist())))
                total_d['{}-{}'.format(mag_[0], mag_[1])][pm_] = total_sources
            except ValueError:
                total_d['{}-{}'.format(mag_[0], mag_[1])][pm_] = 0

    return total_d


def redo_data_d():
    """ Creates a dictionary
    TODO - Automatic number!

    :return: tmp_d
    """
    prfs_d = extract_settings_elvis()
    total_d = splits_by_mag_bin()

    data_d = {}
    mags = [[14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [19, 20],
            [20, 21], [21, 22], [22, 23], [23, 24], [24, 25], [25, 26],
            [26, 27], [27, 28]]
    pms = prfs_d['pms']
    for mag_ in mags:
        mag_bin = '{}-{}'.format(mag_[0], mag_[1])
        data_d[mag_bin] = {}
        data_d[mag_bin][0] = {'right': 0, 'false': 0,
                              'total': total_d[mag_bin][0]}
        for pm_ in pms:
            data_d[mag_bin][pm_] = {'right': 0, 'false': 0,
                                    'total': total_d[mag_bin][pm_]}

    return data_d


def redo_factors_d():
    """ Creates a dictionary
    TODO - Automatic number!

    :return: factors_d
    """
    prfs_d = extract_settings_elvis()

    factors_d = {}
    mags = [[20, 21], [21, 22], [22, 23], [23, 24],
            [24, 25], [25, 26], [26, 27]]
    pms = prfs_d['pms']
    for mag_ in mags:
        mag_bin = '{}-{}'.format(mag_[0], mag_[1])
        factors_d[mag_bin] = {}
        for pm_ in pms:
            factors_d[mag_bin][pm_] = {'f_pur': 0.0, 'f_dr': 0.0, 'f_com': 0.0,
                                       'n_se': 0, 'n_false': 0, 'n_meas': 0, 'n_true': 0}

    return factors_d


class FactorsScampPerformance:

    def __init__(self):
        """ Me da los valores de salida de todos las estrellas  y galaxias
        presentes en filt 3 obtenidos o no

        """
        self.filter_p_number = 9  # First one with enough data for statistics
        self.prfs_d = extract_settings_elvis()
        self.data_d = redo_data_d()
        factors_d = redo_factors_d()

        self.save = True

        logger_name = 'scamp_performance'  # Set as desired
        self.logger = setting_logger(self.prfs_d, logger_name)

        filt_cat = self.gets_filtered_catalog()  # Gets data from filtered
        input_df = self.gets_data()  # Gets data from catalogs
        stats_df = self.extract_stats(filt_cat, input_df)  # Splits due type
        stats_df.to_csv('stats.csv')
        compute_factors(factors_d, stats_df)

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
            # ssos_cat = 'cat_ssos_{}.csv'.format(key_)
            ssos_cat = 'cats/cat_ssos_{}.csv'.format(key_)
            input_df[key_]['SSOs'] = read_csv(ssos_cat, index_col=0)
            stars_cat = 'tmp_stars/stars.csv'
            input_df[key_]['stars'] = read_csv(stars_cat, index_col=0)
            galaxies_cat = 'tmp_galaxies/galaxies.csv'
            input_df[key_]['galaxies'] = read_csv(galaxies_cat, index_col=0)

        return input_df

    def extract_stats(self, filt_cat, input_df):
        """

        :param filt_cat:
        :param input_df:
        :return:
        """
        # Unique sources (?)
        prfs_d = extract_settings_elvis()
        unique_sources = list(set(filt_cat['SOURCE_NUMBER'].tolist()))

        print('Creating new catalogues from filtered catalogue due type')
        print('Total sources: {}'.format(filt_cat['SOURCE_NUMBER'].size))
        for idx_source_, source_ in enumerate(unique_sources):
            print(idx_source_)
            source_df = filt_cat[filt_cat['SOURCE_NUMBER'].isin([source_])]
            # Takes the first value of MAG Series
            i_mag_bin = get_norm_mag(source_df['MEDIAN_MAG_ISO'].iloc[0])
            # Takes the first value of PM Series
            i_pm_norm = get_norm_speed(source_df['PM'].iloc[0])
            alpha = source_df['ALPHA_J2000'].iloc[0]
            delta = source_df['DELTA_J2000'].iloc[0]

            source_d = {'source': [], 'pm': [], 'mag': []}
            right_detections = 0
            for i, row in enumerate(source_df.itertuples(), 1):
                dither_n = get_dither(int(row.CATALOG_NUMBER))
                # Checks object type
                keys = ['RA', 'DEC']  # Catalogue version 2
                test_sso = check_source(input_df[dither_n]['SSOs'],
                                        alpha, delta, keys)
                if test_sso.empty is not True:
                    right_detections += 1
                    source_d['source'].append(row.SOURCE_NUMBER)
                    source_d['pm'].append(row.PM)
                    source_d['mag'].append(row.MEDIAN_MAG_ISO)

            if right_detections >= 3:
                o_mag_bin = get_norm_mag(source_d['mag'][0])
                o_pm_norm = get_norm_speed(source_d['pm'][0])
                self.data_d[o_mag_bin][o_pm_norm]['right'] += 1
            else:
                # print(source_df['MEDIAN_MAG_ISO'].iloc[0])
                # print(i_mag_bin, i_pm_norm)
                # print(self.data_d[i_mag_bin][i_pm_norm])
                # print(' ')
                self.data_d[i_mag_bin][i_pm_norm]['false'] += 1

        mags = [[14, 15], [15, 16], [16, 17], [17, 18], [18, 19],
                [19, 20], [20, 21], [21, 22], [22, 23], [23, 24],
                [24, 25], [25, 26], [26, 27], [27, 28]]
        for mag_ in mags:
            mag_bin = '{}-{}'.format(mag_[0], mag_[1])
            data_df = DataFrame(self.data_d[mag_bin])
            data_df.to_csv('stats_{}.csv'.format(mag_bin))

        idxs = ['14-15', '15-16', '16-17', '17-18', '18-19',
                '19-20', '20-21', '21-22', '22-23', '23-24',
                '24-25', '25-26', '26-27', '27-28']
        pms = prfs_d['pms']
        stats_d = {'idx': idxs}
        for pm_ in pms:
            stats_d['right-{}'.format(pm_)] = []
            stats_d['false-{}'.format(pm_)] = []
            stats_d['total-{}'.format(pm_)] = []
            for mag_ in mags:
                mag_bin = '{}-{}'.format(mag_[0], mag_[1])
                stats_d['right-{}'.format(pm_)].append(self.data_d[mag_bin][pm_]['right'])
                stats_d['false-{}'.format(pm_)].append(self.data_d[mag_bin][pm_]['false'])
                stats_d['total-{}'.format(pm_)].append(self.data_d[mag_bin][pm_]['total'])

        stats_df = DataFrame(stats_d,
                             columns=['idx', 'right-0', 'false-0', 'total-0',
                                      'right-0.1', 'false-0.1', 'total-0.1',
                                      'right-0.3', 'false-0.3', 'total-0.3',
                                      'right-1.0', 'false-1.0', 'total-1.0',
                                      'right-3.0', 'false-3.0', 'total-3.0',
                                      'right-10.0', 'false-10.0', 'total-10.0'])
        stats_df = stats_df.set_index('idx')
        stats_df.to_csv('stats.csv')

        return stats_df


if __name__ == "__main__":
    FactorsScampPerformance()
