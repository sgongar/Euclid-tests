#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance
Mide las detecciones de cada tipo

Versions:
- 0.1

Todo:
    * Improve log messages

"""
from pandas import concat, read_csv, DataFrame
from numpy import mean, std

from misc import all_same, extract_settings
from misc import speeds_range
from regions import Create_regions


def compute_factors(stats_d):
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
    :return:
    """
    for idx in range(0, len(stats_d['f_dr']), 1):
        n_meas = stats_d['N_meas'][idx]
        n_se = stats_d['N_se'][idx]
        n_true = stats_d['N_true'][idx]
        try:
            f_dr = n_meas / n_true
            f_dr = float("{0:.2f}".format(f_dr))
            stats_d['f_dr'][idx] = f_dr
        except ZeroDivisionError:
            stats_d['f_dr'][idx] = 'nan'
        try:
            f_pur = n_se / n_meas
            f_pur = float("{0:.2f}".format(f_pur))
            stats_d['f_pur'][idx] = f_pur
        except ZeroDivisionError:
            stats_d['f_pur'][idx] = 'nan'
        try:
            f_com = n_se / n_true
            f_com = float("{0:.2f}".format(f_com))
            stats_d['f_com'][idx] = f_com
        except ZeroDivisionError:
            stats_d['f_com'][idx] = 'nan'

    return stats_d


def redo_stats_d():
    """ Creates a dictionary

    :return: tmp_d
    """
    stats_d = {'i_pm': [], 'N_meas': [], 'N_false': [], 'N_se': [],
               'N_true': [], 'f_dr': [], 'f_pur': [], 'f_com': [],
               'alpha_std': [], 'delta_std': []}

    return stats_d


def populate_stats_d(stats_d):
    """

    :param stats_d:
    :return:
    """
    prfs_d = extract_settings()
    initial_int = 0
    initial_float = 0.0

    stats_d['i_pm'].append(initial_int)
    stats_d['N_meas'].append(initial_int)
    stats_d['N_false'].append(initial_int)
    stats_d['N_se'].append(initial_float)
    stats_d['N_true'].append(initial_float)
    stats_d['f_dr'].append(initial_float)
    stats_d['f_pur'].append(initial_float)
    stats_d['f_com'].append(initial_float)
    stats_d['alpha_std'].append(initial_float)
    stats_d['delta_std'].append(initial_float)

    true_sources = [21, 22, 14, 17, 19, 22, 23, 19, 18, 21]

    for idx, pm_ in enumerate(prfs_d['pms']):
        stats_d['i_pm'].append(pm_)
        stats_d['N_meas'].append(initial_int)
        stats_d['N_false'].append(initial_int)
        stats_d['N_se'].append(initial_float)
        stats_d['N_true'].append(true_sources[idx])
        stats_d['f_dr'].append(initial_float)
        stats_d['f_pur'].append(initial_float)
        stats_d['f_com'].append(initial_float)
        stats_d['alpha_std'].append(initial_float)
        stats_d['delta_std'].append(initial_float)

    return stats_d


def redo_tmp_d():
    """ Creates a dictionary

    :return: tmp_d
    """
    tmp_d = {'sex_cf': [], 'scmp_cf': [], 'boolean_l': [], 'catalog': [],
             'CCD': [], 'source': [], 'epoch': [], 'i_pm': [],
             'i_pm_alpha': [], 'i_pm_delta': [], 'i_alpha': [],
             'i_delta': [], 'o_pm_alpha': [], 'o_pm_delta': [],
             'o_pm_alpha_err': [], 'o_pm_delta_err': [], 'o_alpha': [],
             'o_delta': [], 'error_a': [], 'error_b': []}

    return tmp_d


def redo_err_d():
    """ Creates a dictionary

    :return: err_d
    """
    err_d = {'sex_cf': [], 'scmp_cf': [], 'catalog': [],
             'CCD': [], 'source': [], 'scmp_source': [], 'i_pm': [],
             'i_pm_alpha': [], 'i_pm_delta': [], 'o_pm': [], 'o_pm_alpha': [],
             'o_pm_delta': [], 'i_alpha': [], 'i_delta': [], 'o_alpha': [],
             'o_delta': []}

    return err_d


def check_star(catalog_n, i_df, o_alpha, o_delta):
    """

    :param catalog_n:
    :param i_df:0
    :param o_alpha:
    :param o_delta:
    :return:
    """
    tolerance = 0.0001389  # 0.5 arcsecond

    i_df = i_df[i_df['catalog'].isin([catalog_n])]
    i_df = i_df[i_df['alpha_j2000'] + tolerance > o_alpha]
    i_df = i_df[o_alpha > i_df['alpha_j2000'] - tolerance]
    i_df = i_df[i_df['delta_j2000'] + tolerance > o_delta]
    i_df = i_df[o_delta > i_df['delta_j2000'] - tolerance]

    return i_df


def check_mag(i_df, o_alpha, o_delta):
    """

    :param i_df:
    :param o_alpha:
    :param o_delta:
    :return:
    """
    tolerance = 0.0001389  # 0.5 arcsecond

    i_df = i_df[i_df['alpha_j2000'] + tolerance > o_alpha]
    i_df = i_df[o_alpha > i_df['alpha_j2000'] - tolerance]
    i_df = i_df[i_df['delta_j2000'] + tolerance > o_delta]
    i_df = i_df[o_delta > i_df['delta_j2000'] - tolerance]

    if i_df.empty:
        flag_mag = False
    else:
        flag_mag = True

    return flag_mag


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


class ScampPerformanceStars:

    def __init__(self, logger, mag, sex_cf, scmp_cf):
        """

        """
        self.save = True
        self.norm_speed = False  # Detected SSOs are classified according
        self.filter_p_number = '7s'
        # their input pm
        self.logger = logger
        self.mag = mag
        self.sex_cf = sex_cf
        self.scmp_cf = scmp_cf

        self.bypassfilter = True
        self.prfs_d = extract_settings()

        self.check()

    def get_norm_speed(self, o_pm):
        """

        :return:
        """
        speeds_d = speeds_range(self.prfs_d, 50)

        pm_norm = 0
        for key_ in speeds_d.keys():
            low = speeds_d[key_][0]
            high = speeds_d[key_][1]
            if low < o_pm < high:
                pm_norm = key_

        return pm_norm

    def creates_input_dict(self):
        """ Creates an input dictionary. Each key contains SSOs' information
        for each dither.

        :return: input_dict
        """
        input_dict = {}
        # Loops over the four dithers
        for dither in range(1, 5, 1):
            cat_location = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'],
                                                   self.mag)
            cat_name = '{}/Cat_20-21_d{}'.format(cat_location, dither)
            input_dict[dither] = '{}.dat'.format(cat_name)
        input_dict = Create_regions(input_dict).check_ssos(self.mag, True)

        return input_dict

    def creates_input_df(self, input_dict):
        """ Creates an input dataframe from an input dictionary.

        :return: input dataframe
        """
        input_list = []
        for key_ in input_dict.keys():
            input_list.append(input_dict[key_])

        input_df = concat(input_list, axis=0)
        # Look for < 3 coincidences
        input_df = concat(g for _, g in input_df.groupby('source')
                          if len(g) >= 3)
        input_df = input_df.reset_index(drop=True)

        if self.save:
            input_df.to_csv('inputs.csv')

        return input_df

    def check(self):
        """

        :return:
        """
        self.logger.debug('Scamp configuration {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration {}'.format(self.sex_cf))

        # Creates a dictionary for statistics
        stats_d = redo_stats_d()
        stats_d = populate_stats_d(stats_d)

        # Creates a dictionary for mean and populated with lists
        nonsso_d = {'a_image': [], 'b_image': []}
        for idx in range(0, len(self.prfs_d['pms']) + 1, 1):
            nonsso_d['a_image'].append([])
            nonsso_d['b_image'].append([])

        # Creates a dictionary for standard deviation and populated with lists
        sso_d = {'a_image': [], 'b_image': []}
        for idx in range(0, len(self.prfs_d['pms']) + 1, 1):
            sso_d['a_image'].append([])
            sso_d['b_image'].append([])

        """
        # Creates a dictionary for standard deviation and populated with lists
        sso_d = {'theta_image': []}
        for idx in range(0, len(self.prfs_d['pms']) + 1, 1):
            sso_d['a_image'].append([])
            sso_d['b_image'].append([])
        """

        input_dict = self.creates_input_dict()
        input_df = self.creates_input_df(input_dict)

        # Gets the name of filtered file
        filter_n = 'filt_{}_{}_{}.csv'.format(self.scmp_cf, self.mag,
                                              self.filter_p_number)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filter_n)
        self.logger.debug('opens filtered catalog {}'.format(filter_n))
        # Opens filtered file
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)
        # Gets unique sources from filtered file
        o_unique_sources = list(set(o_cat['SOURCE_NUMBER'].tolist()))

        # A catalogue populated by stars
        cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], self.mag)
        cat_name = '{}/Cat_20-21_d1.dat'.format(cat_loc)
        input_stars = Create_regions(cat_name).get_stars(self.mag)
        sources_n = len(o_unique_sources)

        self.logger.debug('unique sources to be analysed {}'.format(sources_n))
        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(o_unique_sources):
            # self.logger.debug('idx position {}'.format(idx))
            # Initiate some values
            tmp_d = {'flag_sso': [], 'source': [], 'theta_image': [],
                     'a_image': [], 'b_image': [], 'i_alpha': [],
                     'i_delta': [], 'o_alpha': [], 'o_delta': [], 'i_pm': [],
                     'i_pm_alpha': [], 'i_pm_delta': [], 'o_pm': [],
                     'o_pm_alpha': [], 'o_pm_delta': [], 'o_pm_norm': []}

            # Check if actual source lies in 20-21 magnitude gap
            # flag_mag = False
            o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                source = row.SOURCE_NUMBER
                catalog_n = row.CATALOG_NUMBER
                # a_image = row.A_IMAGE
                # b_image = row.B_IMAGE
                theta_image = row.THETA_IMAGE
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000
                o_pm = row.PM
                o_pm_alpha = row.PMALPHA
                o_pm_delta = row.PMDELTA

                flag_mag = check_mag(input_stars, o_alpha, o_delta)
                out_df = check_star(catalog_n, input_df, o_alpha, o_delta)

                if out_df.empty and flag_mag:
                    tmp_d['flag_sso'].append(False)
                    tmp_d['source'].append(source)
                    tmp_d['theta_image'].append(theta_image)
                    tmp_d['o_alpha'].append(o_alpha)
                    tmp_d['o_delta'].append(o_delta)
                    tmp_d['o_pm'].append(o_pm)
                    tmp_d['o_pm_alpha'].append(o_pm_alpha)
                    tmp_d['o_pm_delta'].append(o_pm_delta)
                elif out_df.empty is not True:
                    # it's a SSO
                    tmp_d['flag_sso'].append(True)
                    tmp_d['source'].append(source)
                    tmp_d['theta_image'].append(theta_image)
                    i_alpha = out_df['alpha_j2000'].iloc[0]
                    tmp_d['i_alpha'].append(i_alpha)
                    i_delta = out_df['delta_j2000'].iloc[0]
                    tmp_d['i_delta'].append(i_delta)
                    tmp_d['o_alpha'].append(o_alpha)
                    tmp_d['o_delta'].append(o_delta)
                    i_pm = out_df['pm_values'].iloc[0]
                    if i_pm == 0.001:
                        print('source {}'.format(source))
                    tmp_d['i_pm'].append(i_pm)
                    i_pm_alpha = out_df['pm_alpha'].iloc[0]
                    tmp_d['i_pm_alpha'].append(i_pm_alpha)
                    i_pm_delta = out_df['pm_delta'].iloc[0]
                    tmp_d['i_pm_delta'].append(i_pm_delta)
                    tmp_d['o_pm'].append(o_pm)
                    tmp_d['o_pm_alpha'].append(o_pm_alpha)
                    tmp_d['o_pm_delta'].append(o_pm_delta)
                else:
                    pass

            if len(set(tmp_d['flag_sso'])) == 1 and tmp_d['flag_sso'][0] is True:
                flag_sso = True
            else:
                flag_sso = False
            if len(set(tmp_d['o_pm'])) == 1:
                flag_pm = True
            else:
                flag_pm = False

            # stats for non-SSOs
            if flag_sso is not True and flag_pm:
                o_pm_norm = self.get_norm_speed(tmp_d['o_pm'][0])
                for idx_pm in range(0, len(tmp_d['o_pm']), 1):
                    tmp_d['o_pm_norm'].append(o_pm_norm)

                idx = stats_d['i_pm'].index(o_pm_norm)

                """
                if len(tmp_d['source']) >= 3:
                    print('false {}'.format(len(tmp_d['source'])))
                else:
                    print('error')

                # nonsso_d['a_image'][idx].append(a_image)
                # nonsso_d['b_image'][idx].append(b_image)
                """
                stats_d['N_meas'][idx] += 1
                stats_d['N_false'][idx] += 1
            # stats for SSOsls
            elif flag_sso and flag_pm:
                """
                o_pm_norm = self.get_norm_speed(tmp_d['o_pm'][0])
                for idx_ssos in range(0, len(tmp_d['o_pm']), 1):
                    tmp_d['o_pm_norm'].append(o_pm_norm)

                idx = stats_d['i_pm'].index(o_pm_norm)
                """
                """
                if tmp_d['i_pm'][0] == 10.0:
                    print(tmp_d['source'])
                    print(tmp_d['o_pm'])
                    print(tmp_d['i_alpha'], tmp_d['i_delta'])
                    print(tmp_d['o_alpha'], tmp_d['o_delta'])
                    print(" ")
                """
                idx = stats_d['i_pm'].index(tmp_d['i_pm'][0])

                # sso_d['a_image'][idx].append(a_image)
                # sso_d['b_image'][idx].append(b_image)

                stats_d['N_meas'][idx] += 1
                stats_d['N_se'][idx] += 1

        # for key_ in ['a_image', 'b_image']:
        #     print(key_)
        #     for idx_pm, pm_ in enumerate(stats_d['i_pm']):
        #         print(idx_pm, pm_)
        #         print('non-sso mean {}'.format(mean(nonsso_d[key_][idx_pm])))
        #         print('non-sso std {}'.format(std(nonsso_d[key_][idx_pm])))
        #         print('sso mean {}'.format(mean(sso_d[key_][idx_pm])))
        #         print('sso std {}'.format(std(sso_d[key_][idx_pm])))

        stats_d = compute_factors(stats_d)
        stats_df = DataFrame(stats_d, columns=['i_pm', 'N_true', 'N_se',
                                               'N_false', 'N_meas', 'f_pur',
                                               'f_dr', 'f_com'])
        # Saves dataframe to csv file
        stats_df.to_csv('scamp_test_{}.csv'.format(self.filter_p_number))

        return stats_d
