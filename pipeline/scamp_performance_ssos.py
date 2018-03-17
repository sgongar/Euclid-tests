# !/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance
Mide las detecciones de cada tipo

Versions:
- 0.1

Todo:
    * Improve log messages

"""
from pandas import concat, read_csv, DataFrame
from numpy import mean, median, std

from misc import extract_settings, speeds_range
from regions import Create_regions


def redo_stats_d():
    """ Creates a dictionary

    :return: tmp_d
    """
    stats_d = {'i_pm': [], 'mean': [], 'std': [], 'min': [], 'max': []}

    return stats_d


def populate_stats_d(stats_d):
    """

    :param stats_d:
    :return:
    """
    prfs_d = extract_settings()
    initial_float = 0.0

    for idx, pm_ in enumerate(prfs_d['pms'][:-1]):
        stats_d['i_pm'].append(pm_)
        stats_d['mean'].append(initial_float)
        stats_d['std'].append(initial_float)
        stats_d['min'].append(initial_float)
        stats_d['max'].append(initial_float)

    return stats_d


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
        print(i_df['object'].iloc[0])

    return flag_mag


class ScampPerformanceSSOs:

    def __init__(self, logger, mag, sex_cf, scmp_cf):
        """
        Son estadisticas para las velocidades de entrada, no de salida.
        Eso quiere decir que no conozco la dispersion de los valores.
        Si no conozco la dispersion la mejora en los valores cuando ajuste b
        y class star se repartira entre todas las velocidades.

        """
        self.save = False
        self.norm_speed = False  # Detected SSOs are classified according
        self.filter_p_number = '9'
        # their input pm
        self.logger = logger
        self.mag = mag
        self.sex_cf = sex_cf
        self.scmp_cf = scmp_cf

        self.bypassfilter = True
        self.prfs_d = extract_settings()

        stars = True
        ssos = False
        if stars:
            self.check_stars()
        elif ssos:
            self.check_ssos()

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

    def check_ssos(self):
        """

        :return:
        """
        self.logger.debug('Scamp configuration {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration {}'.format(self.sex_cf))

        # Creates a dictionary for standard deviation and populated with lists
        sso_d = {'median_a_image': [], 'median_b_image': [],
                 'median_class_star': [], 'dispersion': [],
                 'ellipticity': []}
        for idx in range(0, len(self.prfs_d['pms']), 1):
            sso_d['median_a_image'].append([])
            sso_d['median_b_image'].append([])
            sso_d['median_class_star'].append([])
            sso_d['dispersion'].append([])
            sso_d['ellipticity'].append([])

        input_dict = self.creates_input_dict()
        input_df = self.creates_input_df(input_dict)

        # Gets the name of filtered file
        filter_n = 'filt_{}_{}_{}.csv'.format(self.scmp_cf, self.mag,
                                              self.filter_p_number)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filter_n)
        self.logger.debug('Opens filtered catalog {}'.format(filter_n))
        # Opens filtered file
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)
        # Gets unique sources from filtered file
        o_unique_sources = list(set(o_cat['SOURCE_NUMBER'].tolist()))

        # A catalogue populated by stars
        cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], self.mag)
        cat_name = '{}/Cat_20-21_d1.dat'.format(cat_loc)
        input_stars = Create_regions(cat_name).get_stars(self.mag)
        sources_n = len(o_unique_sources)

        self.logger.debug('Unique sources to be analysed {}'.format(sources_n))
        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(o_unique_sources):
            self.logger.debug('idx position {}'.format(idx))
            # Initiate some values
            tmp_d = {'flag_sso': [], 'source': [], 'theta_image': [],
                     'median_a_image': [], 'median_b_image': [], 'i_pm': [],
                     'o_pm': [], 'median_class_star': [], 'dispersion': [],
                     'ellipticity': []}

            # Check if actual source lies in 20-21 magnitude gap
            # flag_mag = False
            o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                source = row.SOURCE_NUMBER
                catalog_n = row.CATALOG_NUMBER
                median_a_image = row.MEDIAN_A_IMAGE  # Hardcoded to median
                median_b_image = row.MEDIAN_B_IMAGE
                theta_image = row.THETA_IMAGE
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000
                o_pm = row.PM
                median_class_star = row.MEDIAN_CLASS_STAR
                ellipticity = row.ELLIPTICITY

                flag_mag = check_mag(input_stars, o_alpha, o_delta)
                out_df = check_star(catalog_n, input_df, o_alpha, o_delta)

                if out_df.empty and flag_mag:
                    tmp_d['median_a_image'].append(median_a_image)
                    tmp_d['median_b_image'].append(median_b_image)
                    tmp_d['median_class_star'].append(median_class_star)
                    o_pm_norm = self.get_norm_speed(o_pm)
                    tmp_d['o_pm'].append(o_pm_norm)
                elif out_df.empty is not True:
                    # it's a SSO
                    tmp_d['flag_sso'].append(True)
                    tmp_d['source'].append(source)
                    tmp_d['median_a_image'].append(median_a_image)
                    tmp_d['median_b_image'].append(median_b_image)
                    tmp_d['theta_image'].append(theta_image)
                    tmp_d['median_class_star'].append(median_class_star)
                    i_pm = out_df['pm_values'].iloc[0]
                    tmp_d['i_pm'].append(i_pm)
                    tmp_d['o_pm'].append(o_pm)
                    tmp_d['dispersion'].append(i_pm / o_pm)
                    tmp_d['ellipticity'].append(ellipticity)
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
            # if flag_sso is False and flag_pm is False:
            if flag_sso is False:
                """
                try:
                    if tmp_d['o_pm'][0] == 0:
                        idx = 0
                        for median_a_image_ in tmp_d['median_a_image']:
                            sso_d['median_a_image'][idx].append(median_a_image_)
                        for median_b_image_ in tmp_d['median_b_image']:
                            sso_d['median_b_image'][idx].append(median_b_image_)
                        for median_class_star_ in tmp_d['median_class_star']:
                            sso_d['median_class_star'][idx].append(
                                median_class_star_)
                    else:
                        idx = self.prfs_d['pms'].index(tmp_d['o_pm'][0]) + 1
                        for median_a_image_ in tmp_d['median_a_image']:
                            sso_d['median_a_image'][idx].append(median_a_image_)
                        for median_b_image_ in tmp_d['median_b_image']:
                            sso_d['median_b_image'][idx].append(median_b_image_)
                        for median_class_star_ in tmp_d['median_class_star']:
                            sso_d['median_class_star'][idx].append(median_class_star_)
                except IndexError:
                    pass
                """
                pass
            # stats for SSOs
            elif flag_sso and flag_pm:
                idx = self.prfs_d['pms'].index(tmp_d['i_pm'][0])
                """
                for median_a_image_ in tmp_d['median_a_image']:
                    sso_d['median_a_image'][idx].append(median_a_image_)
                for median_b_image_ in tmp_d['median_b_image']:
                    sso_d['median_b_image'][idx].append(median_b_image_)
                for median_class_star_ in tmp_d['median_class_star']:
                    sso_d['median_class_star'][idx].append(median_class_star_)
                for dispersion_ in tmp_d['dispersion']:
                    sso_d['dispersion'][idx].append(dispersion_)
                for ellipticity_ in tmp_d['ellipticity']:
                    sso_d['ellipticity'][idx].append(ellipticity_)
                """
                sso_d['median_a_image'][idx].append(tmp_d['median_a_image'][0])
                sso_d['median_b_image'][idx].append(tmp_d['median_b_image'][0])
                sso_d['median_class_star'][idx].append(tmp_d['median_class_star'][0])
                sso_d['dispersion'][idx].append(tmp_d['dispersion'][0])
                sso_d['ellipticity'][idx].append(tmp_d['ellipticity'][0])

        for sso_key in sso_d.keys():
            """
            # Creates a dictionary for statistics
            stats_d = redo_stats_d()
            stats_d = populate_stats_d(stats_d)

            for stats_key in ['mean', 'std', 'min', 'max']:
                for idx_pm in range(0, len(stats_d['i_pm']), 1):
                    if stats_key == 'mean':
                        element_mean = mean(sso_d[sso_key][idx_pm])
                        stats_d[stats_key][idx_pm] = element_mean
                    elif stats_key == 'std':
                        element_std = std(sso_d[sso_key][idx_pm])
                        stats_d[stats_key][idx_pm] = element_std
                    elif stats_key == 'min':
                        try:
                            element_min = min(sso_d[sso_key][idx_pm])
                        except ValueError:
                            element_min = 'nan'
                        stats_d[stats_key][idx_pm] = element_min
                    elif stats_key == 'max':
                        try:
                            element_max = max(sso_d[sso_key][idx_pm])
                        except ValueError:
                            element_max = 'nan'
                        stats_d[stats_key][idx_pm] = element_max

            stats_df = DataFrame(stats_d, columns=['i_pm', 'mean', 'std',
                                                   'min', 'max'])
            # print(stats_df.describe())  # useful?
            stats_df.to_csv('stats/{}_{}.csv'.format(sso_key,
                                                     self.filter_p_number))
            """
            df2_keys = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]
            sso_d2 = {0.001: [], 0.003: [], 0.01: [], 0.03: [],
                      0.1: [], 0.3: [], 1: [], 3: [], 10: [], 30: []}

            sso_df = DataFrame(sso_d[sso_key])
            for column_ in sso_df.columns:
                for idx, value_ in enumerate(sso_df[column_]):
                    sso_d2[df2_keys[idx]].append(value_)

            sso_df.to_csv('data/{}_{}.csv'.format(sso_key,
                                                  self.filter_p_number))
            sso_df2 = DataFrame(sso_d2)
            sso_df2.to_csv('data/f_{}_{}.csv'.format(sso_key,
                                                     self.filter_p_number))

    def check_stars(self):
        """

        :return:
        """
        self.logger.debug('Scamp configuration {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration {}'.format(self.sex_cf))

        # Creates a dictionary
        sso_d = {'catalog_n': [], 'alpha_j2000': [], 'delta_j2000': [],
                 'median_a_image': [], 'median_b_image': [],
                 'median_class_star': [], 'dispersion': [], 'ellipticity': []}
        for idx in range(0, len(self.prfs_d['pms']) + 1, 1):
            sso_d['catalog_n'].append([])
            sso_d['alpha_j2000'].append([])
            sso_d['delta_j2000'].append([])
            sso_d['median_a_image'].append([])
            sso_d['median_b_image'].append([])
            sso_d['median_class_star'].append([])
            sso_d['dispersion'].append([])
            sso_d['ellipticity'].append([])

        input_dict = self.creates_input_dict()
        input_df = self.creates_input_df(input_dict)

        # Gets the name of filtered file
        filter_n = 'filt_{}_{}_{}.csv'.format(self.scmp_cf, self.mag,
                                              self.filter_p_number)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filter_n)
        self.logger.debug('Opens filtered catalog {}'.format(filter_n))
        # Opens filtered file
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)
        # Gets unique sources from filtered file
        o_unique_sources = list(set(o_cat['SOURCE_NUMBER'].tolist()))

        # A catalogue populated by stars
        cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'], self.mag)
        cat_name = '{}/Cat_20-21_d1.dat'.format(cat_loc)
        input_stars = Create_regions(cat_name).get_stars(self.mag)
        sources_n = len(o_unique_sources)

        source_l = []
        self.logger.debug('Unique sources to be analysed {}'.format(sources_n))
        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(o_unique_sources):
            self.logger.debug('idx position {}'.format(idx))
            # Initiate some values
            tmp_d = {'flag_sso': [], 'flag_mag': [], 'source': [],
                     'catalog_n': [], 'alpha_j2000': [], 'delta_j2000': [],
                     'theta_image': [], 'median_a_image': [],
                     'median_b_image': [], 'i_pm': [], 'o_pm': [],
                     'median_class_star': [], 'dispersion': [],
                     'ellipticity': []}

            # Check if actual source lies in 20-21 magnitude gap
            # flag_mag = False
            o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                source = row.SOURCE_NUMBER
                catalog_n = row.CATALOG_NUMBER
                median_a_image = row.MEDIAN_A_IMAGE  # Hardcoded to median
                median_b_image = row.MEDIAN_B_IMAGE
                # theta_image = row.THETA_IMAGE
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000
                o_pm = row.PM
                median_class_star = row.MEDIAN_CLASS_STAR
                # ellipticity = row.ELLIPTICITY

                flag_mag = check_mag(input_stars, o_alpha, o_delta)
                tmp_d['flag_mag'].append(flag_mag)
                out_df = check_star(catalog_n, input_df, o_alpha, o_delta)

                if out_df.empty and flag_mag:
                    tmp_d['catalog_n'].append(catalog_n)
                    tmp_d['alpha_j2000'].append(o_alpha)
                    tmp_d['delta_j2000'].append(o_delta)
                    tmp_d['median_a_image'].append(median_a_image)
                    tmp_d['median_b_image'].append(median_b_image)
                    tmp_d['median_class_star'].append(median_class_star)
                    o_pm_norm = self.get_norm_speed(o_pm)
                    tmp_d['o_pm'].append(o_pm_norm)
                    tmp_d['source'].append(source)
                # SSO todo improve definition
                elif out_df.empty is not True:
                    pass
                # todo improve definition
                else:
                    pass

            if len(set(tmp_d['flag_sso'])) == 1 and tmp_d['flag_sso'][0] is True:
                flag_sso = True
            else:
                flag_sso = False

            if len(set(tmp_d['flag_mag'])) == 1 and tmp_d['flag_mag'][0] is True:
                flag_mag = True
            else:
                flag_mag = False

            # stats for non-SSOs
            # if flag_sso is False and flag_pm is False:
            if flag_sso is False and flag_mag:
                if tmp_d['o_pm'][0] == 0:
                    idx = 0
                else:
                    idx = self.prfs_d['pms'].index(tmp_d['o_pm'][0]) + 1
                median_a_image_ = tmp_d['median_a_image'][0]
                sso_d['median_a_image'][idx].append(median_a_image_)
                median_b_image_ = tmp_d['median_b_image'][0]
                sso_d['median_b_image'][idx].append(median_b_image_)
                median_class_star_ = tmp_d['median_class_star'][0]
                sso_d['median_class_star'][idx].append(median_class_star_)
                alpha_j2000_ = tmp_d['alpha_j2000'][0]
                sso_d['alpha_j2000'][idx].append(alpha_j2000_)
                delta_j2000_ = tmp_d['delta_j2000'][0]
                sso_d['delta_j2000'][idx].append(delta_j2000_)
                for catalog_n_ in tmp_d['catalog_n']:
                    sso_d['catalog_n'][idx].append(catalog_n_)
                source_ = tmp_d['source'][0]
                source_l.append(source_)

        for sso_key in sso_d.keys():
            df2_keys = [0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]
            sso_d2 = {0: [], 0.001: [], 0.003: [], 0.01: [], 0.03: [],
                      0.1: [], 0.3: [], 1: [], 3: [], 10: [], 30: []}

            sso_df = DataFrame(sso_d[sso_key])
            for column_ in sso_df.columns:
                for idx, value_ in enumerate(sso_df[column_]):
                    sso_d2[df2_keys[idx]].append(value_)

            sso_df2 = DataFrame(sso_d2)
            sso_df2.to_csv('data_stars/f_{}_{}.csv'.format(sso_key,
                                                           self.filter_p_number))


        print(source_l)
