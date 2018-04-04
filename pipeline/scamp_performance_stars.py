#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance
Obtiene los factores

Versions:
- 0.1

Todo:
    * Improve log messages

"""
from numpy import array, isnan, nanmean, nan
from pandas import concat, read_csv, DataFrame

from misc import extract_settings, speeds_range, get_dither
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


def populate_stats_d(stats_d, mag):
    """

    :param stats_d:
    :param mag:
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

    # re-check!
    true_sources_d = {'20-21': [21, 22, 14, 17, 19, 22, 23, 19, 18, 21],
                      '21-22': [22, 21, 19, 24, 16, 23, 18, 24, 24, 19],
                      '22-23': [18, 24, 21, 22, 16, 22, 20, 20, 21, 26],
                      '23-24': [21, 19, 25, 23, 22, 19, 20, 17, 20, 21],
                      '24-25': [22, 19, 20, 14, 22, 17, 21, 20, 21, 18],
                      # why so many?
                      '25-26': [14, 23, 19, 27, 25, 29, 23, 21, 16, 22],
                      '26-27': [21, 20, 22, 22, 22, 20, 18, 25, 23, 14]}

    for idx, pm_ in enumerate(prfs_d['pms']):
        stats_d['i_pm'].append(pm_)
        stats_d['N_meas'].append(initial_int)
        stats_d['N_false'].append(initial_int)
        stats_d['N_se'].append(initial_float)
        stats_d['N_true'].append(true_sources_d[mag][idx])
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
    tmp_d = {'flag_sso': [], 'source': [], 'catalog_n': [], 'theta_image': [],
             'fwhm_image': [], 'a_image': [], 'a_image_median': [],
             'b_image': [], 'b_image_median': [], 'b_image_med/a_image_med': [],
             'erra_image': [], 'errb_image': [], 'i_alpha': [], 'i_delta': [],
             'o_alpha': [], 'o_delta': [], 'i_pm': [], 'i_pm_alpha': [],
             'i_pm_delta': [], 'o_pm': [], 'o_pm_alpha': [], 'o_pm_delta': [],
             'o_pm_norm': [], 'object': []}

    return tmp_d


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
    tolerance = 0.0001389  # 0.5 arcsecondcp

    i_df = i_df[i_df['alpha_j2000'] + tolerance > o_alpha]
    i_df = i_df[o_alpha > i_df['alpha_j2000'] - tolerance]
    i_df = i_df[i_df['delta_j2000'] + tolerance > o_delta]
    i_df = i_df[o_delta > i_df['delta_j2000'] - tolerance]

    if i_df.empty:
        flag_mag = False
        object_ = ''
    else:
        flag_mag = True
        object_ = i_df['object'].iloc[0]

    return flag_mag, object_


class ScampPerformanceStars:

    def __init__(self, logger, mag, sex_cf, scmp_cf):
        """

        """
        self.save = True
        self.norm_speed = False  # Detected SSOs are classified according
        self.filter_p_number = '9'
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
            cat_name = '{}/Cat_{}_d{}'.format(cat_location, self.mag, dither)
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
            input_df.to_csv('inputs_{}.csv'.format(self.mag))

        return input_df

    def check(self):
        """

        :return:
        """
        self.logger.debug('Magnitude bin: {}'.format(self.mag))
        self.logger.debug('Scamp configuration: {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration: {}'.format(self.sex_cf))

        # Creates a dictionary for statistics
        stats_d = redo_stats_d()
        stats_d = populate_stats_d(stats_d, self.mag)

        # Creates a (dictionary?) list for objects
        sso_d = {'idx': [], 'source': [], 'CCD': [], 'dither': [],
                 'a_image_median': [], 'b_image_median': [], 'theta_image': [],
                 'fwhm_image': [], 'b_image_med/a_image_med': [],
                 'erra_image': [], 'errb_image': [], 'alpha': [], 'delta': [],
                 'i_pm': [], 'o_pm': []}
        objects_d = {'idx': [], 'source': [], 'CCD': [], 'dither': [],
                     'a_image_median': [], 'b_image_median': [],
                     'theta_image': [], 'fwhm_image': [],
                     'b_image_med/a_image_med': [], 'erra_image': [],
                     'errb_image': [], 'alpha': [], 'delta': [], 'n_pm': [],
                     'o_pm': []}

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
        cat_name = '{}/Cat_{}_d1.dat'.format(cat_loc, self.mag)
        input_stars = Create_regions(cat_name).get_stars(self.mag)

        sources_n = len(o_unique_sources)

        self.logger.debug('unique sources to be analysed {}'.format(sources_n))
        # Loops over unique sources of filtered file
        idx_false = 0
        idx_true = 0
        for idx, source_ in enumerate(o_unique_sources):
            self.logger.debug('idx position {}'.format(idx))
            # Initiate some values
            tmp_d = redo_tmp_d()

            # Check if actual source lies in 20-21 magnitude gap
            # flag_mag = False
            o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                source = row.SOURCE_NUMBER
                catalog_n = row.CATALOG_NUMBER
                a_image_median = row.MEDIAN_A_IMAGE
                b_image_median = row.MEDIAN_B_IMAGE
                erra_image = row.ERRA_IMAGE
                errb_image = row.ERRB_IMAGE
                theta_image = row.THETA_IMAGE
                fwhm_image = row.FWHM_IMAGE
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000
                o_pm = row.PM
                o_pm_alpha = row.PMALPHA
                o_pm_delta = row.PMDELTA
                # mag_iso_median = row.MEDIAN_MAG_ISO

                flag_mag, object_ = check_mag(input_stars, o_alpha, o_delta)
                out_df = check_star(catalog_n, input_df, o_alpha, o_delta)

                if out_df.empty and flag_mag:
                    tmp_d['flag_sso'].append(False)
                    tmp_d['catalog_n'].append(catalog_n)
                    tmp_d['source'].append(source_)
                    tmp_d['a_image_median'].append(a_image_median)
                    tmp_d['b_image_median'].append(b_image_median)
                    b_vs_a_median = b_image_median/a_image_median
                    tmp_d['b_image_med/a_image_med'].append(b_vs_a_median)
                    tmp_d['erra_image'].append(erra_image)
                    tmp_d['errb_image'].append(errb_image)
                    tmp_d['theta_image'].append(theta_image)
                    tmp_d['fwhm_image'].append(fwhm_image)
                    tmp_d['o_alpha'].append(o_alpha)
                    tmp_d['o_delta'].append(o_delta)
                    tmp_d['o_pm'].append(o_pm)
                    tmp_d['o_pm_alpha'].append(o_pm_alpha)
                    tmp_d['o_pm_delta'].append(o_pm_delta)
                    tmp_d['object'].append(object_)
                elif out_df.empty is not True:
                    # it's a SSO
                    tmp_d['flag_sso'].append(True)
                    tmp_d['catalog_n'].append(catalog_n)
                    tmp_d['source'].append(source)
                    tmp_d['a_image_median'].append(a_image_median)
                    tmp_d['b_image_median'].append(b_image_median)
                    b_vs_a_median = b_image_median/a_image_median
                    tmp_d['b_image_med/a_image_med'].append(b_vs_a_median)
                    tmp_d['erra_image'].append(erra_image)
                    tmp_d['errb_image'].append(errb_image)
                    tmp_d['theta_image'].append(theta_image)
                    tmp_d['fwhm_image'].append(fwhm_image)
                    i_alpha = out_df['alpha_j2000'].iloc[0]
                    tmp_d['i_alpha'].append(i_alpha)
                    i_delta = out_df['delta_j2000'].iloc[0]
                    tmp_d['i_delta'].append(i_delta)
                    tmp_d['o_alpha'].append(o_alpha)
                    tmp_d['o_delta'].append(o_delta)
                    i_pm = out_df['pm_values'].iloc[0]
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
                # Normalise speed
                o_pm_norm = self.get_norm_speed(tmp_d['o_pm'][0])
                for idx_pm in range(0, len(tmp_d['o_pm']), 1):
                    # tmp_d['o_pm_norm'].append(o_pm_norm)
                    objects_d['n_pm'].append(o_pm_norm)

                idx = stats_d['i_pm'].index(o_pm_norm)
                stats_d['N_meas'][idx] += 1
                stats_d['N_false'][idx] += 1

                for catalog_n_ in tmp_d['catalog_n']:
                    objects_d['idx'].append(idx_false)
                    ccd, dither = get_dither(catalog_n_)
                    objects_d['CCD'].append(ccd)
                    objects_d['dither'].append(dither)
                for source_n in tmp_d['source']:
                    objects_d['source'].append(int(source_n))
                for a_image_median_ in tmp_d['a_image_median']:
                    objects_d['a_image_median'].append(a_image_median_)
                for b_image_median_ in tmp_d['b_image_median']:
                    objects_d['b_image_median'].append(b_image_median_)
                for theta_image_ in tmp_d['theta_image']:
                    objects_d['theta_image'].append(theta_image_)
                for b_vs_a_median_ in tmp_d['b_image_med/a_image_med']:
                    objects_d['b_image_med/a_image_med'].append(b_vs_a_median_)
                for erra_image_ in tmp_d['erra_image']:
                    objects_d['erra_image'].append(erra_image_)
                for errb_image_ in tmp_d['errb_image']:
                    objects_d['errb_image'].append(errb_image_)
                for fwhm_image_ in tmp_d['fwhm_image']:
                    objects_d['fwhm_image'].append(fwhm_image_)
                for alpha_ in tmp_d['o_alpha']:
                    objects_d['alpha'].append(alpha_)
                for delta_ in tmp_d['o_delta']:
                    objects_d['delta'].append(delta_)
                for o_pm_ in tmp_d['o_pm']:
                    objects_d['o_pm'].append(o_pm_)

                idx_false += 1

                # objects_d.append(tmp_d['object'][0])
            # stats for SSO
            elif flag_sso and flag_pm:
                idx = stats_d['i_pm'].index(tmp_d['i_pm'][0])

                stats_d['N_meas'][idx] += 1
                stats_d['N_se'][idx] += 1

                for catalog_n_ in tmp_d['catalog_n']:
                    sso_d['idx'].append(idx_true)
                    ccd, dither = get_dither(catalog_n_)
                    sso_d['CCD'].append(ccd)
                    sso_d['dither'].append(dither)
                for source_n in tmp_d['source']:
                    sso_d['source'].append(int(source_n))
                for a_image_median_ in tmp_d['a_image_median']:
                    sso_d['a_image_median'].append(a_image_median_)
                for b_image_median_ in tmp_d['b_image_median']:
                    sso_d['b_image_median'].append(b_image_median_)
                for theta_image_ in tmp_d['theta_image']:
                    sso_d['theta_image'].append(theta_image_)
                for b_vs_a_median_ in tmp_d['b_image_med/a_image_med']:
                    sso_d['b_image_med/a_image_med'].append(b_vs_a_median_)
                for erra_image_ in tmp_d['erra_image']:
                    sso_d['erra_image'].append(erra_image_)
                for errb_image_ in tmp_d['errb_image']:
                    sso_d['errb_image'].append(errb_image_)
                for fwhm_image_ in tmp_d['fwhm_image']:
                    sso_d['fwhm_image'].append(fwhm_image_)
                for alpha_ in tmp_d['o_alpha']:
                    sso_d['alpha'].append(alpha_)
                for delta_ in tmp_d['o_delta']:
                    sso_d['delta'].append(delta_)
                for i_pm_ in tmp_d['i_pm']:
                    sso_d['i_pm'].append(i_pm_)
                for o_pm_ in tmp_d['o_pm']:
                    sso_d['o_pm'].append(o_pm_)

                idx_true += 1

        stats_d = compute_factors(stats_d)
        stats_df = DataFrame(stats_d, columns=['i_pm', 'N_true', 'N_se',
                                               'N_false', 'N_meas', 'f_pur',
                                               'f_dr', 'f_com'])
        # Saves dataframe to csv file
        stats_df.to_csv('scamp_test_{}_{}.csv'.format(self.mag,
                                                      self.filter_p_number))

        objects_df = DataFrame(objects_d, columns=['idx', 'source', 'CCD',
                                                   'dither', 'a_image_median',
                                                   'b_image_median',
                                                   'theta_image',
                                                   'b_image_med/a_image_med',
                                                   'erra_image', 'errb_image',
                                                   'fwhm_image', 'alpha',
                                                   'delta', 'n_pm', 'o_pm'])
        objects_df.to_csv('false_positives_{}.csv'.format(self.mag))

        sso_df = DataFrame(sso_d, columns=['idx', 'source', 'CCD', 'dither',
                                           'a_image_median', 'b_image_median',
                                           'theta_image',
                                           'b_image_med/a_image_med',
                                           'erra_image', 'errb_image',
                                           'fwhm_image', 'alpha', 'delta',
                                           'i_pm', 'o_pm'])
        sso_df.to_csv('true_positives_{}.csv'.format(self.mag))

        return stats_d


class TotalScampPerformanceStars:
    def __init__(self, logger, mag, sex_cf, scmp_cf):
        """ Me da los valores de salida de todos las estrellas  y galaxias
        presentes en filt 3 obtenidos o no

        """
        self.logger = logger
        self.filter_p_number = 3
        self.prfs_d = extract_settings()

        self.mag = mag
        self.scmp_cf = scmp_cf
        self.sex_cf = sex_cf

        self.save = True

        self.data_d = self.creates_output_dict()
        input_sources = self.check_pm_distribution()
        self.get_data(input_sources)

    def check_source(self, catalog_n, o_cat, i_alpha, i_delta):
        """

        :param catalog_n:
        :param o_cat:
        :param i_alpha:
        :param i_delta:
        :return:
        """
        o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
        o_df = o_df[o_df['ALPHA_J2000'] + self.prfs_d['tolerance'] > i_alpha]
        o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - self.prfs_d['tolerance']]
        o_df = o_df[o_df['DELTA_J2000'] + self.prfs_d['tolerance'] > i_delta]
        o_df = o_df[i_delta > o_df['DELTA_J2000'] - self.prfs_d['tolerance']]

        return o_df

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
            cat_name = '{}/Cat_{}_d{}'.format(cat_location, self.mag, dither)
            input_dict[dither] = '{}.dat'.format(cat_name)
        input_dict = Create_regions(input_dict).check_ssos(self.mag, True)

        return input_dict

    def creates_input_df(self, input_dict):
        """ Creates an input dataframe from an input dictionary.

        :return: input dataframe
        """
        occurrences = 3

        input_list = []
        for key_ in input_dict.keys():
            input_list.append(input_dict[key_])

        input_df = concat(input_list, axis=0)
        # Look for >= 3 coincidences
        input_df = concat(g for _, g in input_df.groupby('source')
                          if len(g) >= occurrences)
        input_df = input_df.reset_index(drop=True)

        if self.save:
            self.logger.debug('Saves input catalog')
            input_df.to_csv('full_stats_ssos/inputs_{}.csv'.format(self.mag))

        return input_df

    def creates_output_dict(self):
        """

        :return: data_d
        """
        data_d = {}
        for pm_ in self.prfs_d['pms']:
            data_d[pm_] = []

        return data_d

    def gets_filtered_catalog(self):
        """

        :return: filtered_cat, filtered catalog
        """
        filter_n = 'filt_{}_{}_{}.csv'.format(self.scmp_cf, self.mag,
                                              self.filter_p_number)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filter_n)

        self.logger.debug('Opens filtered catalog {}'.format(filter_n))
        # Cross with filtered data - Opens datafile
        filtered_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        return filtered_cat

    def check_pm_distribution(self):
        """

        :return:
        """
        # Creates an input dictionary with all input sources
        self.logger.debug('Magnitude bin: {}'.format(self.mag))
        self.logger.debug('Scamp configuration: {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration: {}'.format(self.sex_cf))

        # Creates a dictionary
        sso_d = {'catalog_n': [], 'median_a_image': [],
                 'median_erra_image': [], 'median_b_image': [],
                 'median_errb_image': [], 'median_class_star': [],
                 'ellipticity': [], 'median_mag_iso': [],
                 'median_magerr_iso': [], 'output_pm': [],
                 'median_flux_iso': []}
        for idx in range(0, len(self.prfs_d['pms']), 1):
            sso_d['median_mag_iso'].append([])
            sso_d['median_magerr_iso'].append([])
            sso_d['catalog_n'].append([])
            sso_d['median_a_image'].append([])
            sso_d['median_erra_image'].append([])
            sso_d['median_b_image'].append([])
            sso_d['median_errb_image'].append([])
            sso_d['median_class_star'].append([])
            sso_d['ellipticity'].append([])
            sso_d['output_pm'].append([])
            sso_d['median_flux_iso'].append([])

        input_dict = self.creates_input_dict()
        input_df = self.creates_input_df(input_dict)

        # Open particular file!
        filter_cat = self.gets_filtered_catalog()

        # Gets unique sources from input data
        unique_sources = list(set(input_df['source'].tolist()))
        input_sources = []
        sources_n = len(unique_sources)
        self.logger.debug('Input sources to be analysed {}'.format(sources_n))
        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            # print('source {}'.format(source_))
            tmp_d = redo_tmp_d()
            flag_sso = False
            # i_pm = 0
            # Gets associated data in input catalog
            cat_df = input_df[input_df['source'].isin([source_])]
            # Iterate over each detection of each source
            for i, row in enumerate(cat_df.itertuples(), 1):
                catalog_n = row.catalog
                i_pm = row.pm_values
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # Checks if there is a source closed to input one
                o_df = self.check_source(catalog_n, filter_cat,
                                         i_alpha, i_delta)

                if o_df.empty is not True and o_df['PM'].size == 1:
                    scmp_source = o_df['SOURCE_NUMBER'].iloc[0]
                    input_sources.append(int(scmp_source))
                else:
                    pass

        return list(set(input_sources))

    def get_data(self, input_sources):
        """

        :param input_sources:
        :return:
        """

        star_d = {'catalog_n': [], 'median_a_image': [],
                  'median_erra_image': [], 'median_b_image': [],
                  'median_errb_image': [], 'median_class_star': [],
                  'ellipticity': [], 'median_mag_iso': [],
                  'median_magerr_iso': [], 'output_pm': [],
                  'median_flux_iso': []}
        for idx in range(0, len(self.prfs_d['pms']) + 1, 1):
            star_d['median_mag_iso'].append([])
            star_d['median_magerr_iso'].append([])
            star_d['catalog_n'].append([])
            star_d['median_a_image'].append([])
            star_d['median_erra_image'].append([])
            star_d['median_b_image'].append([])
            star_d['median_errb_image'].append([])
            star_d['median_class_star'].append([])
            star_d['ellipticity'].append([])
            star_d['output_pm'].append([])
            star_d['median_flux_iso'].append([])
        o_df = self.gets_filtered_catalog()
        stars_df = o_df[~o_df['SOURCE_NUMBER'].isin(input_sources)]

        # Gets unique sources from input data
        unique_sources = list(set(stars_df['SOURCE_NUMBER'].tolist()))
        unique_sources = array(unique_sources)
        unique_sources = unique_sources[~isnan(unique_sources)]

        unique_sources = [value for value in unique_sources if value != nan]
        sources_n = len(unique_sources)

        self.logger.debug('Fixed sources to be analysed {}'.format(sources_n))
        # Loops over input data (Luca's catalog)
        for idx_source, source_ in enumerate(unique_sources):
            # print('source {}'.format(source_))
            tmp_d = redo_tmp_d()
            flag_sso = False
            # i_pm = 0
            print(idx_source)
            # Gets associated data in input catalog
            cat_df = stars_df[stars_df['SOURCE_NUMBER'].isin([source_])]

            o_pm = cat_df['PM'].iloc[0]
            o_pm_norm = self.get_norm_speed(o_pm)

            pms = [0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
            idx = pms.index(o_pm_norm)

            median_mag_iso = cat_df['MEDIAN_MAG_ISO'].iloc[0]
            star_d['median_mag_iso'][idx].append(median_mag_iso)
            median_magerr_iso = cat_df['MEDIAN_MAGERR_ISO'].iloc[0]
            star_d['median_magerr_iso'][idx].append(median_magerr_iso)
            median_a_image_ = cat_df['MEDIAN_A_IMAGE'].iloc[0]
            star_d['median_a_image'][idx].append(median_a_image_)
            median_erra_image_ = cat_df['MEDIAN_ERRA_IMAGE'].iloc[0]
            star_d['median_erra_image'][idx].append(median_erra_image_)
            median_b_image_ = cat_df['MEDIAN_B_IMAGE'].iloc[0]
            star_d['median_b_image'][idx].append(median_b_image_)
            median_errb_image_ = cat_df['MEDIAN_ERRB_IMAGE'].iloc[0]
            star_d['median_errb_image'][idx].append(median_errb_image_)
            median_class_star_ = cat_df['MEDIAN_CLASS_STAR'].iloc[0]
            star_d['median_class_star'][idx].append(median_class_star_)
            ellipticity_ = cat_df['ELLIPTICITY'].iloc[0]
            star_d['ellipticity'][idx].append(ellipticity_)
            output_pm_ = o_pm  # really needed?
            star_d['output_pm'][idx].append(output_pm_)
            median_flux_iso_ = cat_df['MEDIAN_FLUX_ISO'].iloc[0]
            star_d['median_flux_iso'][idx].append(median_flux_iso_)

        keys = ['median_mag_iso', 'median_magerr_iso', 'median_a_image',
                'median_erra_image', 'median_b_image', 'median_errb_image',
                'median_flux_iso', 'median_class_star', 'ellipticity',
                'output_pm']

        for star_key in keys:
            df2_keys = [0, 0.001, 0.003, 0.01, 0.03,
                        0.1, 0.3, 1, 3, 10, 30]
            sso_d2 = {0: [], 0.001: [], 0.003: [], 0.01: [], 0.03: [],
                      0.1: [], 0.3: [], 1: [], 3: [], 10: [], 30: []}

            sso_df = DataFrame(star_d[star_key])
            for column_ in sso_df.columns:
                for idx, value_ in enumerate(sso_df[column_]):
                    sso_d2[df2_keys[idx]].append(value_)

            for key_ in sso_d2.keys():
                sso_d2[key_].append('mean {}'.format(nanmean(sso_d2[key_])))

            sso_df2 = DataFrame(sso_d2)
            sso_df2_filename = 'f_{}_{}_{}.csv'.format(self.mag, star_key,
                                                       self.filter_p_number)
            sso_df2.to_csv('full_stats_stars/{}'.format(sso_df2_filename))