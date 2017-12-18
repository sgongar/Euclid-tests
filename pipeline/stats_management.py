#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Todo:
    * Improve log messages

"""

from collections import Counter
from os import listdir, getcwd, path, makedirs

from astropy.io import fits
from astropy.table import Table
from numpy import sqrt
from pandas import concat, DataFrame, read_csv

from misc import extract_settings, speeds_range, all_same

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def search(o_cat, catalog_n, i_alpha, i_delta):
    """

    :param o_cat:
    :param i_alpha:
    :param i_delta:
    :return:
    """
    tolerance = 0.0001

    o_df = o_cat[o_cat['catalog'].isin([catalog_n])]
    o_df = o_df[o_df['alpha_j2000'] + tolerance > i_alpha]
    o_df = o_df[i_alpha > o_df['alpha_j2000'] - tolerance]
    o_df = o_df[o_df['delta_j2000'] + tolerance > i_delta]
    o_df = o_df[i_delta > o_df['delta_j2000'] - tolerance]

    return o_df


def gets_speed(pm, speeds):
    """

    :param pm:
    :param speeds:
    :return:
    """
    output_pm = 0.0

    for key_ in speeds.keys():
        if speeds[key_][1] > pm > speeds[key_][0]:
            output_pm = float(key_)

    return output_pm


class ExtractStats:

    def __init__(self, logger, mag, sex_cf, sex_d, scmp_cf, scmp_d):
        """

        :param logger:
        :param mag:
        :param sex_cf:
        :param scmp_cf:
        """
        self.prfs_d = extract_settings()

        self.mag = mag
        self.sex_cf = sex_cf
        self.sex_d = sex_d
        self.scmp_cf = scmp_cf
        self.scmp_d = scmp_d

        self.logger = logger

        (df_filter, df_input, valid_sources,
         data_full, data_merged) = self.extract_stats()

        self.filter_stats(df_filter, df_input, valid_sources,
                          data_full, data_merged)

    def extract_stats(self):
        """

        @param logger: a logger object.
        @param mag: magnitude gap to be analysed.

        @return stats_d:
        """
        # Gets all sources detected by scamp after filtering process
        # Opens csv file
        filt_cat_dir = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'], self.mag,
                                            self.sex_cf, self.scmp_cf)
        filt_cat_name = 'filt_{}_{}_2.csv'.format(self.scmp_cf, self.mag)
        filt_cat_loc = '{}/{}'.format(filt_cat_dir, filt_cat_name)
        self.logger.debug('reading filtered catalog {}'.format(filt_cat_loc))

        df_filter = read_csv(filt_cat_loc, index_col=0)
        df_filter = concat(g for _, g in df_filter.groupby("SOURCE_NUMBER")
                           if len(g) >= int(self.prfs_d['detections']))

        # Gets all unique sources after filter
        self.logger.debug('getting unique sources from post-filter scamp catalog')
        valid_sources = list(set(df_filter['SOURCE_NUMBER'].tolist()))

        # Gets all sources created from Luca's
        df_input = read_csv('input_sources.csv', index_col=0)

        results_dir = '{}/{}/{}/{}'.format(self.prfs_d['catalogs_dir'],
                                           self.mag, self.sex_cf, self.scmp_cf)
        full_n = '{}/full_{}_{}_1.cat'.format(results_dir, self.scmp_cf,
                                              self.mag)
        full_cat = fits.open(full_n)
        data_full = Table(full_cat[2].data).to_pandas()
        data_full = concat(g for _, g in data_full.groupby('SOURCE_NUMBER')
                           if len(g) >= 3)
        data_full = data_full.reset_index(drop=True)

        merged_n = '{}/merged_{}_{}_1.cat'.format(results_dir, self.scmp_cf,
                                                  self.mag)
        merged_cat = fits.open(merged_n)
        data_merged = Table(merged_cat[2].data).to_pandas()

        return df_filter, df_input, valid_sources, data_full, data_merged

    def filter_stats(self, df_filter, df_input, valid_sources,
                     data_full, data_merged):
        """

        :param df_filter:
        :param df_input:
        :param valid_sources:
        :return:
        """
        # Confidence values to be tested
        confidence_list = self.prfs_d['confidences']
        # Created dictionaries for all outputs
        filt_d, pm_d, stats_d = create_dicts()  # fixme check dictionaries vars

        confidence_list = [100]
        # Looks for sources one by one in sextractor/scamp catalog
        # Try to find if they're an SSO or not
        # source_ number refered to custom catalog, not sextractor/scamp one
        # Confidence values, how pm can be away from real value
        for confidence_ in confidence_list:
            # total_pms is a list with all proper motions
            # total values is a dict with all proper motions
            # and their frequencies
            total_pms = self.speeds_distribution(df_input)
            total_values = dict(Counter(total_pms))

            # values is a list with all different proper motions availables
            values = sorted(total_values.keys())

            # Harcoded TODO
            speeds = speeds_range(self.prfs_d, confidence_)

            empty_value = 0.0
            for value_ in values:
                stats_d['confidence'].append(confidence_)
                position_maxerr = self.scmp_d['position_maxerr']
                stats_d['position_maxerr'].append(position_maxerr)
                posangle_maxerr = self.scmp_d['posangle_maxerr']
                stats_d['posangle_maxerr'].append(posangle_maxerr)
                pixscale_maxerr = self.scmp_d['pixscale_maxerr']
                stats_d['pixscale_maxerr'].append(pixscale_maxerr)
                crossid_radius = self.scmp_d['crossid_radius']
                deblending = self.sex_d['deblending']
                stats_d['deblending'].append(deblending)
                mincount = self.sex_d['mincount']
                stats_d['mincount'].append(mincount)
                threshold = self.sex_d['threshold']
                stats_d['threshold'].append(threshold)
                area = self.sex_d['area']
                stats_d['area'].append(area)
                stats_d['cross_id'].append(crossid_radius)
                stats_d['mag'].append(self.mag)
                stats_d['pm'].append(value_)
                stats_d['N_true'].append(total_values[value_])
                stats_d['N_meas'].append(empty_value)
                stats_d['N_se'].append(empty_value)
                stats_d['N_meas-N_se'].append(empty_value)
                stats_d['f_dr'].append(empty_value)
                stats_d['f_pur'].append(empty_value)
                stats_d['f_com'].append(empty_value)

            sources = dict(Counter(df_filter['SOURCE_NUMBER'].tolist()))
            sources = sources.keys()

            # Mira a todas las fuentes
            for source_ in sources:
                tmp_flag = []
                tmp_pm_in = []
                # Saca los datos del catalogo para cada objeto
                o_source = data_merged[data_merged['SOURCE_NUMBER'].isin([source_])]
                o_pm_alpha = float(o_source['PMALPHA_J2000']) / 8.75e6
                o_pm_delta = float(o_source['PMDELTA_J2000']) / 8.75e6

                # Saca el proper motion de cada objeto
                o_pm = sqrt(o_pm_alpha ** 2 + o_pm_delta ** 2)
                o_pm_norm = gets_speed(o_pm, speeds)

                # Saca el catalog, alpha y delta para cada instante del objeto
                # df_source = data_full[data_full['SOURCE_NUMBER'].isin([source_])]
                df_source = df_filter[df_filter['SOURCE_NUMBER'].isin([source_])]
                for idx, row in df_source.iterrows():
                    catalog_n = int(row.CATALOG_NUMBER)
                    # Saca el alpha y el delta de cada objeto
                    alpha = float(row.ALPHA_J2000)
                    delta = float(row.DELTA_J2000)
                    # Busca si el objeto es un sso o no
                    i_df = search(df_input, catalog_n, alpha, delta)
                    if i_df.empty is not True:
                        tmp_flag.append('True')
                        i_pm = float(i_df['pm_values'].iloc[0])
                        tmp_pm_in.append(i_pm)
                    else:
                        tmp_flag.append('False')

                if o_pm_norm != 0.0:  # scamp has detected some movement
                    idx = values.index(o_pm_norm)
                    stats_d['N_meas'][idx] += 1
                    # checks if dataframes for each dither are not empty
                    flag_right, flag_number_right = all_same(tmp_flag)
                    if flag_right:  # flag True - object is an SSO
                        # checks if PM in all input dataframes are equal
                        flag_pm, flag_number_pm = all_same(tmp_pm_in)
                        if flag_pm:
                            idx = values.index(o_pm_norm)  # gets idx for extracted pm
                            if tmp_pm_in[0] == o_pm_norm:  # right pm
                                stats_d['N_se'][idx] += 1
                                flag = True
                            else:  # wrong pm
                                stats_d['N_meas-N_se'][idx] += 1
                                flag = False
                        else:
                            pass  # fixme strange error, should raise something
                    else:  # flag False - object is not an SSO
                        idx = values.index(o_pm_norm) # gets idx for extracted pm
                        stats_d['N_meas-N_se'][idx] += 1
                else:
                    pass  # star!

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
        # N_meas: number of all detected sources(including false detections)
        # N_se: number of simulated sources recovered by source extraction
        # N_true: number of simulated input sources
        # f_dr: detection rate f_pur: purity
        # f_com: completeness

        # f_dr = N_meas / N_true = (N_se + N_false) / N_true
        # f_pur = N_se / N_meas = N_se / (N_se + N_false)
        # f_com = f_dr * f_pur = N_se / N_true

        df_stats = DataFrame(stats_d,
                             columns=['confidence', 'cross_id',
                                      'pixscale_maxerr', 'position_maxerr',
                                      'posangle_maxerr', 'deblending',
                                      'mincount', 'threshold', 'area', 'mag',
                                      'pm', 'N_true', 'N_se', 'N_meas-N_se',
                                      'N_meas', 'f_dr', 'f_pur', 'f_com'])
        df_stats.to_csv('stats_{}_{}_{}.csv'.format(self.sex_cf, self.scmp_cf,
                                                    self.mag))

    def speeds_distribution(self, df_input):
        """

        @param df_input:

        @return total_pms:
        """
        # Creates a new DataFrame just to get ssos' speed distribution
        # With this dataframe we can know how objects' speed is distribuited
        ssos_created = concat(g for _,
                              g in df_input.groupby('source') if len(g) >= 3)

        unique_ssos = list(set(ssos_created['source'].tolist()))
        total_pms = []

        for source_ in unique_ssos:
            df_ssos = ssos_created[ssos_created['source'].isin([source_])]
            total_pms.append(float(df_ssos['pm_values'].iloc[0]))

        return total_pms


def merge_stats(logger, prfs_d):
    """ Merges stats file by confidence, pm and magnitude values

    @param prfs_d:

    @return True: if everything goes alright.
    """

    list_files = listdir(getcwd() + '/stats/')
    stats_d = {}

    # TODO Improve description
    # Merge stats files by confidence and pm
    for mag_ in prfs_d['mags']:
        for confidence_ in prfs_d['confidences']:
            for pm_ in prfs_d['pms']:
                conf = '{}_{}_{}'.format(mag_, confidence_, pm_)
                stats_d[conf] = []
                for file_ in list_files:
                    if 'stats_' in file_ and '_.csv' in file_:
                        try:
                            f_df = read_csv(getcwd() + '/stats/' + file_,
                                            index_col=0)
                            df = f_df[f_df['pm'].isin([pm_])]
                            df = df[df['confidence'].isin([confidence_])]
                            stats_d[conf].append(df)
                        except Exception as e:
                            print e, file_

    for mag_ in prfs_d['mags']:
        for confidence_ in prfs_d['confidences']:
            for pm_ in prfs_d['pms']:
                conf = '{}_{}_{}'.format(mag_, confidence_, pm_)
                concated_df = concat(stats_d[conf])
                if not path.exists('{}'.format(confidence_)):
                    makedirs('{}'.format(confidence_))
                concated_df.to_csv('{}/{}.csv'.format(confidence_, conf))

    return True


def create_dicts():
    """

    """
    filt_keys = ['CCD', 'i_alpha', 'i_delta', 'o_alpha', 'o_delta',
                 'm_alpha', 'm_delta', 'dither', 'pm_input', 'pm_output',
                 'source_scamp', 'source_merged', 'dispersion', 'confidence',
                 'source_input']
    filt_dict = {}
    for key_ in filt_keys:
        filt_dict[key_] = []

    pm_keys = ['cross_id', 'pixscale_maxerr', 'posangle_maxerr',
               'position_maxerr', 'confidence', 'CCD', 'alpha_source',
               'delta_source', 'alpha_detected', 'delta_detected',
               'alpha_difference', 'delta_difference', 'total_difference',
               'mag', 'dither', 'pm_input', 'pm_output', 'dispersion',
               'source_luca', 'source_scamp', 'flag_dispersion']
    pm_dict = {}
    for key_ in pm_keys:
        pm_dict[key_] = []

    stats_keys = ['cross_id', 'pixscale_maxerr', 'posangle_maxerr',
                  'position_maxerr', 'confidence', 'deblending', 'mincount',
                  'threshold', 'area','mag', 'pm', 'N_true', 'N_meas', 'N_se',
                  'N_meas-N_se', 'f_dr', 'f_pur', 'f_com']
    stats_dict = {}
    for key_ in stats_keys:
        stats_dict[key_] = []

    return filt_dict, pm_dict, stats_dict
