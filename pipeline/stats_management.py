#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * Improve log messages

"""

from collections import Counter
from math import hypot
from numpy import isclose
from os import listdir, getcwd, path, makedirs

from pandas import concat, DataFrame, read_csv

from misc import extract_settings, speeds_range

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
"""
__license__ = "GPL"
"""
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class ExtractStats:

    def __init__(self, logger, mag, scmp_d, f_conf, filt):
        """

        @param logger:
        @param mag:
        @param scmp_d:
        @param f_conf:
        @filt: a boolean variable

        """
        prfs_d = extract_settings()

        (df_filter, df_input,
         df_output, valid_sources) = self.extract_stats(logger, prfs_d,
                                                        mag, scmp_d, f_conf)

        if filt:
            self.filtered_stats(prfs_d, df_filter, df_input,
                                df_output, valid_sources, scmp_d, mag, f_conf)
        else:
            self.total_stats()

    def extract_stats(self, logger, prfs_d, mag, scmp_d, f_conf):
        """

        @param logger: a logger object.
        @param prfs_d: a dictionary with all script parameters.
        @param mag: magnitude gap to be analysed.
        @param scmp_d: a dictionary with sextractor/scamp parameters.

        @return stats_d:
        """
        # Gets all sources detected by scamp after filtering process
        # Opens csv file
        filtered_cat_n = '{}/filt_{}_{}__5.csv'.format(prfs_d['results_dir'],
                                                       f_conf, mag)
        logger.debug('reading filtered catalog {}'.format(filtered_cat_n))
        df_filter = read_csv(filtered_cat_n, index_col=0)
        df_filter = concat(g for _, g in df_filter.groupby("SOURCE_NUMBER")
                           if len(g) >= int(prfs_d['detections']))

        # Gets all unique sources after filter
        logger.debug('getting unique sources from post-filter scamp catalog')
        valid_sources = list(set(df_filter['SOURCE_NUMBER'].tolist()))

        # Gets all sources created from Luca's
        df_input = read_csv('{}_{}_ssos.csv'.format(f_conf, mag), index_col=0)

        # Gets all sources detected by sextractor and scamp
        # This is useful for know how many objects could be
        # detected without filter TODO Check!
        df_output = read_csv('{}_{}_sso_cat.csv'.format(f_conf, mag),
                             index_col=0)

        # TODO Check if it's necessary!
        # Gets all sources detected by scamp before filtering process.
        # Useful for obtained scamp calculated proper motion.
        # merged_n = '{}/merged_{}_1.cat'.format(prfs_d['results_dir'], f_conf)
        # merged_cat = fits.open(merged_n)
        # data_merged = Table(merged_cat[2].data).to_pandas()

        return (df_filter, df_input, df_output, valid_sources)

    def filtered_stats(self, prfs_d, df_filter, df_input, df_output,
                       valid_sources, scmp_d, mag, f_conf):
        """

        @param prfs_d:
        @param df_filter: a DataFrame
        @param df_input:
        @param df_output:
        @param valid_sources: a list with all sources found by sextractor and
        scamp and filtered.
        @param scmp_d:
        @param mag:
        @param f_conf:

        Important parameters:
        @df_filter: a DataFrame object which contains all ssos filtered.
        @df_input: a DataFrame object which contains all ssos created.
        @df_output: a DataFrame object which contains all ssos detected by
                    sextractor and scamp. These sources were checked against
                    positions obtained from ssos.csv.
        """
        # Confidence values to be tested
        confidence_list = prfs_d['confidences']
        # Created dictionaries for all outputs
        filt_d, pm_d, stats_d = create_dicts()

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
            speeds = speeds_range(prfs_d, confidence_)
            empty_value = 0
            for value_ in values:
                stats_d['confidence'].append(confidence_)
                stats_d['position_maxerr'].append(scmp_d['position_maxerr'])
                stats_d['posangle_maxerr'].append(scmp_d['posangle_maxerr'])
                stats_d['pixscale_maxerr'].append(scmp_d['pixscale_maxerr'])
                stats_d['cross_id'].append(scmp_d['crossid_radius'])
                stats_d['mag'].append(mag)
                stats_d['pm'].append(value_)
                stats_d['input_number'].append(total_values[value_])
                stats_d['total_filter'].append(empty_value)
                stats_d['right_filter'].append(empty_value)
                # stats_d['f_dr'].append(empty_value)
                # stats_d['f_com'].append(empty_value)

            # valid_sources comes from post-filter
            for i, source_ in enumerate(valid_sources):
                # Source presents at filter's output
                str_source = str(source_)
                int_source = int(source_)

                # A ver, en df_output tengo las fuentes de entrada asociadas a
                # las fuentes de salida. Asi que busco si la source de mi
                # estudio se encuentra alli
                # TODO Change name at input file!!
                d_o = df_output[df_output['scamp_source_number'].isin([str_source])]
                # If it's empty there is not a real SSO in post-filter scamp source
                # so this sources goes to total sources but not to real ones

                # Source is a real SSO
                # Obtengo la velocidad que he obtenido del catalogo post-filter
                d_f = df_filter[df_filter['SOURCE_NUMBER'].isin([int_source])]
                pm_output = float(d_f['PM'].iloc[0])
                # y el valor de la fuente segun el catalogo de luca

                if d_o.empty is not True:
                    for ccd_ in range(0, 3, 1):
                        position_maxerr = scmp_d['position_maxerr']
                        pm_d['position_maxerr'].append(position_maxerr)
                        posangle_maxerr = scmp_d['posangle_maxerr']
                        pm_d['posangle_maxerr'].append(posangle_maxerr)
                        pixscale_maxerr = scmp_d['pixscale_maxerr']
                        pm_d['pixscale_maxerr'].append(pixscale_maxerr)
                        crossid_radius = scmp_d['crossid_radius']
                        pm_d['cross_id'].append(crossid_radius)
                        pm_d['mag'].append(mag)

                        # Source values
                        pm_d['source_scamp'].append(int_source)
                        luca_source = int(d_o['sources'].iloc[0])
                        pm_d['source_luca'].append(luca_source)

                        # Obtengo la velocidad asociada
                        # Para ello utilizo el catalogo ssos.csv
                        d_i = df_input[df_input['source'].isin([luca_source])]
                        alpha_source = float(d_i['alpha_j2000'].iloc[ccd_])
                        pm_d['alpha_source'].append(alpha_source)
                        delta_source = float(d_i['delta_j2000'].iloc[ccd_])
                        pm_d['delta_source'].append(delta_source)

                        # pm_input, pm_output
                        pm_input = float(d_i['pm_values'].iloc[0])
                        pm_d['pm_input'].append(pm_input)
                        pm_d['pm_output'].append(pm_output)

                        # alpha, delta, present at filter output
                        alpha_detected = float(d_f['ALPHA_J2000'].iloc[ccd_])
                        pm_d['alpha_detected'].append(alpha_detected)
                        delta_detected = float(d_f['DELTA_J2000'].iloc[ccd_])
                        pm_d['delta_detected'].append(delta_detected)

                        # Difference between input/output
                        alpha_difference = alpha_source - alpha_detected
                        alpha_difference = alpha_difference * 3600  # to secs
                        pm_d['alpha_difference'].append(alpha_difference)
                        delta_difference = delta_source - delta_detected
                        delta_difference = delta_difference * 3600  # to secs
                        pm_d['delta_difference'].append(delta_difference)
                        total_difference = hypot(alpha_difference,
                                                 delta_difference)
                        pm_d['total_difference'].append(total_difference)

                        ccd = d_i['CCD'].iloc[ccd_]
                        pm_d['CCD'].append(ccd)
                        dither = d_i['dither_values'].iloc[ccd_]
                        pm_d['dither'].append(dither)

                    # Once proper motions are obtained dispersion should be
                    # calculated
                    if pm_input < pm_output:
                        dispersion = (pm_output / pm_input) * 100
                    elif pm_input > pm_output:
                        dispersion = (pm_output / pm_input) * 100

                    for ccd_ in range(0, 3, 1):
                        pm_d['confidence'].append(confidence_)
                        pm_d['dispersion'].append(dispersion)

                    # If dispersion is lower than confidence the SSO is right
                    if isclose(dispersion, 100, atol=confidence_):
                        for ccd_ in range(0, 3, 1):
                            pm_d['flag_dispersion'].append(False)
                        idx = stats_d['pm'].index(pm_input)
                        idx = idx + (12 * confidence_list.index(confidence_))
                        stats_d['right_filter'][idx] += 1
                        stats_d['total_filter'][idx] += 1
                    # else the object is process as a false positive
                    else:
                        for ccd_ in range(0, 3, 1):
                            pm_d['flag_dispersion'].append(True)
                        # TODO As these lines are equal to next ones,
                        # a new function should be created?
                        # There is a "null" pm for cases that cannot be
                        # organizados
                        pm_normalized = 0
                        for key_ in speeds.keys():
                            # print "fuera", type(speeds[key_][0]), pm_output
                            # print "fuera_2", speeds[key_][0], speeds[key_][1]
                            if speeds[key_][0] < pm_output < speeds[key_][1]:
                                pm_normalized = key_
                                # print "dentro", pm_normalized, key_

                        if pm_normalized is not 0:
                            idx = stats_d['pm'].index(pm_normalized)
                            # idx_tmp should be changing across results
                            idx_tmp = 12 * confidence_list.index(confidence_)
                            idx = idx + idx_tmp
                            stats_d['total_filter'][idx] += 1

                elif d_o.empty is True:
                    # Source is not an SSO
                    # Obten la velocidad que te dio para ver donde lo metes
                    # de errores vamos
                    # There is a "null" pm for cases that cannot be
                    # organizados
                    pm_normalized = 0
                    for key_ in speeds.keys():
                        # print "fuera", type(speeds[key_][0]), pm_output
                        # print "fuera_2", speeds[key_][0], speeds[key_][1]
                        if speeds[key_][0] < pm_output < speeds[key_][1]:
                            pm_normalized = key_
                            # print "dentro", pm_normalized, key_

                    if pm_normalized is not 0:
                        idx = stats_d['pm'].index(pm_normalized)
                        idx = idx + (12 * confidence_list.index(confidence_))
                        stats_d['total_filter'][idx] += 1

        for i, value_ in enumerate(stats_d['input_number']):
            total_number = float(stats_d['total_filter'][i])
            input_number = float(stats_d['input_number'][i])
            right_number = float(stats_d['right_filter'][i])
            f_dr = total_number / input_number
            stats_d['f_dr'].append(f_dr)
            f_com = right_number / input_number
            stats_d['f_com'].append(f_com)

        df_stats = DataFrame(stats_d)
        df_stats.to_csv('stats/stats_{}.csv'.format(f_conf))

        df_pm = DataFrame(pm_d)
        df_pm.to_csv('stats/pm_{}.csv'.format(f_conf))

    def total_stats(self, prfs_d, df_filter, df_input, df_output,
                    valid_sources, scmp_d, mag, f_conf):
        """

        """
        # Looks for sources one by one in sextractor/scamp catalog
        # Try to find if they're an SSO or not
        # source_ number refered to custom catalog, not sextractor/scamp one
        # Confidence values, how pm can be away from real value
        # Confidence values to be tested

        confidence_list = prfs_d['confidences']
        # Created dictionaries for all outputs
        filt_d, pm_d, stats_d = create_dicts()

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
            speeds = speeds_range(prfs_d, confidence_)
            empty_value = 0
            for value_ in values:
                stats_d['confidence'].append(confidence_)
                stats_d['position_maxerr'].append(scmp_d['position_maxerr'])
                stats_d['posangle_maxerr'].append(scmp_d['posangle_maxerr'])
                stats_d['pixscale_maxerr'].append(scmp_d['pixscale_maxerr'])
                stats_d['cross_id'].append(scmp_d['crossid_radius'])
                stats_d['mag'].append(mag)
                stats_d['pm'].append(value_)
                stats_d['input_number'].append(total_values[value_])
                stats_d['total_filter'].append(empty_value)
                stats_d['right_filter'].append(empty_value)
                # stats_d['f_dr'].append(empty_value)
                # stats_d['f_com'].append(empty_value)

            # valid_sources comes from post-filter
            for i, source_ in enumerate(valid_sources):
                # Source presents at filter's output
                str_source = str(source_)
                int_source = int(source_)

                # A ver, en df_output tengo las fuentes de entrada asociadas a
                # las fuentes de salida. Asi que busco si la source de mi
                # estudio se encuentra alli
                # TODO Change name at input file!!
                d_o = df_output[df_output['scamp_source_number'].isin([str_source])]
                # If it's empty there is not a real SSO in post-filter scamp source
                # so this sources goes to total sources but not to real ones

                # Source is a real SSO
                # Obtengo la velocidad que he obtenido del catalogo post-filter
                d_f = df_filter[df_filter['SOURCE_NUMBER'].isin([int_source])]
                pm_output = float(d_f['PM'].iloc[0])
                # y el valor de la fuente segun el catalogo de luca

                if d_o.empty is not True:
                    for ccd_ in range(0, 3, 1):
                        position_maxerr = scmp_d['position_maxerr']
                        pm_d['position_maxerr'].append(position_maxerr)
                        posangle_maxerr = scmp_d['posangle_maxerr']
                        pm_d['posangle_maxerr'].append(posangle_maxerr)
                        pixscale_maxerr = scmp_d['pixscale_maxerr']
                        pm_d['pixscale_maxerr'].append(pixscale_maxerr)
                        crossid_radius = scmp_d['crossid_radius']
                        pm_d['cross_id'].append(crossid_radius)
                        pm_d['mag'].append(mag)

                        # Source values
                        pm_d['source_scamp'].append(int_source)
                        luca_source = int(d_o['sources'].iloc[0])
                        pm_d['source_luca'].append(luca_source)

                        # Obtengo la velocidad asociada
                        # Para ello utilizo el catalogo ssos.csv
                        d_i = df_input[df_input['source'].isin([luca_source])]
                        alpha_source = float(d_i['alpha_j2000'].iloc[ccd_])
                        pm_d['alpha_source'].append(alpha_source)
                        delta_source = float(d_i['delta_j2000'].iloc[ccd_])
                        pm_d['delta_source'].append(delta_source)

                        # pm_input, pm_output
                        pm_input = float(d_i['pm_values'].iloc[0])
                        pm_d['pm_input'].append(pm_input)
                        pm_d['pm_output'].append(pm_output)

                        # alpha, delta, present at filter output
                        alpha_detected = float(d_f['ALPHA_J2000'].iloc[ccd_])
                        pm_d['alpha_detected'].append(alpha_detected)
                        delta_detected = float(d_f['DELTA_J2000'].iloc[ccd_])
                        pm_d['delta_detected'].append(delta_detected)

                        ccd = d_i['CCD'].iloc[ccd_]
                        pm_d['CCD'].append(ccd)
                        dither = d_i['dither_values'].iloc[ccd_]
                        pm_d['dither'].append(dither)

                    # Once proper motions are obtained dispersion should be
                    # calculated
                    if pm_input < pm_output:
                        dispersion = pm_output / pm_input
                    elif pm_input > pm_output:
                        dispersion = pm_output / pm_input
                    else:
                        raise Exception  # TODO what value will be get?

                    for ccd_ in range(0, 3, 1):
                        pm_d['confidence'].append(confidence_)
                        pm_d['dispersion'].append(dispersion)

                    # If dispersion is lower than confidence the SSO is right
                    if dispersion < confidence_:
                        for ccd_ in range(0, 3, 1):
                            pm_d['flag_dispersion'].append(False)
                        idx = stats_d['pm'].index(pm_input)
                        idx = idx + (12 * confidence_list.index(confidence_))
                        stats_d['right_filter'][idx] += 1
                        stats_d['total_filter'][idx] += 1
                    # else the object is process as a false positive
                    else:
                        for ccd_ in range(0, 3, 1):
                            pm_d['flag_dispersion'].append(True)
                        # TODO As these lines are equal to next ones,
                        # a new function should be created?
                        # There is a "null" pm for cases that cannot be
                        # organizados
                        pm_normalized = 0
                        for key_ in speeds.keys():
                            # print "fuera", type(speeds[key_][0]), pm_output
                            # print "fuera_2", speeds[key_][0], speeds[key_][1]
                            if speeds[key_][0] < pm_output < speeds[key_][1]:
                                pm_normalized = key_
                                # print "dentro", pm_normalized, key_

                        if pm_normalized is not 0:
                            idx = stats_d['pm'].index(pm_normalized)
                            # idx_tmp should be changing across results
                            idx_tmp = 12 * confidence_list.index(confidence_)
                            idx = idx + idx_tmp
                            stats_d['total_filter'][idx] += 1

                elif d_o.empty is True:
                    # Source is not an SSO
                    # Obten la velocidad que te dio para ver donde lo metes
                    # de errores vamos
                    # There is a "null" pm for cases that cannot be
                    # organizados
                    pm_normalized = 0
                    for key_ in speeds.keys():
                        # print "fuera", type(speeds[key_][0]), pm_output
                        # print "fuera_2", speeds[key_][0], speeds[key_][1]
                        if speeds[key_][0] < pm_output < speeds[key_][1]:
                            pm_normalized = key_
                            # print "dentro", pm_normalized, key_

                    if pm_normalized is not 0:
                        idx = stats_d['pm'].index(pm_normalized)
                        idx = idx + (12 * confidence_list.index(confidence_))
                        stats_d['total_filter'][idx] += 1
                else:
                    raise Exception

        for i, value_ in enumerate(stats_d['input_number']):
            total_number = float(stats_d['total_filter'][i])
            input_number = float(stats_d['input_number'][i])
            right_number = float(stats_d['right_filter'][i])
            f_dr = total_number / input_number
            stats_d['f_dr'].append(f_dr)
            f_com = right_number / input_number
            stats_d['f_com'].append(f_com)

        df_stats = DataFrame(stats_d)
        df_stats.to_csv('stats/stats_{}_.csv'.format(f_conf))

        df_pm = DataFrame(pm_d)
        df_pm.to_csv('stats/pm_{}_.csv'.format(f_conf))

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
                  'position_maxerr', 'confidence', 'mag', 'pm', 'input_number',
                  'total_filter', 'right_filter', 'f_dr', 'f_com']
    stats_dict = {}
    for key_ in stats_keys:
        stats_dict[key_] = []

    return filt_dict, pm_dict, stats_dict
