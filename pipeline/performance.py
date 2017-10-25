#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for

Versions:
- 0.2 Now supports confidence intervals

Todo:
    * Improve log messages

"""

from pandas import concat, read_csv, DataFrame

from misc import all_same, speeds_range
from regions import Create_regions


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
"""
__license__ = "GPL"
"""
__version__ = "0.2"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


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
        # i_df.to_csv('input_sources.csv')

        """
        alpha_df = i_df['alpha_j2000']
        delta_df = i_df['delta_j2000']

        df = concat([alpha_df, delta_df], axis=1)
        df.to_csv('input_sources.reg')
        """

        # Open particular file!
        filt_n = 'filt_{}_{}_4.csv'.format(scmp_cf, mag)
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
                if o_df.empty is not True:
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
        """ shows in red input sources
        shows in blue extracted ones

        @param i_alpha:
        @param i_delta
        """
        i_dict = {'i_alpha': i_alpha_l, 'i_delta': i_delta_l}
        i_df = DataFrame(i_dict)
        i_df.to_csv('{}.csv'.format(source_), index=False,
                    header=False, sep=" ")

    """
    def show_regions(self, input_fits, regions, num):


        cmd_11 = 'ds9 {} -zoom to fit -histequ '.format(input_fits)
        cmd_12 = '-regions -format xy '
        # load inputs regions
        cmd_13 = '{}pipeline/output_{}.reg '.format(prfs_d['fed_home'], num)
        cmd_14 = '-saveimage jpeg {} -exit'.format(output_image)
        cmd_1 = cmd_11 + cmd_12 + cmd_13 + cmd_14

        ds9_process = Popen(cmd_1, shell=True)
        ds9_process.wait()

        return True
    """

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
        stats_d['PM'] = [0.0003, 0.001, 0.003, 0.01, 0.03,
                         0.1, 0.3, 3, 10, 30, 100, 300]

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
