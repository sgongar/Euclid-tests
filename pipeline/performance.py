#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for

Todo:
    * Improve log messages

"""

from multiprocessing import Process
from os import listdir
from subprocess import Popen

from pandas import concat, read_csv, DataFrame

from misc import create_configurations, get_cats, extract_settings
from misc import setting_logger, check_distance
from regions import Create_regions


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


class ScampPerformance:

    def __init__(self):
        pass

    def check(self, logger, prfs_d, mag, scmp_cf, sex_cf, idx_file):
        """

        @param logger:
        @param prfs_d:

        @return True: if everything goes alright.
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

        # Open particular file!

        filt_n = 'filt_{}_{}_5.csv'.format(scmp_cf, mag)
        filter_o_n = '{}/{}/{}/{}'.format(prfs_d['filter_dir'],
                                          sex_cf, scmp_cf, filt_n)

        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)
        stats_d, out_d = self.create_dict()

        unique_sources = list(set(i_df['source'].tolist()))

        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            cat = i_df[i_df['source'].isin([source_])]
            boolean_l = []
            tmp_catalog = []
            tmp_pm = []
            tmp_alpha = []
            tmp_delta = []

            # Iterate over each source
            for i, row in enumerate(cat.itertuples(), 1):
                catalog_n = row.catalog
                pm = row.pm_values

                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # Checks if there is a source closed to input one
                o_df = self.check_source(catalog_n, o_cat, i_alpha, i_delta)

                if o_df.empty is not True:
                    boolean_l.append(True)
                    tmp_catalog.append(catalog_n)
                    tmp_pm.append(pm)
                    tmp_alpha.append(i_alpha)
                    tmp_delta.append(i_delta)

                    print "True"
                else:
                    # Appends alpha and delta to a temp list
                    # Just to check what happens in that area

                    # Create a big region file with all faults
                    out_d['catalog'].append(catalog_n)
                    out_d['PM'].append(pm)
                    out_d['alpha_j2000'].append(i_alpha)
                    out_d['delta_j2000'].append(i_delta)

                    boolean_l.append(False)

            # Total number
            idx = stats_d['PM'].index(pm)
            stats_d['total'][idx] += 1

            if len(list(set(boolean_l))) == 1 and list(set(boolean_l))[0] == True:
                idx = stats_d['PM'].index(pm)
                stats_d['right'][idx] += 1
            else:
                print "catalog", tmp_catalog
                print "pm", tmp_pm
                print "tmp_alpha", tmp_alpha
                print "tmp_delta", tmp_delta
                # self.create_regions(tmp_alpha, tmp_delta, source_)
                # self.show_regions()
                pass

        """
        stats_df = DataFrame(stats_d)
        stats_df.to_csv('stats_{}.csv'.format(idx_file))

        out_df = DataFrame(out_d)
        # Saves output to easily readable files
        out_df.to_csv('errors_{}.csv'.format(idx_file))
        out_df.to_csv('errors_{}.reg'.format(idx_file),
                      index=False, header=False, sep=" ")
        """

        return stats_d

    def check_source(self, catalog_n, o_cat, i_alpha, i_delta):
        """

        @param catalog_n:
        @param o_cat:
        @param i
        """
        o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
        o_df = o_df[o_df['ALPHA_J2000'] + 0.001 > i_alpha]
        o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - 0.001]
        o_df = o_df[o_df['DELTA_J2000'] + 0.001 > i_delta]
        o_df = o_df[i_delta > o_df['DELTA_J2000'] - 0.001]

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

    def create_dict(self):
        """

        """
        stats_keys = ['total', 'right', 'false',
                      'f_dr', 'f_pur', 'f_com']

        stats_d = {}
        stats_d['PM'] = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3,
                         1, 3, 10, 30, 100, 300]

        for key_ in stats_keys:
            stats_d[key_] = []
            for value_ in range(len(stats_d['PM'])):
                stats_d[key_].append(0)

        out_keys = ['alpha_j2000', 'delta_j2000',
                    'catalog', 'PM']
        out_d = {}

        for key_ in out_keys:
            out_d[key_] = []

        return stats_d, out_d


if __name__ == '__main__':
    performance = ScampPerformance()
