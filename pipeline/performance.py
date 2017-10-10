#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for

Todo:
    * Improve log messages

"""

from multiprocessing import Process
from os import listdir

from pandas import concat, read_csv

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
        logger = setting_logger()
        prfs_d = extract_settings()

        if not self.check(logger, prfs_d):
            raise Exception

    def check(self, logger, prfs_d):
        """

        @param logger:
        @param prfs_d:

        @return True: if everything goes alright.
        """
        confs = [2, 0.1, 5, 4, 'models/gauss_2.0_5x5.conv']

        analysis_d = {'deblend_mincount': 0.1,
                      'analysis_thresh': 5,
                      'detect_thresh': 5,
                      'deblend_nthresh': 2, 'detect_minarea': 4,
                      'filter': 'models/gauss_2.0_5x5.conv'}

        # Creates an input dictionary with all input sources
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

        # Cross with filtered data - Opens datafile
        filter_n = 'filt_10_1.2_5_0.033_20-21__5.csv'
        o_cat = read_csv('{}/{}'.format(prfs_d['filter_dir'], filter_n),
                                        index_col=0)

        stats_d = self.create_dict()

        tmp_sources_i = []
        tmp_sources_l = []

        unique_sources = list(set(i_df['source'].tolist()))

        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            cat = i_df[i_df['source'].isin([source_])]
            boolean_l = []

            for i, row in enumerate(cat.itertuples(), 1):
                catalog_n = row.catalog
                pm = row.pm_values

                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
                o_df = o_df[o_df['ALPHA_J2000'] + 0.001 > i_alpha]
                o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - 0.001]
                o_df = o_df[o_df['DELTA_J2000'] + 0.001> i_delta]
                o_df = o_df[i_delta > o_df['DELTA_J2000'] - 0.001]

                if o_df.empty is not True:
                    print catalog_n, "True"
                    boolean_l.append(True)
                else:
                    boolean_l.append(False)
            
            print boolean_l
            # Total number
            idx = stats_d['PM'].index(pm)
            stats_d['total'][idx] += 1

            if len(list(set(boolean_l))) == 1 and list(set(boolean_l))[0] == True:
                idx = stats_d['PM'].index(pm)
                stats_d['right'][idx] += 1
            else:
                pass

        print stats_d

        return True

    def create_dict(self):
        """

        """
        stats_keys = ['total', 'right', 'false']

        stats_d = {}
        stats_d['PM'] = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3,
                         1, 3, 10, 30, 100, 300]

        for key_ in stats_keys:
            stats_d[key_] = []
            for value_ in range(len(stats_d['PM'])):
                stats_d[key_].append(0)

        return stats_d


if __name__ == '__main__':
    performance = ScampPerformance()
