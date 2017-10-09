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
from misc import setting_logger
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

        input_df = concat(input_l, axis=0)
        # Look for < 3 coincidences
        input_df = concat(g for _, g in input_df.groupby('source')
                          if len(g) >= 3)
        input_df = input_df.reset_index(drop=True)
        # input_df.to_csv('test.csv')  # Saves DataFrame to csv file

        # Cross with filtered data - Opens datafile
        filter_n = 'filt_10_1.2_5_0.033_20-21__5.csv'
        filter_cat = read_csv('{}/{}'.format(prfs_d['filter_dir'], filter_n),
                              index_col=0)

        for source_ in filter_cat['SOURCE_NUMBER'].tolist():
            cat = filter_cat[filter_cat['SOURCE_NUMBER'].isin([source_])]

            print cat
            """
            for cat in range(0, 3, 1):
                cat_number = 
            """

        return True


if __name__ == '__main__':
    performance = ScampPerformance()
