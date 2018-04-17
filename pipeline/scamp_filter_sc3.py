#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Versions:
- 0.1: Initial release.

Todo:
    * Improve documentation
    * Speed-up areas process

*GNU Terry Pratchett*
"""

from multiprocessing import Process

from astropy.io import fits
from astropy.table import Table
from numpy import mean, median
from pandas import concat, read_csv, Series

from misc import pm_compute, extract_settings_sc3, check_source
from misc import b_filter, create_folder, get_dither, confidence_filter

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class ScampFilterSC3:  # TODO Split scamp_filter method into single methods

    def __init__(self, logger, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_d:
        """
        print('scamp_filter')
        # Filter variables
        self.class_star_limit = 0.97
        self.proper_motion = 1.5
        self.proper_motion_dects = 1.5

        # Analysis variables
        self.prfs_d = extract_settings_sc3()
        self.logger = logger

        self.save = True

        # Filtered catalog dir
        self.filt_n = 'filt_'
        self.filter_o_n = '{}/{}'.format(self.prfs_d['filtered'], self.filt_n)

        # Saves _1.csv
        (merged_db, full_db) = self.scamp_filter()
        # Saves _2.csv
        self.compute_pm(merged_db, full_db)
        # Saves _3.csv
        # self.get_areas()

        """
        full_db = read_csv('{}_3.csv'.format(self.filter_o_n), index_col=0)
        self.slow_db, self.fast_db = self.filter_class(full_db)
        if self.save:
            self.save_message('4s')
            self.slow_db.to_csv('{}_4s.csv'.format(self.filter_o_n))
            # self.save_message('4f')
            self.fast_db.to_csv('{}_4f.csv'.format(self.filter_o_n))

        filter_j = []
        for pipeline in ['slow', 'fast']:
            filter_p = Process(target=self.choose_pipeline, args=(pipeline,))
            filter_j.append(filter_p)
            filter_p.start()

        active_filter = list([job.is_alive() for job in filter_j])
        while True in active_filter:
            active_filter = list([job.is_alive() for job in filter_j])
            pass
        # Merges catalogs
        slow_db = read_csv('{}_6s.csv'.format(self.filter_o_n), index_col=0)
        fast_db = read_csv('{}_6f.csv'.format(self.filter_o_n), index_col=0)

        full_db = concat([slow_db, fast_db])

        slow_db = full_db[full_db['PM'] < self.proper_motion_dects]
        try:
            slow_db = self.filter_detections(slow_db, 3)
            if self.save:
                self.save_message('7s')
                slow_db.to_csv('{}_7s.csv'.format(self.filter_o_n))
        except ValueError:
            slow = False

        fast_db = full_db[full_db['PM'] > self.proper_motion_dects]
        fast_db = self.filter_detections(fast_db, 1)
        if self.save:
            self.save_message('7f')
            fast_db.to_csv('{}_7f.csv'.format(self.filter_o_n))

        if slow:
            full_db = concat([slow_db, fast_db])
        else:
            full_db = fast_db

        if self.save:
            self.save_message('8')
            full_db.to_csv('{}_8.csv'.format(self.filter_o_n))

        full_db = full_db[full_db['PM'] > 0.001]

        if self.save:
            self.save_message('9')
            full_db.to_csv('{}_9.csv'.format(self.filter_o_n))
        """

    def save_message(self, order):
        """

        :param order:
        :return:
        """
        self.logger.debug('Saves data to: ')
        # self.logger.debug('Dir: {}'.format(self.filter_dir))
        self.logger.debug('Name: {}_{}.csv'.format(self.filt_n, order))

    def get_cat(self, cat_n):
        """

        :param cat_n:
        :return: cat_file
        """
        """
        dither = ['1', '2', '3', '4']
        cat_part_2 = '000000-0000000'
        cat_times = ['20170630T011437.3Z', ]

        for d in dither:
            cats_final.append('{}-{}-{}__{}.cat'.format(cat_part_1, d,
                                                        cat_part_2,
                                                        cat_times[idx_time]))
        """
        cat_part_1 = 'EUC_VIS_SWL-DET-00'

        cats_final = []
        cats_partial = ['', '1-000000-0000000__20170630T011437.3Z_00.00_0',
                        '1-000000-0000000__20170630T011437.3Z_00.00_10',
                        '1-000000-0000000__20170630T011437.3Z_00.00_11',
                        '1-000000-0000000__20170630T011437.3Z_00.00_12',
                        '1-000000-0000000__20170630T011437.3Z_00.00_13',
                        '1-000000-0000000__20170630T011437.3Z_00.00_14',
                        '1-000000-0000000__20170630T011437.3Z_00.00_15',
                        '1-000000-0000000__20170630T011437.3Z_00.00_16',
                        '1-000000-0000000__20170630T011437.3Z_00.00_17',
                        '1-000000-0000000__20170630T011437.3Z_00.00_18',
                        '1-000000-0000000__20170630T011437.3Z_00.00_19',
                        '1-000000-0000000__20170630T011437.3Z_00.00_1',
                        '1-000000-0000000__20170630T011437.3Z_00.00_20',
                        '1-000000-0000000__20170630T011437.3Z_00.00_21',
                        '1-000000-0000000__20170630T011437.3Z_00.00_22',
                        '1-000000-0000000__20170630T011437.3Z_00.00_23',
                        '1-000000-0000000__20170630T011437.3Z_00.00_24',
                        '1-000000-0000000__20170630T011437.3Z_00.00_25',
                        '1-000000-0000000__20170630T011437.3Z_00.00_26',
                        '1-000000-0000000__20170630T011437.3Z_00.00_27',
                        '1-000000-0000000__20170630T011437.3Z_00.00_28',
                        '1-000000-0000000__20170630T011437.3Z_00.00_29',
                        '1-000000-0000000__20170630T011437.3Z_00.00_2',
                        '1-000000-0000000__20170630T011437.3Z_00.00_30',
                        '1-000000-0000000__20170630T011437.3Z_00.00_31',
                        '1-000000-0000000__20170630T011437.3Z_00.00_32',
                        '1-000000-0000000__20170630T011437.3Z_00.00_33',
                        '1-000000-0000000__20170630T011437.3Z_00.00_34',
                        '1-000000-0000000__20170630T011437.3Z_00.00_35',
                        '1-000000-0000000__20170630T011437.3Z_00.00_3',
                        '1-000000-0000000__20170630T011437.3Z_00.00_4',
                        '1-000000-0000000__20170630T011437.3Z_00.00_5',
                        '1-000000-0000000__20170630T011437.3Z_00.00_6',
                        '1-000000-0000000__20170630T011437.3Z_00.00_7',
                        '1-000000-0000000__20170630T011437.3Z_00.00_8',
                        '1-000000-0000000__20170630T011437.3Z_00.00_9',
                        '2-000000-0000000__20170630T011642.0Z_00.00_0',
                        '2-000000-0000000__20170630T011642.0Z_00.00_10',
                        '2-000000-0000000__20170630T011642.0Z_00.00_11',
                        '2-000000-0000000__20170630T011642.0Z_00.00_12',
                        '2-000000-0000000__20170630T011642.0Z_00.00_13',
                        '2-000000-0000000__20170630T011642.0Z_00.00_14',
                        '2-000000-0000000__20170630T011642.0Z_00.00_15',
                        '2-000000-0000000__20170630T011642.0Z_00.00_16',
                        '2-000000-0000000__20170630T011642.0Z_00.00_17',
                        '2-000000-0000000__20170630T011642.0Z_00.00_18',
                        '2-000000-0000000__20170630T011642.0Z_00.00_19',
                        '2-000000-0000000__20170630T011642.0Z_00.00_1',
                        '2-000000-0000000__20170630T011642.0Z_00.00_20',
                        '2-000000-0000000__20170630T011642.0Z_00.00_21',
                        '2-000000-0000000__20170630T011642.0Z_00.00_22',
                        '2-000000-0000000__20170630T011642.0Z_00.00_23',
                        '2-000000-0000000__20170630T011642.0Z_00.00_24',
                        '2-000000-0000000__20170630T011642.0Z_00.00_25',
                        '2-000000-0000000__20170630T011642.0Z_00.00_26',
                        '2-000000-0000000__20170630T011642.0Z_00.00_27',
                        '2-000000-0000000__20170630T011642.0Z_00.00_28',
                        '2-000000-0000000__20170630T011642.0Z_00.00_29',
                        '2-000000-0000000__20170630T011642.0Z_00.00_2',
                        '2-000000-0000000__20170630T011642.0Z_00.00_30',
                        '2-000000-0000000__20170630T011642.0Z_00.00_31',
                        '2-000000-0000000__20170630T011642.0Z_00.00_32',
                        '2-000000-0000000__20170630T011642.0Z_00.00_33',
                        '2-000000-0000000__20170630T011642.0Z_00.00_34',
                        '2-000000-0000000__20170630T011642.0Z_00.00_35',
                        '2-000000-0000000__20170630T011642.0Z_00.00_3',
                        '2-000000-0000000__20170630T011642.0Z_00.00_4',
                        '2-000000-0000000__20170630T011642.0Z_00.00_5',
                        '2-000000-0000000__20170630T011642.0Z_00.00_6',
                        '2-000000-0000000__20170630T011642.0Z_00.00_7',
                        '2-000000-0000000__20170630T011642.0Z_00.00_8',
                        '2-000000-0000000__20170630T011642.0Z_00.00_9',
                        '3-000000-0000000__20170630T011848.6Z_00.00_0',
                        '3-000000-0000000__20170630T011848.6Z_00.00_10',
                        '3-000000-0000000__20170630T011848.6Z_00.00_11',
                        '3-000000-0000000__20170630T011848.6Z_00.00_12',
                        '3-000000-0000000__20170630T011848.6Z_00.00_13',
                        '3-000000-0000000__20170630T011848.6Z_00.00_14',
                        '3-000000-0000000__20170630T011848.6Z_00.00_15',
                        '3-000000-0000000__20170630T011848.6Z_00.00_16',
                        '3-000000-0000000__20170630T011848.6Z_00.00_17',
                        '3-000000-0000000__20170630T011848.6Z_00.00_18',
                        '3-000000-0000000__20170630T011848.6Z_00.00_19',
                        '3-000000-0000000__20170630T011848.6Z_00.00_1',
                        '3-000000-0000000__20170630T011848.6Z_00.00_20',
                        '3-000000-0000000__20170630T011848.6Z_00.00_21',
                        '3-000000-0000000__20170630T011848.6Z_00.00_22',
                        '3-000000-0000000__20170630T011848.6Z_00.00_23',
                        '3-000000-0000000__20170630T011848.6Z_00.00_24',
                        '3-000000-0000000__20170630T011848.6Z_00.00_25',
                        '3-000000-0000000__20170630T011848.6Z_00.00_26',
                        '3-000000-0000000__20170630T011848.6Z_00.00_27',
                        '3-000000-0000000__20170630T011848.6Z_00.00_28',
                        '3-000000-0000000__20170630T011848.6Z_00.00_29',
                        '3-000000-0000000__20170630T011848.6Z_00.00_2',
                        '3-000000-0000000__20170630T011848.6Z_00.00_30',
                        '3-000000-0000000__20170630T011848.6Z_00.00_31',
                        '3-000000-0000000__20170630T011848.6Z_00.00_32',
                        '3-000000-0000000__20170630T011848.6Z_00.00_33',
                        '3-000000-0000000__20170630T011848.6Z_00.00_34',
                        '3-000000-0000000__20170630T011848.6Z_00.00_35',
                        '3-000000-0000000__20170630T011848.6Z_00.00_3',
                        '3-000000-0000000__20170630T011848.6Z_00.00_4',
                        '3-000000-0000000__20170630T011848.6Z_00.00_5',
                        '3-000000-0000000__20170630T011848.6Z_00.00_6',
                        '3-000000-0000000__20170630T011848.6Z_00.00_7',
                        '3-000000-0000000__20170630T011848.6Z_00.00_8',
                        '3-000000-0000000__20170630T011848.6Z_00.00_9',
                        '4-000000-0000000__20170630T012050.1Z_00.00_0',
                        '4-000000-0000000__20170630T012050.1Z_00.00_10',
                        '4-000000-0000000__20170630T012050.1Z_00.00_11',
                        '4-000000-0000000__20170630T012050.1Z_00.00_12',
                        '4-000000-0000000__20170630T012050.1Z_00.00_13',
                        '4-000000-0000000__20170630T012050.1Z_00.00_14',
                        '4-000000-0000000__20170630T012050.1Z_00.00_15',
                        '4-000000-0000000__20170630T012050.1Z_00.00_16',
                        '4-000000-0000000__20170630T012050.1Z_00.00_17',
                        '4-000000-0000000__20170630T012050.1Z_00.00_18',
                        '4-000000-0000000__20170630T012050.1Z_00.00_19',
                        '4-000000-0000000__20170630T012050.1Z_00.00_1',
                        '4-000000-0000000__20170630T012050.1Z_00.00_20',
                        '4-000000-0000000__20170630T012050.1Z_00.00_21',
                        '4-000000-0000000__20170630T012050.1Z_00.00_22',
                        '4-000000-0000000__20170630T012050.1Z_00.00_23',
                        '4-000000-0000000__20170630T012050.1Z_00.00_24',
                        '4-000000-0000000__20170630T012050.1Z_00.00_25',
                        '4-000000-0000000__20170630T012050.1Z_00.00_26',
                        '4-000000-0000000__20170630T012050.1Z_00.00_27',
                        '4-000000-0000000__20170630T012050.1Z_00.00_28',
                        '4-000000-0000000__20170630T012050.1Z_00.00_29',
                        '4-000000-0000000__20170630T012050.1Z_00.00_2',
                        '4-000000-0000000__20170630T012050.1Z_00.00_30',
                        '4-000000-0000000__20170630T012050.1Z_00.00_31',
                        '4-000000-0000000__20170630T012050.1Z_00.00_32',
                        '4-000000-0000000__20170630T012050.1Z_00.00_33',
                        '4-000000-0000000__20170630T012050.1Z_00.00_34',
                        '4-000000-0000000__20170630T012050.1Z_00.00_35',
                        '4-000000-0000000__20170630T012050.1Z_00.00_3',
                        '4-000000-0000000__20170630T012050.1Z_00.00_4',
                        '4-000000-0000000__20170630T012050.1Z_00.00_5',
                        '4-000000-0000000__20170630T012050.1Z_00.00_6',
                        '4-000000-0000000__20170630T012050.1Z_00.00_7',
                        '4-000000-0000000__20170630T012050.1Z_00.00_8',
                        '4-000000-0000000__20170630T012050.1Z_00.00_9']

        for cat_ in cats_partial:
            cats_final.append('{}{}.cat'.format(cat_part_1, cat_))

        cat_file = cats_final[cat_n]

        return cat_file

    def scamp_filter(self):
        """

        :return:
        """
        self.logger.info("Filtering scamp's output")

        # create_folder(self.logger, self.filter_dir)

        # Getting full catalog
        # Full catalog name
        full_cat_n = 'full_1.cat'  # todo - hardoced
        full_n = '{}/{}'.format(self.prfs_d['output_cats'], full_cat_n)

        # Just for logging reasons
        self.logger.debug('Opens full catalog')
        self.logger.debug('Dir: {}'.format(self.prfs_d['output_cats']))
        self.logger.debug('Name: {}'.format(full_cat_n))
        full_cat = fits.open(full_n)  # Opens full catalog
        full_db = Table(full_cat[2].data)  # Converts it to Astropy Table
        self.logger.debug('Converts full Astropy catalog to Pandas format')
        full_db = full_db.to_pandas()  # Converts it to Pandas format

        # Getting merge catalog
        # Merged catalog name
        mrgd_cat_n = 'merged_1.cat'
        mrgd_n = '{}/{}'.format(self.prfs_d['output_cats'], mrgd_cat_n)

        # Just for logging reasons
        self.logger.debug('Opens merged catalog')
        self.logger.debug('Dir: {}'.format(self.prfs_d['references']))
        self.logger.debug('Name: merged_1.cat')  # Be careful, hardcoded!
        merged_cat = fits.open(mrgd_n)  # Opens merged catalog
        self.logger.debug('Converts merged Astropy catalog to Pandas format')
        # todo - checks if Pandas format gonna be useful!
        merged_db = Table(merged_cat[2].data)

        merged_db = merged_db.to_pandas()

        # Removing 0 catalog detections
        self.logger.debug('Removes 0 catalog detections')
        full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

        full_db = self.filter_detections(full_db, 3)

        if self.save:
            self.save_message('1')
            full_db.to_csv('{}_1.csv'.format(self.filter_o_n))

        return merged_db, full_db

    def compute_pm(self, merged_db, full_db):
        """

        :param merged_db:
        :param full_db:
        :return: full_db
        """
        # Computing pm
        self.logger.debug('Computes proper motion')
        full_db = pm_compute(self.logger, merged_db, full_db)
        if self.save:
            self.save_message('2')
            full_db.to_csv('{}_2.csv'.format(self.filter_o_n))

    def choose_pipeline(self, pipeline):
        """

        :param pipeline:
        :return:
        """
        if pipeline == 'slow':
            # Starts pipeline for slow objects
            self.slow_pipeline()
        elif pipeline == 'fast':
            # Stars pipeline for fast objects
            self.fast_pipeline()
        else:
            raise Exception

    def slow_pipeline(self):
        """

        :return:
        """
        slow_db = self.slow_db[self.slow_db['PM'] < self.proper_motion]
        if self.save:
            self.save_message('5s')
            slow_db.to_csv('{}_5s.csv'.format(self.filter_o_n))

        slow_db = self.filter_b(slow_db)
        slow_db = slow_db[slow_db['A_IMAGE'] < 2]
        if self.save:
            self.save_message('6s')
            slow_db.to_csv('{}_6s.csv'.format(self.filter_o_n))

    def fast_pipeline(self):
        """

        :return:
        """
        fast_db = self.fast_db[self.fast_db['PM'] > self.proper_motion]
        if self.save:
            self.save_message('5f')
            fast_db.to_csv('{}_5f.csv'.format(self.filter_o_n))
        fast_db = self.filter_coherence(fast_db)
        if self.save:
            self.save_message('6f')
            fast_db.to_csv('{}_6f.csv'.format(self.filter_o_n))

    def get_areas(self):
        """

        :return: full_db
        """
        self.logger.debug('Populates filtered catalog with Sextractor data')

        # Opens filtered file
        filter_cat = read_csv('{}_2.csv'.format(self.filter_o_n), index_col=0)

        # Gets unique sources from filtered file
        unique_sources = list(set(filter_cat['SOURCE_NUMBER'].tolist()))
        l_sourcs = len(unique_sources)  # Just to not break 79 characters
        self.logger.debug('Unique sources to be analysed {}'.format(l_sourcs))

        dict_keys = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'EXTENSION',
                     'ASTR_INSTRUM', 'PHOT_INSTRUM', 'X_IMAGE', 'Y_IMAGE',
                     'ISOAREA_IMAGE', 'A_IMAGE', 'MEDIAN_A_IMAGE',
                     'MEAN_A_IMAGE', 'ERRA_IMAGE',  'MEDIAN_ERRA_IMAGE',
                     'MEAN_ERRA_IMAGE', 'B_IMAGE', 'MEDIAN_B_IMAGE',
                     'MEAN_B_IMAGE', 'ERRB_IMAGE', 'MEDIAN_ERRB_IMAGE',
                     'MEAN_ERRB_IMAGE', 'THETA_IMAGE', 'ERRTHETA_IMAGE',
                     'ALPHA_J2000', 'DELTA_J2000', 'ERRA_WORLD', 'ERRB_WORLD',
                     'ERRTHETA_WORLD', 'EPOCH', 'FWHM_IMAGE', 'CLASS_STAR',
                     'MEDIAN_CLASS_STAR', 'MEAN_CLASS_STAR', 'FLUX_ISO',
                     'MEDIAN_FLUX_ISO', 'MEAN_FLUX_ISO', 'FLUXERR_ISO',
                     'MEDIAN_FLUXERR_ISO', 'MEAN_FLUXERR_ISO', 'FLUX_RADIUS',
                     'ELONGATION', 'ELLIPTICITY', 'MEDIAN_ELLIPTICITY',
                     'MEAN_ELLIPTICITY', 'MAG', 'MAGERR', 'MAG_ISO',
                     'MEDIAN_MAG_ISO', 'MEAN_MAG_ISO', 'MAGERR_ISO',
                     'MEDIAN_MAGERR_ISO', 'MEAN_MAGERR_ISO',
                     'FLAGS_EXTRACTION', 'FLAGS_SCAMP', 'FLAGS_IMA', 'PM',
                     'PMERR', 'PMALPHA', 'PMDELTA', 'PMALPHAERR', 'PMDELTAERR']

        stats_keys = ['MEAN_A_IMAGE', 'MEAN_B_IMAGE', 'MEAN_CLASS_STAR',
                      'MEDIAN_A_IMAGE', 'MEDIAN_B_IMAGE', 'MEDIAN_CLASS_STAR',
                      'MEAN_ERRA_IMAGE', 'MEAN_ERRB_IMAGE', 'MEDIAN_ERRA_IMAGE',
                      'MEDIAN_ERRB_IMAGE', 'MEDIAN_FLUX_ISO', 'MEAN_FLUX_ISO',
                      'MEDIAN_FLUXERR_ISO', 'MEAN_FLUXERR_ISO',
                      'MEDIAN_ELLIPTICITY', 'MEAN_ELLIPTICITY',
                      'MEDIAN_MAG_ISO', 'MEAN_MAG_ISO', 'MEDIAN_MAGERR_ISO',
                      'MEAN_MAGERR_ISO']
        extra_keys = ['A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'ISOAREA_IMAGE',
                      'FWHM_IMAGE', 'FLUX_ISO', 'FLUXERR_ISO', 'FLUX_RADIUS',
                      'MAG_ISO', 'MAGERR_ISO', 'ELONGATION', 'ELLIPTICITY',
                      'CLASS_STAR']
        keys_l = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'EXTENSION',
                  'ASTR_INSTRUM', 'PHOT_INSTRUM', 'X_IMAGE', 'Y_IMAGE',
                  'ERRA_IMAGE', 'ERRB_IMAGE', 'ERRTHETA_IMAGE', 'ALPHA_J2000',
                  'DELTA_J2000', 'ERRA_WORLD', 'ERRB_WORLD', 'ERRTHETA_WORLD',
                  'EPOCH', 'MAG', 'MAGERR', 'FLAGS_EXTRACTION', 'FLAGS_SCAMP',
                  'FLAGS_IMA', 'PM', 'PMERR', 'PMALPHA', 'PMDELTA',
                  'PMALPHAERR', 'PMDELTAERR']

        sub_list_size = len(unique_sources) / self.prfs_d['cores_number']

        sub_list_l = []
        for idx_sub_list in range(0, self.prfs_d['cores_number'], 1):
            if idx_sub_list != (self.prfs_d['cores_number'] - 1):
                idx_down = sub_list_size * idx_sub_list
                idx_up = sub_list_size * (idx_sub_list + 1)
                sub_list_l.append(unique_sources[idx_down:idx_up])
            else:
                idx_down = sub_list_size * idx_sub_list
                sub_list_l.append(unique_sources[idx_down:])

        areas_j = []
        for idx_l in range(0, self.prfs_d['cores_number'], 1):
            areas_p = Process(target=self.get_areas_thread,
                              args=(dict_keys, stats_keys, extra_keys, keys_l,
                                    sub_list_l[idx_l], filter_cat, idx_l,))
            areas_j.append(areas_p)
            areas_p.start()

        active_areas = list([job.is_alive() for job in areas_j])
        while True in active_areas:
            active_areas = list([job.is_alive() for job in areas_j])
            pass

        # Merges areas
        # Merges catalogs
        csv_list = []
        for idx_csv in range(0, self.prfs_d['cores_number'], 1):
            csv_ = read_csv('{}_3_{}.csv'.format(self.filter_o_n, idx_csv),
                            index_col=0)
            csv_list.append(csv_)

        full_db = concat(csv_list)

        if self.save:
            self.save_message('3')
            full_db.to_csv('{}_3.csv'.format(self.filter_o_n))



        """
        sub_list_1_size = len(unique_sources) / 2
        sub_list_1 = unique_sources[:sub_list_1_size]
        sub_list_2 = unique_sources[sub_list_1_size:]
        sub_list_l = [sub_list_1, sub_list_2]

        areas_j = []
        for idx_l in range(0, 2, 1):
            areas_p = Process(target=self.get_areas_thread,
                              args=(dict_keys, stats_keys, extra_keys, keys_l,
                                    sub_list_l[idx_l], filter_cat, idx_l,))
            areas_j.append(areas_p)
            areas_p.start()

        active_areas = list([job.is_alive() for job in areas_j])
        while True in active_areas:
            active_areas = list([job.is_alive() for job in areas_j])
            pass

        # Merges areas
        # Merges catalogs
        list_1 = read_csv('{}_3_0.csv'.format(self.filter_o_n), index_col=0)
        list_2 = read_csv('{}_3_1.csv'.format(self.filter_o_n), index_col=0)

        full_db = concat([list_1, list_2])
        
        if self.save:
            self.save_message('3')
            full_db.to_csv('{}_3.csv'.format(self.filter_o_n))
        """

    def get_areas_thread(self, dict_keys, stats_keys, extra_keys, keys_l,
                         unique_sources_thread, filter_cat, idx_l):
        """

        :param dict_keys:
        :param stats_keys:
        :param extra_keys:
        :param keys_l:
        :param unique_sources_thread:
        :param filter_cat:
        :param idx_l:
        :return:
        """
        tmp_d = {}
        for key_ in keys_l:
            tmp_d[key_] = []

        for key_ in extra_keys:
            tmp_d[key_] = []

        for key_ in stats_keys:
            tmp_d[key_] = []

        cats_number = 144
        cat_d = {}
        for cat_n in range(1, cats_number + 1, 1):
            cat_dir = '/media/sf_CarpetaCompartida/luca_data/VIS_SC3/CCDs/cats/'
            cat_file = self.get_cat(cat_n)
            cat_data = fits.open('{}{}'.format(cat_dir, cat_file))

            ccd_df = Table(cat_data[2].data)
            self.logger.debug('converting CCD catalog to Pandas format')
            cat_d[cat_n] = ccd_df.to_pandas()

        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(unique_sources_thread):
            print(idx, idx_l)
            source_d = {'A_IMAGE': [], 'B_IMAGE': [], 'ERRA_IMAGE': [],
                        'ERRB_IMAGE': [], 'CLASS_STAR': [],
                        'FLUX_ISO': [], 'FLUXERR_ISO': [], 'FLUX_RADIUS': [],
                        'MAG_ISO': [], 'MAGERR_ISO': [], 'ELLIPTICITY': []}
            o_df = filter_cat[filter_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                # Populates temporal dictionary
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000
                cat_n = row.CATALOG_NUMBER
                # self.logger.debug('opening CCD catalog {}'.format(cat_file))

                cat_df = check_source(cat_d[cat_n], o_alpha, o_delta)
                if cat_df.empty:
                    # FIXME problemas entre sextractor y scamp
                    # Como el objeto no esta en las imagenes originales de
                    # sextractor lo borro del catalogo de scamp
                    # print('cat_empty')
                    for idx_k, key_ in enumerate(keys_l):
                        tmp_d[key_].append('nan')
                    for idx_k, key_ in enumerate(extra_keys):
                        tmp_d[key_].append('nan')
                else:
                    # print('cat not empty')
                    for idx_k, key_ in enumerate(keys_l):
                        tmp_d[key_].append(row[idx_k + 1])
                    for idx_k, key_ in enumerate(extra_keys):
                        tmp_d[key_].append(cat_df[key_].iloc[0])

                    source_d['A_IMAGE'].append(cat_df['A_IMAGE'].iloc[0])
                    source_d['B_IMAGE'].append(cat_df['B_IMAGE'].iloc[0])
                    source_d['ERRA_IMAGE'].append(cat_df['ERRA_IMAGE'].iloc[0])
                    source_d['ERRB_IMAGE'].append(cat_df['ERRB_IMAGE'].iloc[0])
                    source_d['CLASS_STAR'].append(cat_df['CLASS_STAR'].iloc[0])
                    source_d['FLUX_ISO'].append(cat_df['FLUX_ISO'].iloc[0])
                    source_d['FLUXERR_ISO'].append(cat_df['FLUXERR_ISO'].iloc[0])
                    source_d['FLUX_RADIUS'].append(cat_df['FLUX_RADIUS'].iloc[0])
                    source_d['MAG_ISO'].append(cat_df['MAG_ISO'].iloc[0])
                    source_d['MAGERR_ISO'].append(cat_df['MAGERR_ISO'].iloc[0])
                    source_d['ELLIPTICITY'].append(cat_df['ELLIPTICITY'].iloc[0])

            if 'nan' in source_d['A_IMAGE']:
                mean_a_image = 'nan'
                median_a_image = 'nan'
            else:
                mean_a_image = mean(source_d['A_IMAGE'])
                median_a_image = median(source_d['A_IMAGE'])

            if 'nan' in source_d['B_IMAGE']:
                mean_b_image = 'nan'
                median_b_image = 'nan'
            else:
                mean_b_image = mean(source_d['B_IMAGE'])
                median_b_image = median(source_d['B_IMAGE'])

            if 'nan' in source_d['ERRA_IMAGE']:
                mean_erra_image = 'nan'
                median_erra_image = 'nan'
            else:
                mean_erra_image = mean(source_d['ERRA_IMAGE'])
                median_erra_image = median(source_d['ERRA_IMAGE'])

            if 'nan' in source_d['ERRB_IMAGE']:
                mean_errb_image = 'nan'
                median_errb_image = 'nan'
            else:
                mean_errb_image = mean(source_d['ERRB_IMAGE'])
                median_errb_image = median(source_d['ERRB_IMAGE'])

            if 'nan' in source_d['CLASS_STAR']:
                mean_class_star = 'nan'
                median_class_star = 'nan'
            else:
                mean_class_star = mean(source_d['CLASS_STAR'])
                median_class_star = median(source_d['CLASS_STAR'])

            if 'nan' in source_d['FLUX_ISO']:
                mean_flux_iso = 'nan'
                median_flux_iso = 'nan'
            else:
                mean_flux_iso = mean(source_d['FLUX_ISO'])
                median_flux_iso = median(source_d['FLUX_ISO'])

            if 'nan' in source_d['FLUXERR_ISO']:
                mean_fluxerr_iso = 'nan'
                median_fluxerr_iso = 'nan'
            else:
                mean_fluxerr_iso = mean(source_d['FLUXERR_ISO'])
                median_fluxerr_iso = median(source_d['FLUXERR_ISO'])

            if 'nan' in source_d['MAG_ISO']:
                mean_mag_iso = 'nan'
                median_mag_iso = 'nan'
            else:
                mean_mag_iso = mean(source_d['MAG_ISO'])
                median_mag_iso = median(source_d['MAG_ISO'])

            if 'nan' in source_d['MAGERR_ISO']:
                mean_magerr_iso = 'nan'
                median_magerr_iso = 'nan'
            else:
                mean_magerr_iso = mean(source_d['MAGERR_ISO'])
                median_magerr_iso = median(source_d['MAGERR_ISO'])

            if 'nan' in source_d['ELLIPTICITY']:
                mean_ellipticity = 'nan'
                median_ellipticiy = 'nan'
            else:
                mean_ellipticity = mean(source_d['ELLIPTICITY'])
                median_ellipticiy = median(source_d['ELLIPTICITY'])

            # Saves mean and median values from sources
            for i_stats in range(0, len(o_df['SOURCE_NUMBER']), 1):
                tmp_d['MEDIAN_A_IMAGE'].append(median_a_image)
                tmp_d['MEDIAN_B_IMAGE'].append(median_b_image)
                tmp_d['MEDIAN_ERRA_IMAGE'].append(median_erra_image)
                tmp_d['MEDIAN_ERRB_IMAGE'].append(median_errb_image)
                tmp_d['MEDIAN_CLASS_STAR'].append(median_class_star)
                tmp_d['MEDIAN_FLUX_ISO'].append(median_flux_iso)
                tmp_d['MEDIAN_FLUXERR_ISO'].append(median_fluxerr_iso)
                tmp_d['MEDIAN_MAG_ISO'].append(median_mag_iso)
                tmp_d['MEDIAN_MAGERR_ISO'].append(median_magerr_iso)
                tmp_d['MEDIAN_ELLIPTICITY'].append(median_ellipticiy)

                tmp_d['MEAN_A_IMAGE'].append(mean_a_image)
                tmp_d['MEAN_B_IMAGE'].append(mean_b_image)
                tmp_d['MEAN_ERRA_IMAGE'].append(mean_erra_image)
                tmp_d['MEAN_ERRB_IMAGE'].append(mean_errb_image)
                tmp_d['MEAN_CLASS_STAR'].append(mean_class_star)
                tmp_d['MEAN_FLUX_ISO'].append(mean_flux_iso)
                tmp_d['MEAN_FLUXERR_ISO'].append(mean_fluxerr_iso)
                tmp_d['MEAN_MAG_ISO'].append(mean_mag_iso)
                tmp_d['MEAN_MAGERR_ISO'].append(mean_magerr_iso)
                tmp_d['MEAN_ELLIPTICITY'].append(mean_ellipticity)

        series_l = []
        series_d = {}
        for key_ in tmp_d.keys():
            series_d[key_] = Series(tmp_d[key_], name=key_)
            series_l.append(series_d[key_])

        full_db = concat(series_l, axis=1)

        if self.save:
            self.save_message('3_{}'.format(idx_l))
            full_db.to_csv('{}_3_{}.csv'.format(self.filter_o_n, idx_l),
                           columns=dict_keys)

    def filter_class(self, full_db):
        """

        :param full_db:
        :return:
        """
        self.logger.debug('Runs class star filter')
        slow_db = full_db[full_db['MEAN_CLASS_STAR'] > self.class_star_limit]
        fast_db = full_db[full_db['MEAN_CLASS_STAR'] < self.class_star_limit]

        return slow_db, fast_db

    def filter_b(self, full_db):
        """

        :return: full_db
        """
        self.logger.debug('Runs B_Image size filter')

        # Gets unique sources from filtered file
        unique_sources = list(set(full_db['SOURCE_NUMBER'].tolist()))
        l_sourcs = len(unique_sources)  # Just to not break 79 characters
        self.logger.debug('Unique sources to be analysed {}'.format(l_sourcs))

        dict_keys = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'EXTENSION',
                     'ASTR_INSTRUM', 'PHOT_INSTRUM', 'X_IMAGE', 'Y_IMAGE',
                     'ISOAREA_IMAGE', 'A_IMAGE', 'MEDIAN_A_IMAGE',
                     'MEAN_A_IMAGE', 'ERRA_IMAGE',  'MEDIAN_ERRA_IMAGE',
                     'MEAN_ERRA_IMAGE', 'B_IMAGE', 'MEDIAN_B_IMAGE',
                     'MEAN_B_IMAGE', 'ERRB_IMAGE', 'MEDIAN_ERRB_IMAGE',
                     'MEAN_ERRB_IMAGE', 'THETA_IMAGE', 'ERRTHETA_IMAGE',
                     'ALPHA_J2000', 'DELTA_J2000', 'ERRA_WORLD', 'ERRB_WORLD',
                     'ERRTHETA_WORLD', 'EPOCH', 'FWHM_IMAGE', 'CLASS_STAR',
                     'MEDIAN_CLASS_STAR', 'MEAN_CLASS_STAR', 'FLUX_ISO',
                     'MEDIAN_FLUX_ISO', 'MEAN_FLUX_ISO', 'FLUXERR_ISO',
                     'MEDIAN_FLUXERR_ISO', 'MEAN_FLUXERR_ISO', 'FLUX_RADIUS',
                     'ELONGATION', 'ELLIPTICITY', 'MEDIAN_ELLIPTICITY',
                     'MEAN_ELLIPTICITY', 'MAG', 'MAGERR', 'MAG_ISO',
                     'MEDIAN_MAG_ISO', 'MEAN_MAG_ISO', 'MAGERR_ISO',
                     'MEDIAN_MAGERR_ISO', 'MEAN_MAGERR_ISO',
                     'FLAGS_EXTRACTION', 'FLAGS_SCAMP', 'FLAGS_IMA', 'PM',
                     'PMERR', 'PMALPHA', 'PMDELTA', 'PMALPHAERR', 'PMDELTAERR']



        # 0.5, hasta 0.3
        # filter_params = {'central_fit': [-0.14273233, 4.56113342],
        #                  'lower_fit': [-0.14288414, 4.49851371],
        #                  'upper_fit': [-0.14258053, 4.62375314]}

        # 0.4 hasta 0.1
        # filter_params = {'central_fit': [-0.14417005, 4.59286226],
        #                  'lower_fit': [-0.14411595, 4.53657401],
        #                  'upper_fit': [-0.14422416, 4.64915051]}

        # 0.4 hasta 0.3 mags = [21-25]
        # filter_params = {'central_fit': [-0.13489287, 4.39921033],
        #                  'lower_fit': [-0.13488375, 4.31335665],
        #                  'upper_fit': [-0.134902, 4.485064]}

        # 0.4 hasta 1 mags = [21-25]
        # filter_params = {'central_fit': [-0.1352622, 4.41103641],
        #                  'lower_fit': [-0.13525472, 4.32794007],
        #                  'upper_fit': [-0.13526968, 4.49413276]}

        # 0.4 hasta 1 mags = [21-27] not lineal response
        # a-b relation without error
        # original sextractor configuration
        # filter_params = {'central_b_fit': [-0.13816623, 4.47217154],
        #                  'lower_b_fit': [-0.1381669, 4.38973617],
        #                  'upper_b_fit': [-0.13816556, 4.55460691],
        #                  'central_a_fit': [0.93667408, 0.2035752],
        #                  'lower_a_fit': [0.93665283, 0.01074982],
        #                  'upper_a_fit': [0.93669534, 0.39640057]}

        # 0.4 hasta 1 mags = [21-23] lineal response
        # a-b relation without error
        # new sextractor configuration
        filter_params = {'central_b_fit': [-0.16207364, 4.96576972],
                         'lower_b_fit': [-0.16145345, 4.70865816],
                         'upper_b_fit': [-0.16269382, 5.22288128],
                         'central_a_fit': [1.01032694, 0.09742901],
                         'lower_a_fit': [1.00744122, -0.09041792],
                         'upper_a_fit': [1.01321265, 0.28527593]}

        sub_list_1_size = len(unique_sources) / 2
        sub_list_1 = unique_sources[:sub_list_1_size]
        sub_list_2 = unique_sources[sub_list_1_size:]
        sub_list_l = [sub_list_1, sub_list_2]

        areas_j = []
        for idx_l in range(0, 2, 1):
            areas_p = Process(target=self.filter_b_thread,
                              args=(dict_keys, sub_list_l[idx_l], full_db,
                                    filter_params, idx_l,))
            areas_j.append(areas_p)
            areas_p.start()

        active_areas = list([job.is_alive() for job in areas_j])
        while True in active_areas:
            active_areas = list([job.is_alive() for job in areas_j])
            pass

        # Merges areas
        # Merges catalogs
        list_1 = read_csv('{}_6s_0.csv'.format(self.filter_o_n), index_col=0)
        list_2 = read_csv('{}_6s_1.csv'.format(self.filter_o_n), index_col=0)

        full_db = concat([list_1, list_2])

        return full_db

    def filter_b_thread(self, dict_keys, unique_sources_thread, full_db,
                        filter_params, idx_l):
        """

        :param dict_keys:
        :param unique_sources_thread:
        :param full_db:
        :param filter_params:
        :param idx_l:
        :return:
        """
        accepted = []
        rejected = []

        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(unique_sources_thread):
            print('filter_b - thread {} - source {}'.format(idx_l, idx))

            o_df = full_db[full_db['SOURCE_NUMBER'].isin([source_])].iloc[0]
            mag = float(o_df['MEDIAN_MAG_ISO'])

            # b test
            b_image = float(o_df['MEDIAN_B_IMAGE'])
            upper_b_test = (filter_params['upper_b_fit'][0] * mag) + \
                           filter_params['upper_b_fit'][1]
            lower_b_test = (filter_params['lower_b_fit'][0] * mag) + \
                           filter_params['lower_b_fit'][1]
            b_test = float(lower_b_test) < b_image < float(upper_b_test)

            # a test
            a_image = float(o_df['MEDIAN_A_IMAGE'])
            upper_a_test = (filter_params['upper_a_fit'][0] * b_image) + \
                           filter_params['upper_a_fit'][1]
            lower_a_test = (filter_params['lower_a_fit'][0] * b_image) + \
                           filter_params['lower_a_fit'][1]
            a_test = float(lower_a_test) < a_image < float(upper_a_test)

            if b_test and a_test:
                accepted.append(source_)
            else:
                rejected.append(source_)

        full_db = full_db[full_db['SOURCE_NUMBER'].isin(accepted)]

        if self.save:
            self.save_message('6s_{}'.format(idx_l))
            full_db.to_csv('{}_6s_{}.csv'.format(self.filter_o_n, idx_l),
                           columns=dict_keys)

    def filter_coherence(self, full_db):
        """

        :param full_db:
        :return: full_db
        """
        self.logger.debug('Runs coherence motion filter')
        full_db = confidence_filter(full_db, 0.97)

        return full_db

    def filter_detections(self, full_db, detections):
        """

        :param full_db:
        :param detections:
        :return:
        """
        self.logger.debug('Filter by detections number')
        full_db = concat(g for _, g in full_db.groupby("SOURCE_NUMBER")
                         if len(g) >= int(detections))

        return full_db
