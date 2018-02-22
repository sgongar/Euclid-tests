#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.
- 0.2: Pipeline splitted.

Todo:
    * Improve documentation
    * Speed-up areas process
"""

from multiprocessing import Process

from astropy.io import fits
from astropy.table import Table
from pandas import concat, read_csv, Series

from misc import pm_compute, extract_settings, check_source
from misc import b_filter, create_folder, get_dither, confidence_filter

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class ScampFilter:  # TODO Split scamp_filter method into single methods

    def __init__(self, logger, mag, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_d:
        """
        # Filter variables
        self.class_star_limit = 0.97
        self.proper_motion = 1
        self.proper_motion_dects = 0.1

        # Analysis variables
        self.prfs_d = extract_settings()
        self.logger = logger
        self.mag = mag
        self.scmp_cf = scmp_cf
        self.sex_cf = '{}_{}_{}_{}_{}'.format(sex_d['deblend_nthresh'],
                                              sex_d['analysis_thresh'],
                                              sex_d['detect_thresh'],
                                              sex_d['deblend_mincount'],
                                              sex_d['detect_minarea'])

        self.save = True

        # Filtered catalog dir
        self.filter_dir = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                               self.mag, self.sex_cf,
                                               self.scmp_cf)
        # Filtered catalog name
        self.filt_n = 'filt_{}_{}'.format(self.scmp_cf, self.mag)

        self.filter_o_n = '{}/{}'.format(self.filter_dir, self.filt_n)
        """
        # Saves _1.csv
        (merged_db, full_db) = self.scamp_filter()
        # Saves _2.csv
        full_db = self.compute_pm(merged_db, full_db)
        self.get_areas()
        """
        full_db = read_csv('{}_3.csv'.format(self.filter_o_n))
        self.slow_db, self.fast_db = self.filter_class(full_db)
        if self.save:
            self.save_message('4s')
            self.slow_db.to_csv('{}_4s.csv'.format(self.filter_o_n))
            self.save_message('4f')
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
        slow_db = self.filter_detections(slow_db, 3)
        if self.save:
            self.save_message('7s')
            slow_db.to_csv('{}_7s.csv'.format(self.filter_o_n))

        fast_db = full_db[full_db['PM'] > self.proper_motion_dects]
        fast_db = self.filter_detections(fast_db, 1)
        if self.save:
            self.save_message('7f')
            fast_db.to_csv('{}_7f.csv'.format(self.filter_o_n))

        full_db = concat([slow_db, fast_db])

        if self.save:
            self.save_message('8')
            full_db.to_csv('{}_8.csv'.format(self.filter_o_n))

        full_db = full_db[full_db['PM'] > 0.0001]

        if self.save:
            self.save_message('9')
            full_db.to_csv('{}_9.csv'.format(self.filter_o_n))

    def save_message(self, order):
        """

        :param order:
        :return:
        """
        self.logger.debug('Saves data to: ')
        # self.logger.debug('Dir: {}'.format(self.filter_dir))
        self.logger.debug('Name: {}_{}.csv'.format(self.filt_n, order))

    def get_cat(self, catalog_n):
        """

        :param catalog_n:
        :return: cat_file
        """
        ccd, dither = get_dither(catalog_n)

        cat_loc = '{}/{}/CCDs/{}/'.format(self.prfs_d['fits_dir'], self.mag,
                                          self.sex_cf)
        cat_name = 'mag_{}_CCD_{}_d{}.cat'.format(self.mag, ccd, dither)
        cat_file = '{}{}'.format(cat_loc, cat_name)

        return cat_file

    def scamp_filter(self):
        """

        :return:
        """
        self.logger.info("Filtering scamp's output")

        filter_dir = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'], self.mag,
                                          self.sex_cf, self.scmp_cf)
        create_folder(self.logger, filter_dir)

        self.logger.debug('Scamp configuration: {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration: {}'.format(self.sex_cf))

        # Full catalog name
        full_n_dir = '{}/{}/{}/{}/'.format(self.prfs_d['catalogs_dir'],
                                           self.mag, self.sex_cf, self.scmp_cf)
        full_n_cat = 'full_{}_{}_1.cat'.format(self.scmp_cf, self.mag)
        full_n = '{}{}'.format(full_n_dir, full_n_cat)
        # Filtered catalog name

        self.logger.debug('Opens full catalog')
        self.logger.debug('Dir: {}'.format(full_n_dir))
        self.logger.debug('Name: {}'.format(full_n_cat))
        full_cat = fits.open(full_n)
        full_db = Table(full_cat[2].data)
        self.logger.debug('Converts full Astropy catalog to Pandas format')
        full_db = full_db.to_pandas()

        # Getting merge catalog
        mrgd_n_dir = '{}/{}/{}/{}/'.format(self.prfs_d['catalogs_dir'],
                                           self.mag, self.sex_cf, self.scmp_cf)
        mrgd_n_cat = 'merged_{}_{}_1.cat'.format(self.scmp_cf, self.mag)
        mrgd_n = '{}{}'.format(mrgd_n_dir, mrgd_n_cat)

        self.logger.debug('Opens merged catalog')
        self.logger.debug('Dir: {}'.format(mrgd_n_dir))
        self.logger.debug('Name: {}'.format(mrgd_n_cat))
        merged_cat = fits.open(mrgd_n)
        self.logger.debug('Converts merged Astropy catalog to Pandas format')
        merged_db = Table(merged_cat[2].data)

        # Removing 0 catalog detections
        self.logger.debug('Removes 0 catalog detections')
        full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

        full_db = self.filter_detections(full_db)

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
        self.logger.debug('computing proper motion')
        full_db = pm_compute(self.logger, merged_db, full_db)
        if self.save:
            self.save_message('2')
            full_db.to_csv('{}_2.csv'.format(self.filter_o_n))

        return full_db

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
        self.logger.debug('runs areas filter')

        # Opens filtered file
        filter_cat = read_csv('{}_2.csv'.format(self.filter_o_n), index_col=0)
        # Gets unique sources from filtered file
        unique_sources = list(set(filter_cat['SOURCE_NUMBER'].tolist()))

        extra_keys = ['A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'ISOAREA_IMAGE',
                      'FWHM_IMAGE', 'FLUX_RADIUS', 'MAG_ISO', 'ELONGATION',
                      'ELLIPTICITY', 'CLASS_STAR']
        keys_l = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'X_IMAGE', 'Y_IMAGE',
                  'ERRA_IMAGE', 'ERRB_IMAGE', 'ERRTHETA_IMAGE', 'ALPHA_J2000',
                  'DELTA_J2000', 'ERRA_WORLD', 'ERRB_WORLD', 'ERRTHETA_WORLD',
                  'EPOCH', 'MAG', 'MAGERR', 'FLAGS_EXTRACTION', 'FLAGS_SCAMP',
                  'FLAGS_IMA', 'PM', 'PMERR', 'PMALPHA', 'PMDELTA',
                  'PMALPHAERR', 'PMDELTAERR']

        sub_list_1_size = len(unique_sources) / 2
        sub_list_1 = unique_sources[:sub_list_1_size]
        sub_list_2 = unique_sources[sub_list_1_size:]
        sub_list_l = [sub_list_1, sub_list_2]

        areas_j = []
        for idx_l in range(0, 2, 1):
            areas_p = Process(target=self.get_areas_thread,
                              args=(extra_keys, keys_l, sub_list_l[idx_l],
                                    filter_cat, idx_l,))
            areas_j.append(areas_p)
            areas_p.start()

        active_areas = list([job.is_alive() for job in areas_j])
        while True in active_areas:
            active_areas = list([job.is_alive() for job in areas_j])
            pass

        # Merges areas
        # Merges catalogs
        list_1 = read_csv('{}_3_1.csv'.format(self.filter_o_n), index_col=0)
        list_2 = read_csv('{}_3_2.csv'.format(self.filter_o_n), index_col=0)

        full_db = concat([list_1, list_2])

        if self.save:
            self.save_message('3')
            full_db.to_csv('{}_3.csv'.format(self.filter_o_n))

    def get_areas_thread(self, extra_keys, keys_l, unique_sources_thread,
                         filter_cat, idx_l):
        """

        :param extra_keys:
        :param keys_l:
        :param unique_sources_thread:
        :param filter_cat:
        :param idx:
        :return:
        """
        tmp_d = {}
        for key_ in keys_l:
            tmp_d[key_] = []

        for key_ in extra_keys:
            tmp_d[key_] = []

        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(unique_sources_thread):
            print(idx)
            o_df = filter_cat[filter_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                # Populates temporal dictionary
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000

                catalog_n = tmp_d['CATALOG_NUMBER']
                cat_file = self.get_cat(catalog_n)
                # self.logger.debug('opening CCD catalog {}'.format(cat_file))

                ccd_cat = fits.open(cat_file)
                ccd_df = Table(ccd_cat[2].data)
                # self.logger.debug('converting CCD catalog to Pandas format')
                ccd_df = ccd_df.to_pandas()

                cat_df = check_source(ccd_df, o_alpha, o_delta)
                if cat_df.empty:
                    # FIXME problemas entre sextractor y scamp
                    # Como el objeto no esta en las imagenes originales de
                    # sextractor lo borro del catalogo de scamp
                    for idx_k, key_ in enumerate(keys_l):
                        tmp_d[key_].append('nan')
                    for idx_k, key_ in enumerate(extra_keys):
                        tmp_d[key_].append('nan')
                else:
                    for idx_k, key_ in enumerate(keys_l):
                        tmp_d[key_].append(row[idx_k + 1])
                    for idx_k, key_ in enumerate(extra_keys):
                        tmp_d[key_].append(cat_df[key_].iloc[0])

        series_l = []
        series_d = {}
        for key_ in tmp_d.keys():
            series_d[key_] = Series(tmp_d[key_], name=key_)
            series_l.append(series_d[key_])

        full_db = concat([series_l], axis=1)

        if self.save:
            self.save_message('3_{}'.format(idx_l))
            full_db.to_csv('{}_3_{}.csv'.format(self.filter_o_n, idx_l))

    def filter_class(self, full_db):
        """

        :param full_db:
        :return:
        """
        self.logger.debug('Runs class star filter')
        slow_db = full_db[full_db['CLASS_STAR'] > self.class_star_limit]
        fast_db = full_db[full_db['CLASS_STAR'] < self.class_star_limit]

        return slow_db, fast_db

    def filter_b(self, full_db):
        """
        1.57 - 1.8
        :param full_db:
        :return: full_db
        """
        self.logger.debug('Runs B_Image size filter')
        full_db = b_filter(full_db, 1.57, 1.7)

        return full_db

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
