#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.
- 0.2: Functions reorganised to two different classes.
- 0.3: New organization for filter Class

"""

from os import path, makedirs
from subprocess import Popen

from astropy.io import fits
from astropy.table import Table
from pandas import concat

from misc import pm_compute, pm_filter, extract_settings
from misc import sn_filter, create_folder
from misc import motion_filter, confidence_filter

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.3"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class Scamp:

    def __init__(self, logger, mag, scmp_d, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_d:
        :param scmp_cf:
        :param sex_d:
        """
        self.prfs_d = extract_settings()

        self.scamp_process(logger, mag, scmp_d, scmp_cf, sex_d)

    def scamp_process(self, logger, mag, scmp_d, scmp_cf, sex_d):
        """

        @param logger: a logger object
        @param mag:
        @param scmp_d:

        @return True: if everything goes alright
        """
        sex_cf = '{}_{}_{}_{}_{}'.format(sex_d['deblend_nthresh'],
                                         sex_d['analysis_thresh'],
                                         sex_d['detect_thresh'],
                                         sex_d['deblend_mincount'],
                                         sex_d['detect_minarea'])

        logger.info('scamp process for magnitude {}'.format(mag))
        logger.info('sextractor configuration {}'.format(sex_cf))
        logger.info('scamp configuration {}'.format(scmp_cf))

        scmp_p_1 = 'scamp -c {}'.format(self.prfs_d['conf_scamp'])
        sex_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], mag,
                                         sex_cf)
        sex_output = 'mag_{}_CCD_x?_y?_d?.cat'.format(mag)
        scmp_p_2 = ' {}/{}'.format(sex_loc, sex_output)
        scmp_p_3 = ' -ASTREFCAT_NAME'
        cat_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], mag,
                                         sex_cf)
        cat_input = 'catalog_{}.cat'.format(mag)
        scmp_p_4 = ' {}/{}'.format(cat_loc, cat_input)
        scmp_p_5 = ' -PIXSCALE_MAXERR {}'.format(scmp_d['pixscale_maxerr'])
        scmp_p_6 = ' -POSANGLE_MAXERR {}'.format(scmp_d['posangle_maxerr'])
        scmp_p_7 = ' -POSITION_MAXERR {}'.format(scmp_d['position_maxerr'])
        scmp_p_8 = ' -CROSSID_RADIUS {}'.format(scmp_d['crossid_radius'])
        # Output catalogs location
        cats_dir = '{}/{}/{}/{}'.format(self.prfs_d['catalogs_dir'], mag,
                                        sex_cf, scmp_cf)
        merged_cat = '{}/merged_{}_{}.cat'.format(cats_dir, scmp_cf, mag)
        scmp_p_9 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat)
        full_cat = '{}/full_{}_{}.cat'.format(cats_dir, scmp_cf, mag)
        scmp_p_10 = ' -FULLOUTCAT_NAME {}'.format(full_cat)
        scmp_p = scmp_p_1 + scmp_p_2 + scmp_p_3 + scmp_p_4 + scmp_p_5
        scmp_p = scmp_p + scmp_p_6 + scmp_p_7 + scmp_p_8 + scmp_p_9
        scmp_p = scmp_p + scmp_p_10

        # print('scmp_p {}'.format(scmp_p))

        create_folder(logger, cats_dir)

        process_scamp = Popen(scmp_p, shell=True)
        process_scamp.wait()

        logger.info('Scamp process finished. Data is ready to be analysed.')

        return True


class ScampFilter:  # TODO Split scamp_filter method into single methodss

    def __init__(self, logger, mag, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_d:
        """
        self.prfs_d = extract_settings()

        self.save = True
        (merged_db, full_db,
         filter_o_n) = self.scamp_filter(logger, mag, scmp_cf, sex_d)
        full_db = self.compute_pm(logger, merged_db, full_db, filter_o_n)
        self.filter_pm(logger, full_db, filter_o_n)

    def scamp_filter(self, logger, mag, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_d:
        :return:
        """
        logger.info("Filtering scamp's output")

        sex_cf = '{}_{}_{}_{}_{}'.format(sex_d['deblend_nthresh'],
                                         sex_d['analysis_thresh'],
                                         sex_d['detect_thresh'],
                                         sex_d['deblend_mincount'],
                                         sex_d['detect_minarea'])

        filter_dir = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'], mag,
                                          sex_cf, scmp_cf)
        create_folder(logger, filter_dir)

        logger.debug('running filter for scamp {} and {} sex'.format(scmp_cf,
                                                                     sex_cf))

        # Full catalog name
        full_n = '{}/{}/{}/{}/full_{}_{}_1.cat'.format(self.prfs_d['catalogs_dir'],
                                                       mag, sex_cf, scmp_cf,
                                                       scmp_cf, mag)
        # Filtered catalog name
        filt_n = 'filt_{}_{}'.format(scmp_cf, mag)

        logger.debug('opening full catalog {}'.format(full_n))
        full_cat = fits.open(full_n)
        full_db = Table(full_cat[2].data)
        logger.debug('converting full catalog to Pandas format')
        full_db = full_db.to_pandas()

        # Getting merge catalog
        mrgd_n = '{}/{}/{}/{}/merged_{}_{}_1.cat'.format(self.prfs_d['catalogs_dir'],
                                                         mag, sex_cf, scmp_cf,
                                                         scmp_cf, mag)

        logger.debug('opening merged catalog {}'.format(mrgd_n))
        merged_cat = fits.open(mrgd_n)
        logger.debug('converting merged catalog to Pandas format')
        merged_db = Table(merged_cat[2].data)

        filter_o_n = '{}/{}'.format(filter_dir, filt_n)

        # Removing 0 catalog detections
        logger.debug('removing 0 catalog detections')
        full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

        full_db = concat(g for _, g in full_db.groupby("SOURCE_NUMBER")
                         if len(g) >= int(self.prfs_d['detections']))

        if self.save:
            full_db.to_csv('{}_1.csv'.format(filter_o_n))

        return merged_db, full_db, filter_o_n

    def compute_pm(self, logger, merged_db, full_db, filter_o_n):
        """

        :param logger:
        :param merged_db:
        :param full_db:
        :param filter_o_n:
        :return:
        """
        # Computing pm
        logger.debug('computing proper motion')
        full_db = pm_compute(logger, merged_db, full_db)
        if self.save:
            logger.debug('saving output to {}_2.csv'.format(filter_o_n))
            full_db.to_csv('{}_2.csv'.format(filter_o_n))

        return full_db

    def filter_pm(self, logger, full_db, filter_o_n):
        """

        :param logger:
        :param full_db:
        :param filter_o_n:
        :return:
        """
        logger.debug('after filtering detections')
        full_db = pm_filter(full_db, self.prfs_d['pm_low'],
                            self.prfs_d['pm_up'])
        if self.save:
            logger.debug('saving output to {}_3.csv'.format(filter_o_n))
            full_db.to_csv('{}_3.csv'.format(filter_o_n))

        full_db = sn_filter(full_db, self.prfs_d['pm_sn'])
        if self.save:
            logger.debug('saving output to {}_4.csv'.format(filter_o_n))
            full_db.to_csv('{}_4.csv'.format(filter_o_n))

        logger.debug('after proper motion')
        full_db = motion_filter(full_db, self.prfs_d['r_fit'])
        if self.save:
            logger.debug('saving output to {}_5.csv'.format(filter_o_n))
            full_db.to_csv('{}_5.csv'.format(filter_o_n))

        logger.debug('after first filter')
        full_db = confidence_filter(full_db, self.prfs_d['r_fit'])
        if self.save:
            logger.debug('saving output to {}_6.csv'.format(filter_o_n))
            full_db.to_csv('{}_6.csv'.format(filter_o_n))
