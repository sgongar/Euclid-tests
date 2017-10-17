#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

Versions:
- 0.1: Initial release.
- 0.2: Functions reorganised to two different classes.

"""

from os import path, makedirs
from subprocess import Popen

from astropy.io import fits
from astropy.table import Table
from pandas import concat

from cats_management import merge_catalog
from misc import pm_compute, pm_filter, extract_settings
from misc import motion_filter, confidence_filter, get_fits

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.2"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class Scamp:

    def __init__(self, logger, mag, scmp_d, scmp_cf, sex_d):
        """
        for now scamp's mode is hardcoded

        """
        mode = {'type': 'scamp'}
        prfs_d = extract_settings()

        self.scamp_process(logger, prfs_d, mode, mag,
                           scmp_d, scmp_cf, sex_d)

    def scamp_process(self, logger, prfs_d, mode, mag,
                      scmp_d, scmp_cf, sex_d):
        """

        @param logger: a logger object
        @param prfs_d:
        @param mode:
        @param mag:
        @param scmp_d:

        @return True: if everything goes alright
        """

        sex_cf = '{}_{}_{}_{}_{}'.format(sex_d['deblend_nthresh'],
                                         sex_d['analysis_thresh'],
                                         sex_d['detect_thresh'],
                                         sex_d['deblend_mincount'],
                                         sex_d['detect_minarea'])

        scmp_cf = '{}_{}_{}_{}'.format(scmp_d['crossid_radius'],
                                       scmp_d['pixscale_maxerr'],
                                       scmp_d['posangle_maxerr'],
                                       scmp_d['position_maxerr'])

        logger.info('scamp process for magnitude {}'.format(mag))
        logger.info('sextractor configuration {}'.format(sex_cf))
        logger.info('scamp configuration {}'.format(scmp_cf))

        if mode['type'] == 'scamp':
            scmp_p_1 = "scamp -c %s" % (prfs_d['conf_scamp'])
            scmp_p_2 = ' {}/{}/mag_{}_CCD_x?_y?_d?.cat'.format(prfs_d['fits_dir'],
                                                               sex_cf, mag)
            scmp_p_3 = ' -ASTREFCAT_NAME'
            scmp_p_4 = ' {}/{}/catalog_{}.cat'.format(prfs_d['output_cats'],
                                                      sex_cf, mag)
            scmp_p_5 = ' -PIXSCALE_MAXERR {}'.format(scmp_d['pixscale_maxerr'])
            scmp_p_6 = ' -POSANGLE_MAXERR {}'.format(scmp_d['posangle_maxerr'])
            scmp_p_7 = ' -POSITION_MAXERR {}'.format(scmp_d['position_maxerr'])
            scmp_p_8 = ' -CROSSID_RADIUS {}'.format(scmp_d['crossid_radius'])
            cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'], sex_cf,
                                         scmp_cf)
            merged_cat = '{}/merged_{}_{}.cat'.format(cats_dir, scmp_cf, mag)
            scmp_p_9 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat)
            full_cat = '{}/full_{}_{}.cat'.format(cats_dir, scmp_cf, mag)
            scmp_p_10 = ' -FULLOUTCAT_NAME {}'.format(full_cat)
            scmp_p = scmp_p_1 + scmp_p_2 + scmp_p_3 + scmp_p_4 + scmp_p_5
            scmp_p = scmp_p + scmp_p_6 + scmp_p_7 + scmp_p_8 + scmp_p_9
            scmp_p = scmp_p + scmp_p_10

            # Creates output dir for desired configuration
            output_dir = '{}/catalogs/{}/{}'.format(prfs_d['results_dir'],
                                                    sex_cf, scmp_cf)

            if not path.exists(output_dir):
                makedirs(output_dir)

            print scmp_p

            process_scamp = Popen(scmp_p, shell=True)
            process_scamp.wait()

        elif mode['type'] == 'custom':

            scmp_p_1 = "scamp -c %s" % (prfs_d['conf_scamp'])
            scmp_p_4 = ' -PIXSCALE_MAXERR {}'.format(scmp_d['pixscale_maxerr'])
            scmp_p_5 = ' -POSANGLE_MAXERR {}'.format(scmp_d['posangle_maxerr'])
            scmp_p_6 = ' -POSITION_MAXERR {}'.format(scmp_d['position_maxerr'])
            scmp_p_7 = ' -CROSSID_RADIUS {}'.format(scmp_d['crossid_radius'])

            # cats_dir value is associated to configuration, not to the file
            # so it should be generated once
            cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'],
                                         sex_cf, scmp_cf)

            fits_files = get_fits(unique=False)
            for idx, fits_ in enumerate(fits_files):
                merged_cat = '{}/m_{}.cat'.format(cats_dir, fits_[2:-5])
                scmp_p_8 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat)

                full_cat = '{}/f_{}.cat'.format(cats_dir, fits_[2:-5])
                scmp_p_9 = ' -FULLOUTCAT_NAME {}'.format(full_cat)

                # gets fits name
                fits_name = ' {}/{}/{}.cat'.format(prfs_d['fits_dir'],
                                                   sex_cf, fits_[:-5])
                scmp_p_2 = ' {}'.format(fits_name)
                # gets cat name
                reference_cat = ' {}/{}/{}1.cat'.format(prfs_d['fits_ref'],
                                                        sex_cf, fits_[:-6])
                scmp_p_3 = ' -ASTREFCAT_NAME {}'.format(reference_cat)

                scmp_p = scmp_p_1 + scmp_p_2 + scmp_p_3 + scmp_p_4 + scmp_p_5
                scmp_p = scmp_p + scmp_p_6 + scmp_p_7 + scmp_p_8 + scmp_p_9

                output_dir = '{}/catalogs/{}/{}'.format(prfs_d['results_dir'],
                                                        sex_cf, scmp_cf)

                if not path.exists(output_dir):
                    makedirs(output_dir)

                process_scamp = Popen(scmp_p, shell=True)
                process_scamp.wait()

        logger.info('Scamp process finished. Data is ready to be analysed.')

        return True


class ScampFilter:  # TODO Split scamp_filter method into single methodss

    def __init__(self, logger, mag, scmp_d, scmp_cf, sex_d):
        """

        @param mag:
        @param scmp_d:
        @param scmp_cf:
        """
        prfs_d = extract_settings()

        self.save = True  # "scmp_p",  save flag - set as True for catalogs saving
        (merged_db, full_db,
         filter_o_n) = self.scamp_filter(logger, prfs_d, mag,
                                         scmp_d, scmp_cf, sex_d)
        full_db = self.compute_pm(logger, prfs_d, merged_db,
                                  full_db, filter_o_n)
        self.filter_pm(logger, prfs_d, full_db, filter_o_n)

    def scamp_filter(self, logger, prfs_d, mag, scmp_d, scmp_cf, sex_d):
        """

        @param logger:
        @param prfs_d:
        @param mag:
        @param scmp_d:
        """
        logger.info("Filtering scamp's output")

        sex_cf = '{}_{}_{}_{}_{}'.format(sex_d['deblend_nthresh'],
                                         sex_d['analysis_thresh'],
                                         sex_d['detect_thresh'],
                                         sex_d['deblend_mincount'],
                                         sex_d['detect_minarea'])

        if not self.check_dir(logger, prfs_d, sex_cf, scmp_cf):
            raise Exception

        logger.debug('running filter for scamp {} and {} sex'.format(scmp_cf,
                                                                     sex_cf))

        # Full catalog name
        full_n = '{}/{}/{}/full_{}_{}_1.cat'.format(prfs_d['catalogs_dir'],
                                                    sex_cf, scmp_cf,
                                                    scmp_cf, mag)
        # Filtered catalog name
        filt_n = 'filt_{}_{}'.format(scmp_cf, mag)

        logger.debug('opening full catalog {}'.format(full_n))
        full_cat = fits.open(full_n)
        full_db = Table(full_cat[2].data)
        logger.debug('converting full catalog to Pandas format')
        full_db = full_db.to_pandas()

        # Getting merge catalog
        mrgd_n = '{}/{}/{}/merged_{}_{}_1.cat'.format(prfs_d['catalogs_dir'],
                                                      sex_cf, scmp_cf,
                                                      scmp_cf, mag)

        logger.debug('opening merged catalog {}'.format(mrgd_n))
        merged_cat = fits.open(mrgd_n)
        logger.debug('converting merged catalog to Pandas format')
        merged_db = Table(merged_cat[2].data)

        filter_o_n = '{}/{}/{}/{}'.format(prfs_d['filter_dir'],
                                          sex_cf, scmp_cf, filt_n)

        # Removing 0 catalog detections
        logger.debug('removing 0 catalog detections')
        full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

        full_db = concat(g for _, g in full_db.groupby("SOURCE_NUMBER")
                         if len(g) >= int(prfs_d['detections']))

        if self.save:
            full_db.to_csv('{}_1.csv'.format(filter_o_n))

        return merged_db, full_db, filter_o_n

    def compute_pm(self, logger, prfs_d, merged_db, full_db, filter_o_n):
        """

        @param logger:
        @param prfs_d:
        @param merged_db:
        @param full_db:
        @param filter_o_n:

        @return full_db:
        """
        # Computing pm
        logger.debug('computing proper motion')
        full_db = pm_compute(logger, merged_db, full_db)
        if self.save:
            logger.debug('saving output to {}_2.csv'.format(filter_o_n))
            full_db.to_csv('{}_2.csv'.format(filter_o_n))

        return full_db

    def filter_pm(self, logger, prfs_d, full_db, filter_o_n):
        """

        @param logger:
        @param prfs_d:
        @param full_db:
        """
        logger.debug('after filtering detections')
        full_db = pm_filter(full_db, prfs_d['pm_low'],
                            prfs_d['pm_up'], prfs_d['pm_sn'])
        if self.save:
            logger.debug('saving output to {}_3.csv'.format(filter_o_n))
            full_db.to_csv('{}_3.csv'.format(filter_o_n))

        logger.debug('after proper motion')
        full_db = motion_filter(logger, full_db, prfs_d['r_fit'])
        if self.save:
            logger.debug('saving output to {}_4.csv'.format(filter_o_n))
            full_db.to_csv('{}_4.csv'.format(filter_o_n))

        logger.debug('after first filter')
        full_db = confidence_filter(logger, full_db, prfs_d['r_fit'])
        if self.save:
            logger.debug('saving output to {}_5.csv'.format(filter_o_n))
            full_db.to_csv('{}_5.csv'.format(filter_o_n))

    def check_dir(self, logger, prfs_d, sex_cf, scmp_cf):
        """

        @param logger:
        @param sex_cf:
        @param scmp_cf:

        @return True: if everything goes alright.
        """
        filter_dir = '{}/{}/{}'.format(prfs_d['filter_dir'], sex_cf, scmp_cf)

        if not path.exists(filter_dir):
            makedirs(filter_dir)
        else:
            logger.debug('directory {} already exists'.format(filter_dir))

        return True
