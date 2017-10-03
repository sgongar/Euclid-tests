#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

"""

from os import path, makedirs
from subprocess import Popen

from astropy.io import fits
from astropy.table import Table

from cats_management import merge_catalog
from misc import pm_compute, pm_filter, extract_settings
from misc import motion_filter, confidence_filter, get_fits

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


class Scamp:

    def __init__(self, logger, mag, scmp_d, f_conf, sex_d):
        """
        for now scamp's mode is hardcoded

        """
        mode = {'type': 'scamp'}
        prfs_d = extract_settings()

        self.scamp_process(logger, prfs_d, mode, mag,
                           scmp_d, f_conf, sex_d)

    def scamp_process(self, logger, prfs_d, mode, mag,
                      scmp_d, f_conf, sex_d):
        """

        @param logger: a logger object
        @param prfs_d:
        @param mode:
        @param mag:
        @param scmp_d:

        @return True: if everything goes alright
        """

        folder_sex = '{}_{}_{}_{}_{}'.format(sex_d['deblend_nthresh'],
                                             sex_d['analysis_thresh'],
                                             sex_d['detect_thresh'],
                                             sex_d['deblend_mincount'],
                                             sex_d['detect_minarea'])

        f_conf = '{}_{}_{}_{}'.format(scmp_d['crossid_radius'],
                                      scmp_d['pixscale_maxerr'],
                                      scmp_d['posangle_maxerr'],
                                      scmp_d['position_maxerr'])

        logger.info('scamp process for magnitude {}'.format(mag))

        if mode['type'] == 'scamp':
            scmp_p_1 = "scamp -c %s" % (prfs_d['conf_scamp'])
            scmp_p_2 = ' {}/{}/mag_{}_CCD_x?_y?_d?.cat'.format(prfs_d['fits_dir'],
                                                               folder_sex, mag)
            scmp_p_3 = ' -ASTREFCAT_NAME'
            scmp_p_4 = ' {}/{}/catalog_{}.cat'.format(prfs_d['output_cats'],
                                                      folder_sex, mag)
            scmp_p_5 = ' -PIXSCALE_MAXERR {}'.format(scmp_d['pixscale_maxerr'])
            scmp_p_6 = ' -POSANGLE_MAXERR {}'.format(scmp_d['posangle_maxerr'])
            scmp_p_7 = ' -POSITION_MAXERR {}'.format(scmp_d['position_maxerr'])
            scmp_p_8 = ' -CROSSID_RADIUS {}'.format(scmp_d['crossid_radius'])
            cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'], folder_sex,
                                         f_conf)
            merged_cat = '{}/merged_{}_{}_.cat'.format(cats_dir, f_conf, mag)
            scmp_p_9 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat)
            full_cat = '{}/full_{}_{}.cat'.format(cats_dir, f_conf, mag)
            scmp_p_10 = ' -FULLOUTCAT_NAME {}'.format(full_cat)
            scmp_p = scmp_p_1 + scmp_p_2 + scmp_p_3 + scmp_p_4 + scmp_p_5
            scmp_p = scmp_p + scmp_p_6 + scmp_p_7 + scmp_p_8 + scmp_p_9
            scmp_p = scmp_p + scmp_p_10

            # Creates output dir for desired configuration
            output_dir = '{}/catalogs/{}/{}'.format(prfs_d['results_dir'],
                                                    folder_sex, f_conf)

            print "output_dir", output_dir

            if not path.exists(output_dir):
                makedirs(output_dir)

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
            cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'], folder_sex,
                                         f_conf)

            fits_files = get_fits(unique=False)
            for idx, fits_ in enumerate(fits_files):
                merged_cat = '{}/m_{}.cat'.format(cats_dir, fits_[2:-5])
                scmp_p_8 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat)

                full_cat = '{}/f_{}.cat'.format(cats_dir, fits_[2:-5])
                scmp_p_9 = ' -FULLOUTCAT_NAME {}'.format(full_cat)

                # gets fits name
                fits_name = ' {}/{}/{}.cat'.format(prfs_d['fits_dir'],
                                                   folder_sex, fits_[:-5])
                scmp_p_2 = ' {}'.format(fits_name)
                # gets cat name
                reference_cat = ' {}/{}/{}1.cat'.format(prfs_d['fits_ref'],
                                                        folder_sex, fits_[:-6])
                scmp_p_3 = ' -ASTREFCAT_NAME {}'.format(reference_cat)

                scmp_p = scmp_p_1 + scmp_p_2 + scmp_p_3 + scmp_p_4 + scmp_p_5
                scmp_p = scmp_p + scmp_p_6 + scmp_p_7 + scmp_p_8 + scmp_p_9

                output_dir = '{}/catalogs/{}/{}'.format(prfs_d['results_dir'],
                                                        folder_sex, f_conf)

                if not path.exists(output_dir):
                    makedirs(output_dir)

                process_scamp = Popen(scmp_p, shell=True)
                process_scamp.wait()

        logger.info('Scamp process finished. Data is ready to be analysed.')

        return True


class ScampFilter:  # TODO Split scamp_filter method into single methodss

    def __init__(self, logger, mag, scmp_d, f_conf):
        """

        @param mag:
        @param scmp_d:
        @param f_conf:
        """
        prfs_d = extract_settings()

        self.save = True  #"scmp_p",  save flag - set as True for catalogs saving
        self.scamp_filter(logger, prfs_d, mag, scmp_d, f_conf)

    def scamp_filter(self, logger, prfs_d, mag, scmp_d, f_conf):
        """

        @param logger:
        @param prfs_d:
        @param mag:
        @param scmp_d:
        """
        logger.info("Filtering scamp's output")

        # Full catalog name
        full_n = '{}/full_{}_{}_1.cat'.format(prfs_d['results_dir'],
                                              f_conf, mag)
        # Filtered catalog name
        filt_n = 'filt_{}_{}_'.format(f_conf, mag)

        logger.debug('opening full catalog {}'.format(full_n))
        full_cat = fits.open(full_n)
        full_db = Table(full_cat[2].data)
        logger.debug('converting full catalog to Pandas format')
        full_db = full_db.to_pandas()

        # Getting merge catalogue
        mrgd_n = '{}/merged_{}_{}_1.cat'.format(prfs_d['results_dir'],
                                                f_conf, mag)
        logger.debug('opening merged catalog {}'.format(mrgd_n))
        merged_cat = fits.open(mrgd_n)
        logger.debug('converting merged catalog to Pandas format')
        merged_db = Table(merged_cat[2].data)

        # Removing 0 catalogue detections
        logger.debug('removing 0 catalog detections')
        full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]
        if self.save:
            full_db.to_csv('{}/{}_1.csv'.format(prfs_d['results_dir'], filt_n))

        # Computing pm
        logger.debug('computing proper motion')
        full_db = pm_compute(logger, merged_db, full_db)
        if self.save:
            full_db.to_csv('{}/{}_3.csv'.format(prfs_d['results_dir'], filt_n))

        logger.debug('after filtering detections')
        full_db = pm_filter(full_db, prfs_d['pm_low'],
                            prfs_d['pm_up'], prfs_d['pm_sn'])
        if self.save:
            full_db.to_csv('{}/{}_4.csv'.format(prfs_d['results_dir'], filt_n))

        logger.debug('after proper motion')
        full_db = motion_filter(logger, full_db, prfs_d['r_fit'])
        if self.save:
            full_db.to_csv('{}/{}_5.csv'.format(prfs_d['results_dir'], filt_n))

        logger.debug('after first filter')
        full_db = confidence_filter(logger, full_db, prfs_d['r_fit'])
        if self.save:
            full_db.to_csv('{}/{}_6.csv'.format(prfs_d['results_dir'], filt_n))


class CatalogCreation:

    def __init__(self):
        """

        """
        prfs_d = extract_settings()
        folder_sex = '2_1.35_1.35_0.1_4'

        if not self.catalog_creation(prfs_d, folder_sex):
            raise Exception

    def catalog_creation(self, prfs_d, folder_sex):
        """

        """
        merge_catalog(prfs_d, folder_sex)

        return True


if __name__ == '__main__':
    test = CatalogCreation()
