#!/usr/bin/python
# -*- coding: utf-8 -*-


"""Python script for sextractor routines

Conventions:
    * Single n means name.
    * Single d means dictionary.
    * Single j means jobs.
    * Single p means process.
    * loc means lcoation.

Todo:
    * Improve log messages

"""


from os import mkdir, path, remove
from subprocess import Popen

from multiprocessing import Process

from cats_management import rebase_catalogue
from cats_management import rewriting_catalogue
from misc import extract_settings, get_fits


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


class Sextractor:

    def __init__(self, logger, analysis_d, analysis_dir, regular):
        """

        """
        prfs_d = extract_settings()
        self.sextractor_process(logger, prfs_d,
                                analysis_d, analysis_dir, regular)

    def sextractor_process(self, logger, prfs_d, analysis_d,
                           analysis_dir, regular):
        """ runs sextractor over all ccds or fpas images

        @param logger: a logger object
        @param prfs_d: a dictionary with all configuration information
        @param analysis_d:
        @param analysis_dir:
        @param regular: a boolean variable, fits files regular struct?

        @return True: if everything goes alright
        """
        logger.info('Starting sextractor process for fits images')

        mags = prfs_d['mags']
        fits_files = get_fits(unique=False)

        folder_n = '{}_{}_{}_{}_{}'.format(analysis_d['deblend_nthresh'],
                                           analysis_d['analysis_thresh'],
                                           analysis_d['detect_thresh'],
                                           analysis_d['deblend_mincount'],
                                           analysis_d['detect_minarea'])

        # Creates folder for each configuration
        folder_loc = '{}/{}'.format(analysis_dir, folder_n)
        try:
            if not path.isfile(folder_loc):
                mkdir(folder_loc)
        except OSError:
            logger.debug('Folder {} already created'.format(folder_n))

        for image_idx in range(0, len(fits_files), prfs_d['cores_number']):
            for mag_ in mags:  # TODO Implement mag selection
                try:
                    sex_j = []
                    for proc in range(0, prfs_d['cores_number'], 1):
                        idx = image_idx + proc  # index

                        sex_file = fits_files[idx]
                        sex_output = '{}/{}.cat'.format(folder_n,
                                                        fits_files[idx][:-5])

                        sex_p = Process(target=self.sextractor_thread,
                                        args=(logger, prfs_d, sex_file,
                                              sex_output, analysis_d,
                                              analysis_dir,))
                        sex_j.append(sex_p)
                        sex_p.start()

                    active_sex = list([job.is_alive() for job in sex_j])
                    while True in active_sex:
                        active_sex = list([job.is_alive() for job in sex_j])
                        pass
                except IndexError:
                    logger.debug('Extraction finished')

        logger.info('Extraction process of fits images finished')

        return True

    def sextractor_thread(self, logger, prfs_d, sextractor_file,
                          sextractor_output, analysis_d, analysis_dir):
        """ runs sextractor on a single file

        @param logger: a logger object
        @param prfs_d: a dictionary with all configuration information
        @param sextractor_file: file to be "sextracted"
        @param sextractor_output: catalogue to be created by sextractor
        @param analysis_d: parameters dict

        @return True: if everything goes alright
        """
        s_1 = "sex -c %s %s/%s" % (prfs_d['conf_sex'], analysis_dir,
                                   sextractor_file)
        s_2 = " -CATALOG_NAME %s/%s" % (analysis_dir,
                                        sextractor_output)
        s_3 = " -PARAMETERS_NAME %s" % (prfs_d['params_sex'])
        s_4 = " -DETECT_MINAREA %s" % (analysis_d['detect_minarea'])
        s_5 = " -DETECT_THRESH %s" % (analysis_d['detect_thresh'])
        s_6 = " -ANALYSIS_THRESH %s" % (analysis_d['analysis_thresh'])
        s_7 = " -DEBLEND_NTHRESH %s" % (analysis_d['deblend_nthresh'])
        s_8 = " -DEBLEND_MINCONT %s" % (analysis_d['deblend_mincount'])
        s_9 = " -FILTER_NAME %s" % (analysis_d['filter'])

        cmd_3 = s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9

        sextractor_p = Popen(cmd_3, shell=True)
        sextractor_p.wait()

        return True


class CatalogCreation:

    def __init__(self, logger, analysis_d):
        """ creates a single catalogue for all images available

        :param logger: a logger object
        :param analysis_d: a dictionary with all configuration information
        """
        prfs_d = extract_settings()
        mags = prfs_d['mags']

        cat_j = []
        for proc in range(0, len(mags), 1):
            mag = mags[proc]
            cat_p = Process(target=self.catalogue_thread,
                            args=(logger, prfs_d, analysis_d, mag,))
            cat_j.append(cat_p)
            cat_p.start()

        active_cat = list([job.is_alive() for job in cat_j])
        while True in active_cat:
            active_cat = list([job.is_alive() for job in cat_j])
            pass

    def catalogue_thread(self, logger, prfs_d, analysis_d, mag):
        """

        @param logger: a logger object
        @param prfs_d: a dictionary with all configuration information
        @param analysis_d:
        @param mag:

        @return True:
        """
        # Fits files location
        coadd_loc = '{}/coadd_{}.fits'.format(prfs_d['fits_ref'], mag)
        coadd_weight_loc = '{}/coadd_{}.weight.fits'.format(prfs_d['fits_ref'],
                                                            mag)

        # Harcoded for single CCD catalog!
        logger.info('swarp process launched for mag {}'.format(mag))
        logger.info('files {}/mag_{}_*.fits'.format(prfs_d['fits_ref'], mag))

        c_11 = 'swarp {}/mag_{}_CCD_x?_y?_d?.fits'.format(prfs_d['fits_ref'],
                                                          mag)
        c_12 = ' -IMAGEOUT_NAME {} -WEIGHTOUT_NAME {}'.format(coadd_loc,
                                                              coadd_weight_loc)
        c_13 = ' -WEIGHT_TYPE NONE -VERBOSE_TYPE QUIET'
        c_14 = ' -GAIN_KEYWORD GAIN'
        c_1 = c_11 + c_12 + c_13 + c_14

        process_1 = Popen(c_1, shell=True)
        process_1.wait()
        logger.info('Swarp process finished')

        folder_sex = '{}_{}_{}_{}_{}'.format(analysis_d['deblend_nthresh'],
                                             analysis_d['analysis_thresh'],
                                             analysis_d['detect_thresh'],
                                             analysis_d['deblend_mincount'],
                                             analysis_d['detect_minarea'])

        # Check if folder exists
        conf_folder = '{}/{}'.format(prfs_d['fits_ref'], folder_sex)
        try:
            if not path.isfile(conf_folder):
                mkdir(conf_folder)
        except OSError:
            logger.debug('folder {} created'.format(conf_folder))

        cat_loc = '{}/catalog_{}.cat'.format(conf_folder, mag)

        logger.info('Swarped image sextraction process')
        c_21 = 'sex -c {} {}/coadd_{}.fits'.format(prfs_d['conf_sex'],
                                                   prfs_d['fits_ref'], mag)
        c_22 = ' -CATALOG_NAME {}'.format(cat_loc)
        c_23 = ' -PARAMETERS_NAME {}'.format(prfs_d['params_cat'])
        c_24 = ' -DETECT_MINAREA {}'.format(analysis_d['detect_minarea'])
        c_25 = ' -DETECT_THRESH {}'.format(analysis_d['detect_thresh'])
        c_26 = ' -ANALYSIS_THRESH {}'.format(analysis_d['analysis_thresh'])
        c_27 = ' -DEBLEND_NTHRESH {}'.format(analysis_d['deblend_nthresh'])
        c_28 = ' -DEBLEND_MINCONT {}'.format(analysis_d['deblend_mincount'])
        c_2 = c_21 + c_22 + c_23 + c_24 + c_25 + c_26 + c_27 + c_28

        process_2 = Popen(c_2, shell=True)
        process_2.wait()
        logger.info('Sextractor process for catalogue finished')
        logger.info('Catalog {} created'.format(cat_loc))

        remove(coadd_loc)
        logger.debug('{} removed'.format(coadd_loc))
        remove(coadd_weight_loc)
        logger.debug('{} removed'.format(coadd_weight_loc))

        return True
