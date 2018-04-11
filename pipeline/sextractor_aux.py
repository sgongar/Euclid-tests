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

from os import remove
from subprocess import Popen

from multiprocessing import Process

from misc import extract_settings, get_fits, create_folder


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class SextractorSC3:

    def __init__(self, logger, analysis_d):
        """

        """
        self.prfs_d = extract_settings()
        self.sextractor_process(logger, analysis_d)

    def sextractor_process(self, logger, analysis_d):
        """

        :param logger:
        :param analysis_d:
        :return:
        """
        logger.info('Starting sextractor process for fits images')

        # TODO hardcoded!
        dither_names = ['EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00_',
                        'EUC_VIS_SWL-DET-002-000000-0000000__20170630T011642.0Z_00.00_',
                        'EUC_VIS_SWL-DET-003-000000-0000000__20170630T011848.6Z_00.00_',
                        'EUC_VIS_SWL-DET-004-000000-0000000__20170630T012050.1Z_00.00_']
        analysis_dir = '/pcdisk/holly/sgongora/Dev'  # TODO hardcoded!
        cores_number = 24

        fits_files = []
        # flag_files = []
        active_sex = []
        for idx in range(0, 36, 1):
            for dither in dither_names:
                fits_files.append('{}{}.fits'.format(dither, idx))
                # flag_files.append('{}f{}.fits'.format(dither, idx))

        for image_idx in range(0, len(fits_files), cores_number):
            try:
                sex_j = []
                for proc in range(0, cores_number, 1):
                    idx = image_idx + proc  # index

                    sex_file = fits_files[idx]
                    # sex_flag = flag_files[idx]
                    folder_loc = '{}/CCDs'.format(analysis_dir)
                    cat_name = '{}.cat'.format(fits_files[idx][:-5])

                    # sextractor input and output
                    sex_input = '{}/{}'.format(folder_loc, sex_file)
                    sex_output = '{}/{}'.format(folder_loc, cat_name)
                    # sex_flag = '{}/{}'.format(folder_loc, sex_flag)

                    sex_p = Process(target=self.sextractor_thread,
                                    args=(sex_input, sex_output, analysis_d))
                    sex_j.append(sex_p)
                    sex_p.start()

                    active_sex = list([job.is_alive() for job in sex_j])
                while True in active_sex:
                    active_sex = list([job.is_alive() for job in sex_j])
                    pass
            except IndexError:
                print('Extraction finished')

        return True

    def sextractor_thread(self, sextractor_file, sextractor_output,
                          analysis_d):
        """ runs sextractor on a single file
        todo - improve docstring

        :param sextractor_file: file to be 'sextracted'
        :param sextractor_output: catalog to be created by sextractor
        :param analysis_d: parameters dict
        :return: if everything goes alright
        """

        s_1 = 'sex -c {} {}'.format(self.prfs_d['conf_sex'], sextractor_file)
        s_2 = ' -CATALOG_NAME {}'.format(sextractor_output)
        s_3 = ' -PARAMETERS_NAME {}'.format(self.prfs_d['params_sex'])
        s_4 = ' -STARNNW_NAME {}'.format(self.prfs_d['neural_sex'])
        s_5 = ' -DETECT_MINAREA {}'.format(analysis_d['detect_minarea'])
        s_6 = ' -DETECT_THRESH {}'.format(analysis_d['detect_thresh'])
        s_7 = ' -ANALYSIS_THRESH {}'.format(analysis_d['analysis_thresh'])
        s_8 = ' -DEBLEND_NTHRESH {}'.format(analysis_d['deblend_nthresh'])
        s_9 = ' -DEBLEND_MINCONT {}'.format(analysis_d['deblend_mincount'])
        s_10 = ' -FILTER_NAME {}'.format(analysis_d['filter'])

        cmd_3 = s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10

        sextractor_p = Popen(cmd_3, shell=True)
        sextractor_p.wait()

        return True


class Sextractor:

    def __init__(self, logger, analysis_d):
        """

        """
        self.prfs_d = extract_settings()
        self.sextractor_process(logger, analysis_d)

    def sextractor_process(self, logger, analysis_d):
        """

        :param logger:
        :param analysis_d:
        :return:
        """
        logger.info('Starting sextractor process for fits images')

        analysis_dir = self.prfs_d['fits_dir']

        folder_n = '{}_{}_{}_{}_{}'.format(analysis_d['deblend_nthresh'],
                                           analysis_d['analysis_thresh'],
                                           analysis_d['detect_thresh'],
                                           analysis_d['deblend_mincount'],
                                           analysis_d['detect_minarea'])

        logger.info('sextractor configuration {}'.format(folder_n))

        print('mags {}'.format(self.prfs_d['mags']))

        for mag_ in self.prfs_d['mags']:
            fits_files = get_fits(unique=False, mag=mag_)
            for image_idx in range(0, len(fits_files),
                                   self.prfs_d['cores_number']):
                try:
                    sex_j = []
                    for proc in range(0, self.prfs_d['cores_number'], 1):
                        idx = image_idx + proc  # index

                        sex_file = fits_files[idx]
                        folder_loc = '{}/{}/CCDs/{}'.format(analysis_dir, mag_,
                                                            folder_n)
                        create_folder(logger, folder_loc)
                        cat_name = '{}.cat'.format(fits_files[idx][:-5])

                        # sextractor input and output
                        sex_input = '{}/{}/CCDs/{}'.format(analysis_dir, mag_,
                                                           sex_file)
                        sex_output = '{}/{}'.format(folder_loc, cat_name)

                        sex_p = Process(target=self.sextractor_thread,
                                        args=(sex_input, sex_output,
                                              analysis_d,))
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

    def sextractor_thread(self, sextractor_file, sextractor_output,
                          analysis_d):
        """ runs sextractor on a single file
        todo - improve docstring

        :param sextractor_file: file to be 'sextracted'
        :param sextractor_output: catalog to be created by sextractor
        :param analysis_d: parameters dict
        :return: if everything goes alright
        """

        s_1 = 'sex -c {} {}'.format(self.prfs_d['conf_sex'], sextractor_file)
        s_2 = ' -CATALOG_NAME {}'.format(sextractor_output)
        s_3 = ' -PARAMETERS_NAME {}'.format(self.prfs_d['params_sex'])
        s_4 = ' -STARNNW_NAME {}'.format(self.prfs_d['neural_sex'])
        s_5 = ' -DETECT_MINAREA {}'.format(analysis_d['detect_minarea'])
        s_6 = ' -DETECT_THRESH {}'.format(analysis_d['detect_thresh'])
        s_7 = ' -ANALYSIS_THRESH {}'.format(analysis_d['analysis_thresh'])
        s_8 = ' -DEBLEND_NTHRESH {}'.format(analysis_d['deblend_nthresh'])
        s_9 = ' -DEBLEND_MINCONT {}'.format(analysis_d['deblend_mincount'])
        s_10 = ' -FILTER_NAME {}'.format(analysis_d['filter'])

        cmd_3 = s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10

        print(cmd_3)

        # sextractor_p = Popen(cmd_3, shell=True)
        # sextractor_p.wait()

        return True


class CatalogCreation:

    def __init__(self, logger, analysis_d):
        """ creates a single catalogue for all images available

        :param logger: a logger object
        :param analysis_d: a dictionary with all configuration information
        """
        self.prfs_d = extract_settings()

        cat_j = []
        for mag_ in self.prfs_d['mags']:
            cat_p = Process(target=self.catalogue_thread,
                            args=(logger, analysis_d, mag_,))
            cat_j.append(cat_p)
            cat_p.start()

        active_cat = list([job.is_alive() for job in cat_j])
        while True in active_cat:
            active_cat = list([job.is_alive() for job in cat_j])
            pass

    def catalogue_thread(self, logger, analysis_d, mag):
        """

        :param logger:
        :param analysis_d:
        :param mag:
        :return:
        """
        # Fits files location
        coadd_loc = '{}/coadd_{}.fits'.format(self.prfs_d['fits_ref'], mag)
        coadd_w_loc = '{}/coadd_{}.w.fits'.format(self.prfs_d['fits_ref'], mag)

        # Harcoded for single CCD catalog!
        logger.info('swarp process launched for mag {}'.format(mag))
        logger.info('files {}/mag_{}_*.fits'.format(self.prfs_d['fits_ref'],
                                                    mag))

        ccd_names = 'CCD_x?_y?_Stars.fits'.format(mag)
        c_11 = 'swarp {}/{}'.format(self.prfs_d['fits_ref'], ccd_names)
        c_12 = ' -IMAGEOUT_NAME {} -WEIGHTOUT_NAME {}'.format(coadd_loc,
                                                              coadd_w_loc)
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
        conf_folder = '{}/{}'.format(self.prfs_d['fits_ref'], folder_sex)
        create_folder(logger, conf_folder)

        cat_loc = '{}/catalog_{}.cat'.format(conf_folder, mag)

        logger.info('Swarped image sextraction process')
        c_21 = 'sex -c {} {}/coadd_{}.fits'.format(self.prfs_d['conf_sex'],
                                                   self.prfs_d['fits_ref'],
                                                   mag)
        c_22 = ' -CATALOG_NAME {}'.format(cat_loc)
        c_23 = ' -PARAMETERS_NAME {}'.format(self.prfs_d['params_cat'])
        c_24 = ' -DETECT_MINAREA {}'.format(analysis_d['detect_minarea'])
        c_25 = ' -DETECT_THRESH {}'.format(analysis_d['detect_thresh'])
        c_26 = ' -ANALYSIS_THRESH {}'.format(analysis_d['analysis_thresh'])
        c_27 = ' -DEBLEND_NTHRESH {}'.format(analysis_d['deblend_nthresh'])
        c_28 = ' -DEBLEND_MINCONT {}'.format(analysis_d['deblend_mincount'])
        c_2 = c_21 + c_22 + c_23 + c_24 + c_25 + c_26 + c_27 + c_28

        print('c_2 {}'.format(c_2))
        process_2 = Popen(c_2, shell=True)
        process_2.wait()
        logger.info('Sextractor process for catalogue finished')
        logger.info('Catalog {} created'.format(cat_loc))

        remove(coadd_loc)
        logger.debug('{} removed'.format(coadd_loc))
        remove(coadd_w_loc)
        logger.debug('{} removed'.format(coadd_w_loc))

        return True
