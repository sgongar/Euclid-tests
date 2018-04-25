#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for sextractor routines

Conventions:
    * Single n means name.
    * Single d means dictionary.
    * Single j means jobs.
    * Single p means process.
    * loc means location.

Versions:
- 0.1: Initial release.

Todo:
    * Improve documentation
    * Improve log messages

*GNU Terry Pratchett*
"""

from os import listdir
from subprocess import Popen

from multiprocessing import Process

from misc import extract_settings_sc3


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.2"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class SextractorSC3:

    def __init__(self, logger, analysis_d):
        """

        """
        self.prfs_d = extract_settings_sc3()
        self.sextractor_process(logger, analysis_d)

    def sextractor_process(self, logger, analysis_d):
        """

        :param logger:
        :param analysis_d:
        :return:
        """
        logger.info('Starting sextractor process for fits images')

        fits_files = listdir(self.prfs_d['fits_dir'])
        active_sex = []

        for image_idx in range(0, len(fits_files),
                               self.prfs_d['cores_number']):
            try:
                sex_j = []
                for proc in range(0, self.prfs_d['cores_number'], 1):
                    idx = image_idx + proc  # index

                    sex_file = fits_files[idx]
                    cat_name = '{}.cat'.format(fits_files[idx][:-5])

                    # sextractor input and output
                    sex_input = '{}/{}'.format(self.prfs_d['fits_dir'],
                                               sex_file)
                    sex_output = '{}/{}'.format(self.prfs_d['fits_dir'],
                                                cat_name)

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