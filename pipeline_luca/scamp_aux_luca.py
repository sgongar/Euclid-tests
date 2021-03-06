#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.

Todo:
    * Improve documentation
    * sex_cf and scmp_cf should come from main file

*GNU Terry Pratchett*
"""
from subprocess import Popen

from misc import create_folder, extract_settings

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
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

        self.logger = logger
        self.mag = mag
        self.scmp_d = scmp_d
        self.scmp_cf = scmp_cf
        self.sex_d = sex_d  # TODO

        self.scamp_process()

    def scamp_process(self):
        """

        :return:
        """
        sex_cf = '{}_{}_{}_{}_{}'.format(self.sex_d['deblend_nthresh'],
                                         self.sex_d['analysis_thresh'],
                                         self.sex_d['detect_thresh'],
                                         self.sex_d['deblend_mincount'],
                                         self.sex_d['detect_minarea'])

        self.logger.info('Scamp process for magnitude {}'.format(self.mag))
        self.logger.info('Sextractor configuration: {}'.format(sex_cf))
        self.logger.info('Scamp configuration: {}'.format(self.scmp_cf))

        scmp_1 = 'scamp -c {}'.format(self.prfs_d['conf_scamp'])
        sex_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], self.mag,
                                         sex_cf)
        sex_output = 'mag_{}_CCD_x?_y?_d?.cat'.format(self.mag)
        scmp_2 = ' {}/{}'.format(sex_loc, sex_output)
        scmp_3 = ' -ASTREFCAT_NAME'
        cat_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], self.mag,
                                         sex_cf)
        cat_input = 'catalog_{}.cat'.format(self.mag)
        scmp_4 = ' {}/{}'.format(cat_loc, cat_input)
        scmp_5 = ' -PIXSCALE_MAXERR {}'.format(self.scmp_d['pixscale_maxerr'])
        scmp_6 = ' -POSANGLE_MAXERR {}'.format(self.scmp_d['posangle_maxerr'])
        scmp_7 = ' -POSITION_MAXERR {}'.format(self.scmp_d['position_maxerr'])
        scmp_8 = ' -CROSSID_RADIUS {}'.format(self.scmp_d['crossid_radius'])
        # Output catalogs location
        cats_dir = '{}/{}/{}/{}'.format(self.prfs_d['catalogs_dir'], self.mag,
                                        sex_cf, self.scmp_cf)
        merged_cat = '{}/merged_{}_{}.cat'.format(cats_dir, self.scmp_cf,
                                                  self.mag)
        scmp_9 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat)
        full_cat = '{}/full_{}_{}.cat'.format(cats_dir, self.scmp_cf, self.mag)
        scmp_10 = ' -FULLOUTCAT_NAME {}'.format(full_cat)
        scmp_p = scmp_1 + scmp_2 + scmp_3 + scmp_4 + scmp_5
        scmp_p = scmp_p + scmp_6 + scmp_7 + scmp_8 + scmp_9
        scmp_p = scmp_p + scmp_10

        create_folder(self.logger, cats_dir)

        process_scamp = Popen(scmp_p, shell=True)
        process_scamp.wait()

        self.logger.info('Scamp process finished.')

        return True
