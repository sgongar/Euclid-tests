#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.
- 0.2: Functions reorganised to two different classes.
- 0.3: New organization for filter Class
- 0.4: Filter now gets areas
- 0.5: check_source moved to another file
- 0.6: Filter process get out - Class for scamp process in SC3 files added
- 0.7: Scamp process for sc3 images moved to scamp_aux_sc3.py file

Todo:
    * Improve documentation

*GNU Terry Pratchett*

"""
from subprocess import Popen

from misc import extract_settings_sc3

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.7"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class ScampSC3:

    def __init__(self, logger, scmp_d):
        """

        :param logger:
        :param scmp_d:
        """
        self.prfs_d = extract_settings_sc3()

        self.logger = logger
        self.scmp_d = scmp_d

        self.scamp_process()

    def scamp_process(self):
        """

        :return:
        """
        self.logger.info('Scamp process')

        scmp_1 = 'scamp -c {}'.format(self.prfs_d['conf_scamp'])
        scmp_2 = ' {}/*.cat'.format(self.prfs_d['fits_dir'])
        scmp_3 = ' -ASTREFCAT_NAME {}/cat.cat'.format(self.prfs_d['references'])
        scmp_4 = ' -PIXSCALE_MAXERR {}'.format(self.scmp_d['pixscale_maxerr'])
        scmp_5 = ' -POSANGLE_MAXERR {}'.format(self.scmp_d['posangle_maxerr'])
        scmp_6 = ' -POSITION_MAXERR {}'.format(self.scmp_d['position_maxerr'])
        scmp_7 = ' -CROSSID_RADIUS {}'.format(self.scmp_d['crossid_radius'])
        # Output catalogs location
        merged_cat_n = '{}/merged.cat'.format(self.prfs_d['output_cats'])
        scmp_8 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat_n)
        full_cat_n = '{}/full.cat'.format(self.prfs_d['output_cats'])
        scmp_9 = ' -FULLOUTCAT_NAME {}'.format(full_cat_n)
        scmp_p = scmp_1 + scmp_2 + scmp_3 + scmp_4 + scmp_5
        scmp_p = scmp_p + scmp_6 + scmp_7 + scmp_8 + scmp_9
        scmp_p = scmp_p

        process_scamp = Popen(scmp_p, shell=True)
        process_scamp.wait()

        self.logger.info('Scamp process finished.')

        return True
