#!/usr/bin/python
# -*- coding: utf-8 -*-


"""Script principal del paquete.

Versions:
- 0.1: Initial release. Splited from check.py
       Recreated for SC3 analysis pipeline.

Todo:
    *

*GNU Terry Pratchett*

"""
from itertools import product
from multiprocessing import Process
from sys import argv

from images_management_sc3 import extract_images, extract_flags
from misc import setting_logger, extract_settings_sc3
from misc import create_configurations, get_fits_sc3
from misc import create_sextractor_dict, create_scamp_dict
from sextractor_aux import SextractorSC3
from scamp_aux import Scamp
from scamp_filter_sc3 import ScampFilter
from cosmic_sc3 import cosmic

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def scamp_f_name(idx):
    """

    :param idx:
    :return: scmp_d, scmp_cf
    """
    scmp_d, len_confs = create_scamp_dict(idx)

    scmp_cf = '{}_{}_{}_{}'.format(scmp_d['crossid_radius'],
                                   scmp_d['pixscale_maxerr'],
                                   scmp_d['posangle_maxerr'],
                                   scmp_d['position_maxerr'])

    return scmp_d, scmp_cf


class Check:

    def __init__(self):
        """

        """
        self.logger = setting_logger()
        self.prfs_d = extract_settings_sc3()

        # Scamp configurations.
        mode = {'type': 'scamp'}
        self.scamp_confs, self.scamp_confs_n = create_configurations(mode)
        # Sextractor configurations.
        mode = {'type': 'sextractor'}
        self.sex_confs, sex_confs_n = create_configurations(mode)

        if argv[1] == '-full':
            if not self.full_pipeline():
                raise Exception
        if argv[1] == '-catalog':
            if not self.catalog():
                raise Exception
        elif argv[1] == '-split':
            if not self.split():
                raise Exception
        elif argv[1] == '-sextractor':
            if not self.sextractor():
                raise Exception
        elif argv[1] == '-scamp':
            if not self.scamp():
                raise Exception
        elif argv[1] == '-filter':
            self.filt()

    def full_pipeline(self):
        """

        :return:
        """
        if not self.sextractor():
            raise Exception
        if not self.scamp():
            raise Exception
        if not self.filt():
            raise Exception

        return True

    def catalog(self):
        """

        :return:
        """
        # Opens input_catalog

        # Gets simplified header
        # Merges header and data
        # Saves file

    def split(self):
        """

        :return:
        """
        fits_list = get_fits_sc3()

        for fits_ in fits_list:
            extract_images(fits_)
            extract_flags(fits_)
            # extract_rms (?)

        return True

    def clean(self):
        """

        :return: True if everything goes alright
        """

    def sextractor(self):
        """

        :return: True if everything goes alright
        """
        mode = {'type': 'sextractor'}
        confs, total_confs = create_configurations(mode)

        for idx, conf_ in enumerate(confs):
            analysis_d, len_dicts = create_sextractor_dict(idx, False)
            # Just for tests reasons
            if len_dicts != total_confs:
                raise Exception
            # todo - implement a return!
            SextractorSC3(self.logger, analysis_d)

            self.logger.debug('Performs an analysis over a bunch of files')

        return True

    def scamp(self):
        """ todo - improve docstring
            todo - improve return
        Tip. scamp already parallelize its own process so there is no need
        to implement any multiprocess.

        :return:
        """
        mode = {'type': 'sextractor'}
        confs_sex, total_confs = create_configurations(mode)

        for idx_scmp, conf_scmp in enumerate(self.scamp_confs):
            # todo - check if everything is alright
            scmp_d, scmp_cf = scamp_f_name(idx_scmp)
            for idx_sex, conf_sex in enumerate(confs_sex):
                analysis_d, len_dicts = create_sextractor_dict(idx_sex,
                                                               False)
                if not Scamp(self.logger, scmp_d, scmp_cf,
                             analysis_d):
                    raise Exception  # todo improve Exception

        return True

    def filt(self):
        """ todo - improve docstring
            todo - improve return
        :return:
        """

        confs = list(product(self.sex_confs, self.scamp_confs))

        for idx, conf_ in enumerate(confs):
            filt_j = []
            # while len(filt_j) < self.prfs_d['cores_number'] + 1:
            while len(filt_j) < 1:
                sex_d = {'deblend_mincount': conf_[0][1],
                         'analysis_thresh': conf_[0][2],
                         'detect_thresh': conf_[0][2],
                         'deblend_nthresh': conf_[0][0],
                         'detect_minarea': conf_[0][3],
                         'filter': 'models/gauss_2.0_5x5.conv'}

                scmp_cf = '{}_{}_{}_{}'.format(conf_[1][0], conf_[1][1],
                                               conf_[1][2], conf_[1][3])
                filt_p = Process(target=ScampFilter,
                                 args=(self.logger, scmp_cf, sex_d,))
                filt_j.append(filt_p)
                filt_p.start()

            active_filt = list([j.is_alive() for j in filt_j])
            while True in active_filt:
                active_filt = list([j.is_alive() for j in filt_j])
                pass

        return True


if __name__ == '__main__':
    check_process = Check()
