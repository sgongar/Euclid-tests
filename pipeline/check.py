#!/usr/bin/python
# -*- coding: utf-8 -*-


"""Script principal del paquete.

Con este script podremos lanzar las acciones de analsis necesarias.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * Improve log messages

"""

from multiprocessing import Process
from os import path
from sys import argv

from cats_management import look_for_ssos
from cats_management import merge_sso_cat, merge_ssos
from stats_management import merge_stats
from curves import Detections_vs_Tolerance, Detections_vs_Created
from curves import PMs_vs_Detections
from errors import BadSettings
from errors import CatalogueError
from misc import setting_logger, extract_settings
from misc import create_configurations, pipeline_help
from misc import create_sextractor_dict, create_scamp_dict
from performance import SextractorPerformance
from sextractor_aux import Sextractor, catalogue_creation
from scamp_aux import Scamp, scamp_filter
from stats_management import ExtractStats

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


class Check:

    def __init__(self):
        logger = setting_logger()
        prfs_d = extract_settings()
        mode = {'type': 'scamp'}
        confs, total_confs = create_configurations(logger, prfs_d, mode)

        if argv[1] == '-full':
            if not self.full_pipeline(logger, prfs_d, mode,
                                      confs, total_confs):
                raise Exception
        elif argv[1] == '-catalog':
            if not self.catalog(prfs_d):
                raise CatalogueError("Catalog couldn't be created")
        elif argv[1] == '-sextractor':
            self.sextractor(logger, prfs_d)
        elif argv[1] == '-scamp':
            self.scamp(logger, prfs_d, mode, confs)
        elif argv[1] == '-filter':
            self.filt(logger, prfs_d, mode, confs, total_confs)
        elif argv[1] == '-check':
            self.check(logger, prfs_d, confs, total_confs)
        elif argv[1] == '-stats':
            self.stats(logger, prfs_d, confs, total_confs)
        elif argv[1] == '-plot':
            self.plot(logger)
        elif argv[1] == '-help':
            pipeline_help(logger)
        elif argv[1] == '-sextractor_performance':
            self.sextractor_performance(logger, prfs_d)
        elif argv[1] == '-scamp_performance':
            self.scamp_performance(logger, prfs_d)
        else:
            raise BadSettings('not analysis option choosen')

    def full_pipeline(self, logger, prfs_d, mode, confs, total_confs):
        """ run a full pipeline over all files

        @param logger:
        @param prfs_d:
        @param mode:
        @param confs:
        @param total_confs:


        @return True: if everything goes alright
        """
        """
        if not self.catalogue(logger):
            raise Exception
        if not self.sextractor(logger, prfs_d):
            raise Exception
        """
        if not self.scamp(logger, prfs_d, mode, confs):
            raise Exception
        if not self.filt(logger, prfs_d, mode, confs, total_confs):
            raise Exception
        if not self.check(logger, prfs_d, confs, total_confs):
            raise Exception
        if not self.stats(logger, prfs_d, confs, total_confs):
            raise Exception

        return True

    def catalog(self, logger, prfs_d):
        """ creates a catalog.

        @param logger: a logger object.
        @param prfs_d:

        @return True: if everything goes alright.
        """
        analysis_dict, len_dicts = create_sextractor_dict(logger, prfs_d,
                                                          0, True)
        # Catalogue creation. Only created one time.
        if not catalogue_creation(logger, prfs_d, analysis_dict):
            raise CatalogueError("Catalog couldn't be created")

        return True

    def sextractor(self, logger, prfs_d):
        """

        @param logger: a logger object.
        @param prfs_d:

        @return True: if everything goes alright.
        """
        analysis_dir = prfs_d['fits_dir']
        mode = {'type': 'sextractor'}
        confs, total_confs = create_configurations(logger, prfs_d, mode)

        confs = [1]

        # Tries to analyse a custom image with custom output
        # If no argument is passed to performs a regular analysis over
        # all files in FITS folder
        try:
            sextractor_file = argv[2]
            sextractor_output = argv[3]
            """
            logger.debug('perfoming analysis over a single file')
            if not sextractor_thread(logger, prfs_d, sextractor_file,
                                     sextractor_output, analysis_dict,
                                     analysis_dir):
                raise Exception  # TODO Create a custom Exception
            """
        # If no arguments are passed a complete analysis will be perform
        except IndexError:
            for mag in prfs_d['mags']:
                for idx, conf_ in enumerate(confs):
                    (analysis_d,
                     len_dicts) = create_sextractor_dict(logger, prfs_d,
                                                         idx, cat_conf=False)
                    # Just for tests reasons
                    if len_dicts != total_confs:
                        raise Exception

                    # Sextractor extraction parameters - Uncomment if necessary
                    analysis_d = {'deblend_mincount': 0.1,
                                  'analysis_thresh': 1.35,
                                  'detect_thresh': 1.35,
                                  'deblend_nthresh': 2, 'detect_minarea': 4,
                                  'filter': 'models/gauss_2.0_5x5.conv'}

                    Sextractor(logger, prfs_d, analysis_d, analysis_dir,
                               regular=True)

            logger.debug('perfoming analysis over a bunch of files')

        return True

    def scamp(self, logger, prfs_d, mode, confs_scmp):
        """

        @param logger: a logger object.
        @param prfs_d:
        @param mode:
        @param confs:

        Tip. scamp already parallelize its own process so there is no need
        to implement any multiprocess.

        @return True: if everything goes alright.
        """
        mode = {'type': 'sextractor'}
        confs_sex, total_confs = create_configurations(logger, prfs_d, mode)

        confs_sex = [1]

        for mag in prfs_d['mags']:
            for idx_scmp, conf_scmp in enumerate(confs_scmp):
                scmp_d, len_confs = create_scamp_dict(logger, prfs_d, idx_scmp)
                f_conf = '{}_{}_{}_{}'.format(scmp_d['crossid_radius'],
                                              scmp_d['pixscale_maxerr'],
                                              scmp_d['posangle_maxerr'],
                                              scmp_d['position_maxerr'])
                print "scamp"
                for idx_sex, conf_sex in enumerate(confs_sex):
                    print "sextractor"
                    (analysis_d,
                     len_dicts) = create_sextractor_dict(logger, prfs_d,
                                                         idx_sex,
                                                         cat_conf=False)
                    analysis_d = {'deblend_mincount': 0.1,
                                  'analysis_thresh': 1.35,
                                  'detect_thresh': 1.35,
                                  'deblend_nthresh': 2, 'detect_minarea': 4,
                                  'filter': 'models/gauss_2.0_5x5.conv'}
                    if not Scamp(logger, mag, scmp_d, f_conf, analysis_d):
                        raise Exception

        return True

    def filt(self, logger, prfs_d, mode, confs, total_confs):
        """ Performs a complete pipeline to scamp output.

        @param logger:
        @param prfs_d:
        @param mode:
        @param confs:
        @param total_confs:

        @return True if everything goes alright.
        """
        filt_j = []

        for mag in prfs_d['mags']:
            for idx_filt in range(0, total_confs, prfs_d['cores_number']):
                try:
                    filt_j = []
                    idx_proc = 0
                    while len(filt_j) < prfs_d['cores_number']:
                        idx = idx_filt + idx_proc
                        (scmp_d, len_confs) = create_scamp_dict(logger,
                                                                prfs_d, idx)
                        conf = [scmp_d['crossid_radius'],
                                scmp_d['pixscale_maxerr'],
                                scmp_d['posangle_maxerr'],
                                scmp_d['position_maxerr']]
                        f_conf = '{}_{}_{}_{}'.format(conf[0], conf[1],
                                                      conf[2], conf[3])
                        f_name = 'results/filt_{}_20-21__6.csv'.format(f_conf)
                        if not path.isfile(f_name):
                            filt_p = Process(target=scamp_filter,
                                             args=(logger, prfs_d,
                                                   mag, scmp_d,))
                            filt_j.append(filt_p)
                            filt_p.start()
                        idx_proc += 1
                    active_filt = list([j.is_alive() for j in filt_j])
                    while True in active_filt:
                        active_filt = list([j.is_alive() for j in filt_j])
                        pass
                    filt_j = []
                except IndexError:
                    print "finished"

        return True

    def check(self, logger, prfs_d, confs, total_confs):
        """

        @param logger:
        @param prfs_d:
        @param confs:
        @param total_confs:
        """
        check_j = []

        for mag in prfs_d['mags']:
            for idx_conf in range(0, total_confs, 2):
                try:
                    check_j = []
                    idx_proc = 0
                    while len(check_j) < 2:
                        idx = idx_conf + idx_proc
                        (scmp_d, len_confs) = create_scamp_dict(logger,
                                                                prfs_d, idx)

                        conf = [scmp_d['crossid_radius'],
                                scmp_d['pixscale_maxerr'],
                                scmp_d['posangle_maxerr'],
                                scmp_d['position_maxerr']]
                        f_conf = '{}_{}_{}_{}'.format(conf[0], conf[1],
                                                      conf[2], conf[3])
                        """
                        f_name_1 = 'm_20-21_x2_y2'  # TODO Improve!
                        f_name_2 = '20-21_sso_cat.csv'
                        f_name = '{}{}{}'.format(f_name_1, f_conf, f_name_2)

                        if not path.isfile(f_name):
                        """
                        check_p = Process(target=look_for_ssos,
                                          args=(logger, prfs_d,
                                                mag, scmp_d, f_conf, ))
                        check_j.append(check_p)
                        check_p.start()

                        idx_proc += 1

                    active_check = list([j.is_alive() for j in check_j])
                    while True in active_check:
                        active_check = list([j.is_alive() for j in check_j])
                        pass
                    check_j = []
                except IndexError:
                    print "finished"

        # Merge catalog files into a single one and cleans old ones
        for mag in prfs_d['mags']:
            for conf in confs:
                merge_ssos(logger, prfs_d, conf, mag)
                merge_sso_cat(logger, prfs_d, conf, mag)

        return True

    def stats(self, logger, prfs_d, confs, total_confs):
        """

        @param logger:
        @param prfs_d:
        @param confs:
        @param total_confs:


        @return True: if everything goes alright.
        """

        for mag_ in prfs_d['mags']:
            for idx_stats in range(0, total_confs, prfs_d['cores_number']):
                try:
                    stats_j = []
                    idx_proc = 0
                    while len(stats_j) < prfs_d['cores_number']:
                        idx = idx_stats + idx_proc
                        logger.debug('creating configuration parameters')
                        (scmp_d, len_confs) = create_scamp_dict(logger,
                                                                prfs_d, idx)

                        conf = [scmp_d['crossid_radius'],
                                scmp_d['pixscale_maxerr'],
                                scmp_d['posangle_maxerr'],
                                scmp_d['position_maxerr']]

                        f_conf = '{}_{}_{}_{}'.format(conf[0], conf[1],
                                                      conf[2], conf[3])
                        f_name = 'stats_{}_{}.csv'.format(f_conf, mag_)

                        if not path.isfile(f_name):
                            stats_p = Process(target=ExtractStats,
                                              args=(logger, mag_, scmp_d,
                                                    f_conf, filt,))
                            stats_j.append(stats_p)
                            stats_p.start()

                        idx_proc += 1

                    active_stats = list([j.is_alive() for j in stats_j])
                    while True in active_stats:
                        active_stats = list([j.is_alive() for j in stats_j])
                        pass
                    stats_j = []
                except Exception as e:
                    print e

        if not merge_stats(logger, prfs_d):
            raise Exception

        return True

    def plot(self, logger):
        """

        @param logger:
        """

        plot_j = []
        """
        plot_p = Process(target=self.plot_a_thread_)
        plot_j.append(plot_p)
        plot_p.start()

        plot_p = Process(target=self.plot_b_thread_)
        plot_j.append(plot_p)
        plot_p.start()
        """
        plot_p = Process(target=self.plot_c_thread_, args=(logger,))
        plot_j.append(plot_p)
        plot_p.start()

        active_plot = list([j.is_alive() for j in plot_j])
        while True in active_plot:
            active_plot = list([j.is_alive() for j in plot_j])
            pass

        return True

    def plot_a_thread_(self):
        print "dentro de plot_thread_1"
        pdf_plot = Detections_vs_Tolerance()

    def plot_b_thread_(self):
        print "dentro de plot_thread_2"
        pdf_plot = Detections_vs_Created()

    def plot_c_thread_(self, logger):
        logger.debug('plotting plot c')
        pdf_plot = PMs_vs_Detections()

    def sextractor_performance(self, logger, prfs_d):
        sextractor_performance = SextractorPerformance(logger, prfs_d)

    def scamp_perfomance(self, logger, prfs_d):
        pass

if __name__ == '__main__':
    check_process = Check()
