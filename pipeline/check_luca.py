#!/usr/bin/python
# -*- coding: utf-8 -*-


"""Script principal del paquete.

Versions:
- 0.1: Initial release.
- 0.2: Plot functions were removed.
- 0.3:
- 0.4: SextractorSizes added.
- 0.5: Changed reference to ScampFilter

Todo:
    *

*GNU Terry Pratchett*

"""
from itertools import product
from multiprocessing import Process
from sys import argv

from cats_management import look_for_ssos
from cats_management import merge_sso_cat, merge_ssos
from errors import BadSettings
from errors import CatalogueError
from fitting_pm_mag import FitPMMagAgainstSizes
from misc import setting_logger, extract_settings
from misc import create_configurations, pipeline_help
from misc import create_sextractor_dict, create_scamp_dict
from pandas import DataFrame
from sextractor_performance import SextractorPerformance
from sextractor_sizes import SextractorSizes
from scamp_performance_stars import ScampPerformanceStars
from scamp_performance_ssos import ScampPerformance, TotalScampPerformance
from performance import PMPerformance, StatsPerformance
from pm_performance_ssos import SlowPMPerformanceSSOs
from sextractor_aux import Sextractor, CatalogCreation
from scamp_aux import Scamp
from scamp_filter import ScampFilter
from stats_management import ExtractStats

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.4"
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
        self.prfs_d = extract_settings()

        # Scamp configurations.
        mode = {'type': 'scamp'}
        self.scamp_confs, self.scamp_confs_n = create_configurations(mode)
        # Sextractor configurations.
        mode = {'type': 'sextractor'}
        self.sex_confs, sex_confs_n = create_configurations(mode)

        if argv[1] == '-full':
            if not self.full_pipeline():
                raise Exception
        elif argv[1] == '-catalog':
            if not self.catalog():
                raise CatalogueError("Catalog couldn't be created")
        elif argv[1] == '-sextractor':
            if not self.sextractor():
                raise Exception
        elif argv[1] == '-scamp':
            if not self.scamp():
                raise Exception
        elif argv[1] == '-filter':
            self.filt()
        # elif argv[1] == '-check':
        #     self.check(self.scamp_confs, self.scamp_confs_n)
        elif argv[1] == '-stats':
            self.stats(self.scamp_confs_n)
        # elif argv[1] == '-error_performance':
        #     self.error_performance()
        elif argv[1] == '-scamp_performance_ssos':
            self.scamp_performance_ssos()
        elif argv[1] == '-scamp_performance_stars':
            self.scamp_performance_stars()
        elif argv[1] == '-fit_pm_mag':
            self.fit_pm_mag()
        elif argv[1] == '-scamp_performance_fitting':
            self.scamp_performance_fitting()
        elif argv[1] == '-pm_performance_ssos':
            self.pm_performance_ssos()
        elif argv[1] == '-sextractor_size':
            self.sextractor_size()
        elif argv[1] == '-pm_performance':
            self.pm_performance()
        elif argv[1] == '-stats_performance':
            self.stats_performance()
        elif argv[1] == '-help':
            pipeline_help(self.logger)
        else:
            print('to-do options available')
            raise BadSettings('not analysis option choosen')

    def full_pipeline(self):
        """

        :return:
        """
        """
        if not self.catalogue(logger):
            raise Exception
        """
        if not self.sextractor():
            raise Exception
        if not self.scamp():
            raise Exception
        if not self.filt():
            raise Exception
        """
        if not self.check(confs, total_confs):
            raise Exception
        if not self.stats(total_confs):
            raise Exception
        """
        return True

    def catalog(self):
        """ todo - improve docstring
            todo - improve return
        :return:
        """
        # Gets an dictionary with analysis's preferences
        analysis_d, len_dicts = create_sextractor_dict(0, True)

        # Catalogue creation. Only created one time.
        CatalogCreation(self.logger, analysis_d)

        return True

    def sextractor(self):
        """ todo - improve docstring
            todo - improve return
        @return True: if everything goes alright.
        """
        mode = {'type': 'sextractor'}
        confs, total_confs = create_configurations(mode)

        for idx, conf_ in enumerate(confs):
            analysis_d, len_dicts = create_sextractor_dict(idx, False)
            # Just for tests reasons
            if len_dicts != total_confs:
                raise Exception
            # todo - implement a return!
            Sextractor(self.logger, analysis_d)

            self.logger.debug('perfoming analysis over a bunch of files')

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

        for mag in self.prfs_d['mags']:
            for idx_scmp, conf_scmp in enumerate(self.scamp_confs):
                # todo - check if everything is alright
                scmp_d, scmp_cf = scamp_f_name(idx_scmp)
                for idx_sex, conf_sex in enumerate(confs_sex):
                    analysis_d, len_dicts = create_sextractor_dict(idx_sex,
                                                                   False)
                    if not Scamp(self.logger, mag, scmp_d, scmp_cf,
                                 analysis_d):
                        raise Exception  # todo improve Exception

        return True

    def filt(self):
        """ todo - improve docstring
            todo - improve return
        :return:
        """

        confs = list(product(self.sex_confs, self.scamp_confs))

        for mag in self.prfs_d['mags']:
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
                                     args=(self.logger, mag,
                                           scmp_cf, sex_d,))
                    filt_j.append(filt_p)
                    filt_p.start()

                active_filt = list([j.is_alive() for j in filt_j])
                while True in active_filt:
                    active_filt = list([j.is_alive() for j in filt_j])
                    pass

        return True

    def check(self, confs, total_confs):
        """

        :param confs:
        :param total_confs:
        :return:
        """
        for mag in self.prfs_d['mags']:
            for idx_conf in range(0, total_confs, 2):
                try:
                    check_j = []
                    idx_proc = 0
                    while len(check_j) < 2:
                        idx_scmp = idx_conf + idx_proc
                        scmp_d, scmp_cf = scamp_f_name(idx_scmp)

                        check_p = Process(target=look_for_ssos,
                                          args=(self.logger, self.prfs_d,
                                                mag, scmp_d, scmp_cf, ))
                        check_j.append(check_p)
                        check_p.start()

                        idx_proc += 1

                    active_check = list([j.is_alive() for j in check_j])
                    while True in active_check:
                        active_check = list([j.is_alive() for j in check_j])
                        pass
                except IndexError:
                    print("finished")

        # Merge catalog files into a single one and cleans old ones
        for mag in self.prfs_d['mags']:
            for conf in confs:
                merge_ssos(self.logger, self.prfs_d, conf, mag)
                merge_sso_cat(self.logger, self.prfs_d, conf, mag)

        return True

    def stats(self, total_confs):
        """ todo - improve docstring
            todo - improve return

        :param total_confs:
        :return:
        """
        confs = list(product(self.sex_confs, self.scamp_confs))
        mag = '20-21'

        for idx in range(0, len(confs), self.prfs_d['cores_number']):
            filt_j = []
            while len(filt_j) < self.prfs_d['cores_number']:
                # bypassconfidence
                # Runs performance analysis.
                i = idx + len(filt_j)

                sex_cf = '{}_{}_{}_{}_{}'.format(confs[i][0][0], confs[i][0][2],
                                                 confs[i][0][2],
                                                 confs[i][0][1],
                                                 confs[i][0][3])
                scmp_cf = '{}_{}_{}_{}'.format(confs[i][1][0],
                                               confs[i][1][1],
                                               confs[i][1][2],
                                               confs[i][1][3])
                sex_d = {'deblending': confs[i][0][0],
                         'mincount': confs[i][0][1],
                         'threshold': confs[i][0][2],
                         'area': confs[i][0][3]}

                scmp_d = {'crossid_radius': confs[i][1][0],
                          'pixscale_maxerr': confs[i][1][1],
                          'posangle_maxerr': confs[i][1][2],
                          'position_maxerr': confs[i][1][3]}

                # f_name = 'stats_{}_{}.csv'.format(scmp_cf, mag)
                filt_p = Process(target=ExtractStats,
                                 args=(self.logger, mag, sex_cf, sex_d,
                                       scmp_cf, scmp_d, ))
                filt_j.append(filt_p)
                filt_p.start()

            active_filt = list([j.is_alive() for j in filt_j])
            while True in active_filt:
                active_filt = list([j.is_alive() for j in filt_j])
                pass

        return True

    def scamp_performance_ssos(self):
        """ Performs a complete pipeline to scamp output.

        @return True if everything goes alright.
        """
        stats_d = {}

        # idx = 0
        confs = list(product(self.sex_confs, self.scamp_confs))
        # mag = '20-21'

        for mag_ in self.prfs_d['mags']:
            for conf_ in confs:
                sex_cf = '{}_{}_{}_{}_{}'.format(conf_[0][0], conf_[0][2],
                                                 conf_[0][2], conf_[0][1],
                                                 conf_[0][3])
                scmp_cf = '{}_{}_{}_{}'.format(conf_[1][0], conf_[1][1],
                                               conf_[1][2], conf_[1][3])

                # ScampPerformance(self.logger, mag_, sex_cf, scmp_cf)
                TotalScampPerformance(self.logger, mag_, sex_cf, scmp_cf)

        return True

    def scamp_performance_stars(self):
        """ Performs a complete pipeline to scamp output only related to stars.
        No dictionary for statistics will be created.

        todo - improve docstring
        todo - improve return

        @return True if everything goes alright.
        """
        idx = 0
        stats_d = {}
        for mag in self.prfs_d['mags']:
            for idx_scmp, scmp_conf in enumerate(self.scamp_confs):
                for idx_sex, sex_conf in enumerate(self.sex_confs):
                    # Set an index.
                    # Scamp configuration.
                    # Creates a dict from a particular configuration.
                    scmp_d, scmp_cf = scamp_f_name(idx_scmp)

                    # Sextractor configuration.
                    conf = [sex_conf[0], sex_conf[2], sex_conf[2],
                            sex_conf[1], sex_conf[3]]
                    sex_cf = '{}_{}_{}_{}_{}'.format(conf[0], conf[1],
                                                     conf[2], conf[3],
                                                     conf[4])

                    # Runs performance analysis.
                    stats_d[idx] = ScampPerformanceStars(self.logger, mag,
                                                         sex_cf, scmp_cf)
                    idx += 1

        return True

    def pm_performance_ssos(self):
        """ Performs a complete pipeline to scamp output only related to ssos.
        No dictionary for statistics will be created.

        todo - improve docstring
        todo - improve return

        @return True if everything goes alright.
        """
        idx = 0
        stats_d = {}
        for mag in self.prfs_d['mags']:
            for idx_scmp, scmp_conf in enumerate(self.scamp_confs):
                for idx_sex, sex_conf in enumerate(self.sex_confs):
                    # Set an index.
                    # Scamp configuration.
                    # Creates a dict from a particular configuration.
                    scmp_d, scmp_cf = scamp_f_name(idx_scmp)

                    # Sextractor configuration.
                    conf = [sex_conf[0], sex_conf[2], sex_conf[2],
                            sex_conf[1], sex_conf[3]]
                    sex_cf = '{}_{}_{}_{}_{}'.format(conf[0], conf[1],
                                                     conf[2], conf[3],
                                                     conf[4])

                    # Runs performance analysis.
                    stats_d[idx] = SlowPMPerformanceSSOs(self.logger, mag,
                                                         sex_cf, scmp_cf)
                    idx += 1

        return True

    def fit_pm_mag(self):
        """ Performs a complete pipeline to scamp output only related to stars.
        No dictionary for statistics will be created.

        todo - improve docstring
        todo - improve return

        @return True if everything goes alright.
        """
        idx = 0
        stats_d = {}
        for mag in self.prfs_d['mags']:
            for idx_scmp, scmp_conf in enumerate(self.scamp_confs):
                for idx_sex, sex_conf in enumerate(self.sex_confs):
                    # Set an index.
                    # Scamp configuration.
                    # Creates a dict from a particular configuration.
                    scmp_d, scmp_cf = scamp_f_name(idx_scmp)

                    # Sextractor configuration.
                    conf = [sex_conf[0], sex_conf[2], sex_conf[2],
                            sex_conf[1], sex_conf[3]]
                    sex_cf = '{}_{}_{}_{}_{}'.format(conf[0], conf[1],
                                                     conf[2], conf[3],
                                                     conf[4])

                    # Runs performance analysis.
                    stats_d[idx] = FitPMMagAgainstSizes(self.logger, mag,
                                                        sex_cf, scmp_cf)
                    idx += 1

        return True

    def scamp_performance_fitting(self):
        """ Performs a complete pipeline to scamp output only related to stars.
        No dictionary for statistics will be created.

        todo - improve docstring
        todo - improve return

        @return True if everything goes alright.
        """
        idx = 0
        stats_d = {}
        for mag in self.prfs_d['mags']:
            for idx_scmp, scmp_conf in enumerate(self.scamp_confs):
                for idx_sex, sex_conf in enumerate(self.sex_confs):
                    # Set an index.
                    # Scamp configuration.
                    # Creates a dict from a particular configuration.
                    scmp_d, scmp_cf = scamp_f_name(idx_scmp)

                    # Sextractor configuration.
                    conf = [sex_conf[0], sex_conf[2], sex_conf[2],
                            sex_conf[1], sex_conf[3]]
                    sex_cf = '{}_{}_{}_{}_{}'.format(conf[0], conf[1],
                                                     conf[2], conf[3],
                                                     conf[4])

                    # Runs performance analysis.
                    stats_d[idx] = ScampPerformanceStars(self.logger, mag,
                                                         sex_cf, scmp_cf)
                    idx += 1

        return True

    def sextractor_size(self):
        """ plotea el tamaño de los objetos atendiendo a su naturaleza

        todo - improve docstring
        todo - improve return - save an entire file with all valuable data

        :return:
        """
        # Sextractor configurations.
        mode = {'type': 'sextractor'}
        sex_confs, sex_confs_n = create_configurations(mode)

        for mag in self.prfs_d['mags']:
            for idx_sex, sex_conf in enumerate(sex_confs):
                # Set an index.
                # Sextractor configuration.
                conf = [sex_conf[0], sex_conf[2], sex_conf[2],
                        sex_conf[1], sex_conf[3]]
                sex_cf = '{}_{}_{}_{}_{}'.format(conf[0], conf[1],
                                                 conf[2], conf[3],
                                                 conf[4])

                # Runs performance analysis.
                SextractorSizes(self.logger, mag, sex_cf)

        return True

    def error_performance(self):
        """ todo - improve docstring
            todo - improve return

        :return:
        """
        stats_d = {}
        idx = 0

        # Sextractor configurations.
        mode = {'type': 'sextractor'}
        sex_confs, sex_confs_n = create_configurations(mode)

        for mag in self.prfs_d['mags']:
            for idx_sex, sex_conf in enumerate(sex_confs):
                # Set an index.
                # Sextractor configuration.
                conf = [sex_conf[0], sex_conf[2], sex_conf[2],
                        sex_conf[1], sex_conf[3]]
                sex_cf = '{}_{}_{}_{}_{}'.format(conf[0], conf[1],
                                                 conf[2], conf[3],
                                                 conf[4])

                # Runs performance analysis.
                stats_d[idx] = SextractorPerformance().error(self.logger,
                                                             mag, sex_cf)
                idx += 1

        return True

    def pm_performance(self):
        """ input pm vs output pm
        todo - check syntax

        :return:
        """
        idx = 0
        stats_d = {}
        for mag in self.prfs_d['mags']:
            for idx_scmp, scmp_conf in enumerate(self.scamp_confs):
                for idx_sex, sex_conf in enumerate(self.sex_confs):
                    # Set an index.
                    # Scamp configuration.
                    # Creates a dict from a particular configuration.
                    scmp_d, scmp_cf = scamp_f_name(idx_scmp)

                    # Sextractor configuration.
                    conf = [sex_conf[0], sex_conf[2], sex_conf[2],
                            sex_conf[1], sex_conf[3]]
                    sex_cf = '{}_{}_{}_{}_{}'.format(conf[0], conf[1],
                                                     conf[2], conf[3],
                                                     conf[4])

                    # bypassconfidence
                    self.prfs_d['confidences'] = [1]
                    # Runs performance analysis.
                    for confidence_ in self.prfs_d['confidences']:
                        idx += 1
                        stats_d = PMPerformance().check(self.logger, mag,
                                                        scmp_cf, sex_cf,
                                                        confidence_)

        stats_df = DataFrame(stats_d, columns=['pm', 'mean', 'median',
                                               'std', 'max', 'min',
                                               'detected', 'total'])
        stats_df = stats_df.sort_values(stats_df.columns[0], ascending=True)
        stats_df = stats_df.reset_index(drop=True)
        stats_df.to_csv('stats.csv')

        return True

    def stats_performance(self):
        """ Performs a complete pipeline to scamp output.
        todo -

        @return True if everything goes alright.
        """
        stats_d = {}
        idx = 0

        # Sextractor configurations.
        mode = {'type': 'sextractor'}
        sex_confs, sex_confs_n = create_configurations(mode)

        for mag in self.prfs_d['mags']:
            for idx_sex, sex_conf in enumerate(sex_confs):
                # Set an index.
                # Sextractor configuration.
                conf = [sex_conf[0], sex_conf[2], sex_conf[2],
                        sex_conf[1], sex_conf[3]]
                sex_cf = '{}_{}_{}_{}_{}'.format(conf[0], conf[1],
                                                 conf[2], conf[3],
                                                 conf[4])

                # Runs performance analysis.
                stats_d[idx] = StatsPerformance().error(self.logger, mag,
                                                        sex_cf)
                idx += 1

        tmp_d = {'conf': [], 'm_a_stars': [], 'm_b_stars': [], 'm_a_gals': [],
                 'm_b_gals': [], 'm_a_ssos': [], 'm_b_ssos': [],
                 's_a_stars': [], 's_b_stars': [], 's_a_gals': [],
                 's_b_gals': [], 's_a_ssos': [], 's_b_ssos': []}
        for conf_key in stats_d.keys():
            print('keys {}'.format(stats_d.keys()))
            """
            for value_key in stats_d[conf_key].keys():
                for value in stats_d[conf_key][value_key]:
                    tmp_d[value_key].append(value)
            """
        stats_df = DataFrame(tmp_d)
        stats_df.to_csv('std.csv')

        return True


if __name__ == '__main__':
    check_process = Check()
