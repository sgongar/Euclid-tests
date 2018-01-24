#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance
Mide las detecciones de cada tipo

Versions:
- 0.1

Todo:
    * Improve log messages
    * Use same naming convention for all variables

"""
from numpy import array
import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv

from cats_management import check_source, get_input_dicts
from misc import extract_settings, get_dither, get_norm_speed
from scipy.odr import *


def function_returner(p, input_data):
    input_data = array(input_data)
    x = input_data[0]
    z = input_data[1]

    return p[0] * x + p[1] * z + p[2]


def redo_stars_d():
    """ Creates a dictionary

    :return: tmp_d
    """
    stars_d = {'catalog_n': [], 'alpha_j2000': [], 'delta_j2000': [],
               'stars_mag': [], 'stars_mag_err': [], 'stars_pm': [],
               'stars_pm_err': [], 'stars_a': [], 'stars_err_a': [],
               'stars_b': [], 'stars_err_b': [], 'stars_elongation': [],
               'stars_ellipticity': []}

    return stars_d


def redo_galaxies_d():
    """ Creates a dictionary

    :return: galaxies_d
    """
    galaxies_d = {'catalog_n': [], 'alpha_j2000': [], 'delta_j2000': [],
                  'galaxies_mag': [], 'galaxies_mag_err': [],
                  'galaxies_pm': [], 'galaxies_pm_err': [], 'galaxies_a': [],
                  'galaxies_err_a': [], 'galaxies_b': [], 'galaxies_err_b': [],
                  'galaxies_elongation': [], 'galaxies_ellipticity': []}

    return galaxies_d


def redo_ssos_d():
    """

    :return: ssos_d
    """
    ssos_d = {'catalog_n': [], 'alpha_j2000': [], 'delta_j2000': [],
              'ssos_mag': [], 'ssos_mag_err': [], 'ssos_pm': [],
              'ssos_pm_err': [], 'ssos_a': [], 'ssos_err_a': [], 'ssos_b': [],
              'ssos_err_b': [], 'ssos_elongation': [], 'ssos_ellipticity': []}

    return ssos_d


def read_fit_d():
    """ Reads dictionary from

    :return: fit_d
    """
    fit_d = read_csv('fit_df.csv')

    return fit_d


class FitPMMagAgainstSizes:

    def __init__(self, logger, mag, sex_cf, scmp_cf):
        """

        :param logger:
        :param mag:
        :param sex_cf:
        :param scmp_cf:
        """
        self.logger = logger
        self.mag = mag
        self.sex_cf = sex_cf
        self.scmp_cf = scmp_cf
        self.prfs_d = extract_settings()

        self.stars_d = redo_stars_d()  # Creates a dictionary for statistics
        self.galaxies_d = redo_galaxies_d()
        self.ssos_d = redo_ssos_d()
        self.get()  # Gets data from catalogs
        # self.read_df()
        # self.analyse()  # Fits data to obtained data

    def read_df(self):
        """ TODO get-dir

        :return: stars_d, galaxies_d, ssos_d
        """
        stars_df = read_csv('stars_df.csv', index_col=0)
        galaxies_df = read_csv('galaxies_df.csv', index_col=0)
        ssos_df = read_csv('ssos_df.csv', index_col=0)

        self.stars_d = stars_df.to_dict()
        self.galaxies_d = galaxies_df.to_dict()
        self.ssos_d = ssos_df.to_dict()

    def save_df(self):
        """

        :return: True
        """
        stars_df = DataFrame(self.stars_d)
        stars_df.to_csv('stars_df.csv')

        galaxies_df = DataFrame(self.galaxies_d)
        galaxies_df.to_csv('galaxies_df.csv')

        ssos_df = DataFrame(self.ssos_d)
        ssos_df.to_csv('ssos_df.csv')

        return True

    def get(self):
        """

        :return:
        """
        # Creates an input dictionary with all input sources
        self.logger.debug('scamp {} and sextractor {}'.format(self.scmp_cf,
                                                              self.sex_cf))

        # Gets dataframes for each type of source
        i_df = get_input_dicts(self.mag)

        # Gets the name of filtered file
        filter_n = 'filt_{}_{}_3.csv'.format(self.scmp_cf, self.mag)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filter_n)
        # Opens filtered file
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)
        # Gets unique sources from filtered file
        o_uniq_sources = list(set(o_cat['SOURCE_NUMBER'].tolist()))

        # Loops over unique sources of filtered file
        # TODO splits over different threads
        # n_theads = 4
        # lista_n = [o_uniq_sources[i:i + n] for i in
        # xrange(0, len(o_uniq_sources), n)]
        a_image = 0.0
        b_image = 0.0
        err_a_image = 0.0
        err_b_image = 0.0
        mag = 0.0
        mag_err = 0.0
        o_pm = 0.0
        o_pm_err = 0.0
        for idx, source_ in enumerate(o_uniq_sources):
            # Check if actual source lies in 20-21 magnitude gap
            flags = []
            o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                a_image = row.A_IMAGE
                b_image = row.B_IMAGE
                err_a_image = row.ERRA_IMAGE
                err_b_image = row.ERRB_IMAGE
                mag = row.MAG
                mag_err = row.MAGERR
                catalog_n = row.CATALOG_NUMBER
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000
                o_pm = row.PM
                o_pm_err = row.PMERR
                elongation = row.ELONGATION
                ellipticity = row.ELLIPTICITY

                # look for object type!
                dither = get_dither(catalog_n)

                i_stars_df = i_df['stars']
                i_stars_d_df = i_stars_df[
                    i_stars_df['dither_values'].isin([dither])]

                i_galaxies_df = i_df['galaxies']
                i_galaxies_d_df = i_galaxies_df[
                    i_galaxies_df['dither_values'].isin([dither])]

                i_ssos_df = i_df['ssos']
                i_ssos_d_df = i_ssos_df[
                    i_ssos_df['dither_values'].isin([dither])]

                out_df = check_source(catalog_n, i_stars_d_df,
                                      o_alpha, o_delta)
                if out_df.empty is not True:
                    flags.append('star')

                out_df = check_source(catalog_n, i_galaxies_d_df,
                                      o_alpha, o_delta)
                if out_df.empty is not True:
                    flags.append('galaxy')

                out_df = check_source(catalog_n, i_ssos_d_df,
                                      o_alpha, o_delta)
                if out_df.empty is not True:
                    flags.append('sso')

            if len(set(flags)) == 1:
                if flags[0] == 'star':
                    self.stars_d['catalog_n'].append(catalog_n)
                    self.stars_d['alpha_j2000'].append(o_alpha)
                    self.stars_d['delta_j2000'].append(o_delta)
                    self.stars_d['stars_mag'].append(mag)
                    self.stars_d['stars_mag_err'].append(mag_err)
                    # o_pm_norm = get_norm_speed(o_pm)
                    # self.stars_d['stars_pm'].append(o_pm_norm)
                    self.stars_d['stars_pm'].append(o_pm)
                    self.stars_d['stars_pm_err'].append(o_pm_err)
                    self.stars_d['stars_a'].append(a_image)
                    self.stars_d['stars_err_a'].append(err_a_image)
                    self.stars_d['stars_b'].append(b_image)
                    self.stars_d['stars_err_b'].append(err_b_image)
                    self.stars_d['stars_elongation'].append(elongation)
                    self.stars_d['stars_ellipticity'].append(ellipticity)
                elif flags[0] == 'galaxy':
                    self.galaxies_d['catalog_n'].append(catalog_n)
                    self.galaxies_d['alpha_j2000'].append(o_alpha)
                    self.galaxies_d['delta_j2000'].append(o_delta)
                    self.galaxies_d['galaxies_mag'].append(mag)
                    self.galaxies_d['galaxies_mag_err'].append(mag_err)
                    # o_pm_norm = get_norm_speed(o_pm)
                    # self.galaxies_d['galaxies_pm'].append(o_pm_norm)
                    self.galaxies_d['galaxies_pm'].append(o_pm)
                    self.galaxies_d['galaxies_pm_err'].append(o_pm_err)
                    self.galaxies_d['galaxies_a'].append(a_image)
                    self.galaxies_d['galaxies_err_a'].append(err_a_image)
                    self.galaxies_d['galaxies_b'].append(b_image)
                    self.galaxies_d['galaxies_err_b'].append(err_b_image)
                    self.galaxies_d['galaxies_elongation'].append(elongation)
                    self.galaxies_d['galaxies_ellipticity'].append(ellipticity)
                elif flags[0] == 'sso':
                    self.ssos_d['catalog_n'].append(catalog_n)
                    self.ssos_d['alpha_j2000'].append(o_alpha)
                    self.ssos_d['delta_j2000'].append(o_delta)
                    self.ssos_d['ssos_mag'].append(mag)
                    self.ssos_d['ssos_mag_err'].append(mag_err)
                    # o_pm_norm = get_norm_speed(o_pm)
                    # self.ssos_d['ssos_pm'].append(o_pm_norm)
                    self.ssos_d['ssos_pm'].append(o_pm)
                    self.ssos_d['ssos_pm_err'].append(o_pm_err)
                    self.ssos_d['ssos_a'].append(a_image)
                    self.ssos_d['ssos_err_a'].append(err_a_image)
                    self.ssos_d['ssos_b'].append(b_image)
                    self.ssos_d['ssos_err_b'].append(err_b_image)
                    self.ssos_d['ssos_elongation'].append(elongation)
                    self.ssos_d['ssos_ellipticity'].append(ellipticity)
                else:
                    raise Exception

        if not self.save_df():
            raise Exception

    def analyse(self):
        """

        :return:
        """
        stars_mag = []
        for key_ in self.stars_d['stars_mag'].keys():
            stars_mag.append(self.stars_d['stars_mag'][key_])
        stars_mag_err = []
        for key_ in self.stars_d['stars_mag_err'].keys():
            stars_mag_err.append(self.stars_d['stars_mag_err'][key_])

        stars_pm = []
        for key_ in self.stars_d['stars_pm'].keys():
            stars_pm.append(self.stars_d['stars_pm'][key_])
        stars_pm_err = []
        for key_ in self.stars_d['stars_pm_err'].keys():
            stars_pm_err.append(self.stars_d['stars_pm_err'][key_])

        stars_a = []
        for key_ in self.stars_d['stars_a'].keys():
            stars_a.append(self.stars_d['stars_a'][key_])
        stars_err_a = []
        for key_ in self.stars_d['stars_err_a'].keys():
            stars_err_a.append(self.stars_d['stars_err_a'][key_])

        x_data = array(stars_mag)
        x_error = array(stars_mag_err)
        z_data = array(stars_pm)
        z_error = array(stars_pm_err)
        y_data = array(stars_a)
        y_error = array(stars_err_a)

        regression_model = Model(function_returner)
        fit_data = RealData([x_data, z_data], y_data,
                            sx=[x_error, z_error], sy=[y_error])
        odr_model = ODR(fit_data, regression_model, beta0=[0.04, -0.02, 1.75])
        odr_model.set_job(fit_type=0)
        out = odr_model.run()
        # result = out.beta

        print(out.pprint())

        # Plot the data
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(x_data, y_data, z_data)

        plt.show()

        """
        import scipy.optimize as optimize
        import numpy as np

        A = np.array([(19, 20, 24), (10, 40, 28), (10, 50, 31)])

        def func(data, a, b):
            return data[:, 0] * data[:, 1] * a + b

        guess = (1, 1)
        params, pcov = optimize.curve_fit(func, A[:, :2], A[:, 2], guess)
        print(params)
        """