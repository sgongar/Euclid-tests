#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.1

Todo:
    * Improve log messages

"""
from astropy.io import fits
from astropy.table import Table
from numpy import array, std, sqrt, linspace
from pandas import concat, read_csv, DataFrame, Series
import statsmodels.api as sm
from scipy.odr import Model, ODR, RealData
from scipy.stats import linregress

from misc import all_same, extract_settings
from misc import create_folder, speeds_range
from plots import PlotConfidence, PlotError
from plot_fitting import PlotFitting
from regions import Create_regions


def check_source(catalog_n, o_cat, i_alpha, i_delta):
    """

    :param catalog_n:
    :param o_cat:
    :param i_alpha:
    :param i_delta:
    :return:
    """
    # tolerance = 0.0001  # 1.8 arcsecond
    tolerance = 0.0005  # 1.8 arcsecond

    o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
    o_df = o_df[o_df['ALPHA_J2000'] + tolerance > i_alpha]
    o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - tolerance]
    o_df = o_df[o_df['DELTA_J2000'] + tolerance > i_delta]
    o_df = o_df[i_delta > o_df['DELTA_J2000'] - tolerance]

    return o_df


class PMPerformanceSSOs:

    def __init__(self, logger, mag, scmp_cf, sex_cf, confidence_):
        """

        """
        self.logger = logger
        self.bypassfilter = True
        self.checkerror = True
        self.prfs_d = extract_settings()

        self.mag = mag
        self.scmp_cf = scmp_cf
        self.sex_cf = sex_cf
        self.confidence_ = confidence_

        self.save = False

        self.check_pm_distribution()

    def creates_input_dict(self):
        """ Creates an input dictionary. Each key contains SSOs' information
        for each dither.

        :return: input_dict
        """
        input_dict = {}
        # Loops over the four dithers
        for dither in range(1, 5, 1):
            cat_location = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'],
                                                   self.mag)
            cat_name = '{}/Cat_20-21_d{}'.format(cat_location, dither)
            input_dict[dither] = '{}.dat'.format(cat_name)
        input_dict = Create_regions(input_dict).check_ssos(self.mag, True)

        return input_dict

    def creates_input_df(self, input_dict):
        """ Creates an input dataframe from an input dictionary.

        :return: input dataframe
        """
        occurrences = 3

        input_list = []
        for key_ in input_dict.keys():
            input_list.append(input_dict[key_])

        input_df = concat(input_list, axis=0)
        # Look for >= 3 coincidences
        input_df = concat(g for _, g in input_df.groupby('source')
                          if len(g) >= occurrences)
        input_df = input_df.reset_index(drop=True)

        if self.save:
            input_df.to_csv('inputs.csv')

        return input_df

    def check_pm_distribution(self):
        """

        :return:
        """
        # Creates an input dictionary with all input sources
        self.logger.debug('Scamp configuration: {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration: {}'.format(self.sex_cf))

        input_dict = self.creates_input_dict()
        input_df = self.creates_input_df(input_dict)

        # Open particular file!
        filter_n = 'filt_{}_{}_2.csv'.format(self.scmp_cf, self.mag)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filter_n)
        self.logger.debug('opens filtered catalog {}'.format(filter_n))
        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        # stats_d, out_d = create_dict(self.scmp_cf, self.sex_cf,
        #                              self.confidence_)

        # Gets unique sources from input data
        unique_sources = list(set(input_df['source'].tolist()))
        """
        unique_sources = []
        for unique_source_ in unique_sources_:
            if int(unique_source_) > 170:
                unique_sources.append(unique_source_)
        """
        # Loops over input data
        for idx_source, source_ in enumerate(unique_sources):
            # Gets associated data in input catalog
            cat_df = i_df[i_df['source'].isin([source_])]
            # Creates lists for each source
            tmp_d = redo_tmp_d()
            err_d = redo_err_d()
            i_pm = 0.0  # Raises a warning if variable pm it's not created
            # Iterate over each detection of each source
            for i, row in enumerate(cat_df.itertuples(), 1):
                input_source = row.source
                # ccd_ = row.CCD
                # print('CCD {}'.format(ccd_))
                # dither_ = row.dither_values
                # print('dither {}'.format(dither_))

                catalog_n = row.catalog
                i_pm = row.pm_values
                i_pm_alpha = row.pm_alpha  # this is wrong! fixme
                i_pm_delta = row.pm_delta
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # print('i_alpha {}'.format(i_alpha))
                # print('i_delta {}'.format(i_delta))

                # Checks if there is a source closed to input one
                # o_cat contains data from output (filtered) catalog
                o_df = check_source(catalog_n, o_cat, i_alpha, i_delta)
                print('o_df {}'.format(o_df))

                # If there is one saves data from input data
                if o_df.empty is not True and o_df['PM'].size == 1:
                    pm_mask = self.pm_filter(o_df, i_pm, self.confidence_)
                    if pm_mask:
                        if self.checkerror:
                            scmp_source = float(o_df['SOURCE_NUMBER'].iloc[0])
                            # print('source {}'.format(scmp_source))
                            o_alpha = float(o_df['ALPHA_J2000'].iloc[0])
                            o_delta = float(o_df['DELTA_J2000'].iloc[0])
                            o_pm = float(o_df['PM'].iloc[0])
                            o_pm_alpha = float(o_df['PMALPHA'].iloc[0])
                            o_pm_delta = float(o_df['PMDELTA'].iloc[0])
                            # print('o_alpha {}'.format(o_alpha))
                            # print('o_delta {}'.format(o_delta))
                            err_d = self.populate_err_dict(catalog_n, i_pm,
                                                           i_pm_alpha, i_alpha,
                                                           i_pm_delta, i_delta,
                                                           o_pm, o_pm_alpha,
                                                           o_pm_delta,
                                                           cat_df, err_d,
                                                           scmp_source,
                                                           o_alpha, o_delta)
                        else:
                            # Populate a dictionary of temporal data
                            # for data plotting
                            tmp_d = self.populate_tmp_dict(catalog_n, i_pm,
                                                           i_pm_alpha, i_alpha,
                                                           i_pm_delta, i_delta,
                                                           o_df, tmp_d)
                else:
                    tmp_d['boolean_l'].append('False')
                    # cat_df is an input catalog!
                    # should look to this area, where is in that area
                    """
                    scmp_source = float(o_df['SOURCE_NUMBER'].iloc[0])
                    o_alpha = float(o_df['ALPHA_J2000'].iloc[0])
                    o_delta = float(o_df['DELTA_J2000'].iloc[0])
                    o_pm = float(o_df['PM'].iloc[0])
                    o_pm_alpha = float(o_df['PMALPHA'].iloc[0])
                    o_pm_delta = float(o_df['PMDELTA'].iloc[0])
                    err_d = self.populate_err_dict(catalog_n, i_pm, i_pm_alpha,
                                                   i_alpha, i_pm_delta,
                                                   i_delta, o_pm, o_pm_alpha,
                                                   o_pm_delta, cat_df, err_d,
                                                   scmp_source, o_alpha,
                                                   o_delta)
                    """
                    scmp_source = False
                    o_alpha = False
                    o_delta = False
                    o_pm = False
                    o_pm_alpha = False
                    o_pm_delta = False
                    err_d = self.populate_err_dict(catalog_n, i_pm, i_pm_alpha,
                                                   i_alpha, i_pm_delta,
                                                   i_delta, o_pm, o_pm_alpha,
                                                   o_pm_delta, cat_df, err_d,
                                                   scmp_source, o_alpha,
                                                   o_delta)

            logger.debug('input source is {}'.format(input_source))

            # Total number
            idx = stats_d['PM'].index(i_pm)
            stats_d['total'][idx] += 1

            order_mask = check_cat_order(tmp_d['catalog'])

            if not order_mask:
                raise Exception

            # plot errors!
            if False in err_d['scmp_source']:
                false_times = err_d['scmp_source'].count(False)
            else:
                false_times = 0

            if false_times > 1:
                false_flag = False
            else:
                false_flag = True

            if len(err_d['scmp_source']) is not 0 and false_flag:
                flag_detection, sources_number = all_same(err_d['scmp_source'])
                if flag_detection is False:
                    print(err_d['scmp_source'])

                if flag_detection and sources_number >= 3:
                    if self.checkerror:
                        logger.debug('right detection')
                        # Output folder creation
                        plots_dir = self.prfs_d['plots_dir']
                        output_path = '{}/err/{}/{}/{}'.format(plots_dir,
                                                               self.scmp_cf,
                                                               self.sex_cf,
                                                               err_d['i_pm'][0])
                        create_folder(logger, output_path)

                        fits_ = []  # todo - move to dictionary
                        for idx, catalog_ in enumerate(err_d['catalog']):
                            ccd = self.get_ccd(catalog_)
                            fits_.append(ccd)

                        fitted_d = self.confidence(logger,
                                                   err_d['scmp_source'][1],
                                                   self.scmp_cf, self.sex_cf)

                        ok = 'yes'
                        self.plot_err(output_path, err_d, fits_, ok, fitted_d)

                    """
                    idx = stats_d['PM'].index(i_pm)
                    stats_d['right'][idx] += 1
                    fitted_d = self.confidence(tmp_d['source'][0],
                                               self.scmp_cf, self.sex_cf)

                    pm = float(tmp_d['i_pm'][0])

                    # Output folder creation
                    plots_dir = self.prfs_d['plots_dir']
                    self.output_path = '{}/mv/{}/{}/{}'.format(plots_dir,
                                                               self.scmp_cf,
                                                               self.sex_cf,
                                                               pm)
                    create_folder(logger, self.output_path)

                    fits = []  # todo - move to dictionary
                    for idx, catalog_ in enumerate(tmp_d['catalog']):
                        ccd = self.get_ccd(catalog_)
                        fits.append(ccd)

                    self.plot(tmp_d, pm, fitted_d, fits)
                    """
                else:
                    if self.checkerror:
                        logger.debug(
                            'wrong detection - not enough sources number')
                        # Output folder creation
                        plots_dir = self.prfs_d['plots_dir']
                        output_path = '{}/err/{}/{}/{}'.format(plots_dir,
                                                               self.scmp_cf,
                                                               self.sex_cf,
                                                               err_d['i_pm'][0])
                        create_folder(logger, output_path)

                        fits_ = []  # todo - move to dictionary
                        for idx, catalog_ in enumerate(err_d['catalog']):
                            ccd = self.get_ccd(catalog_)
                            fits_.append(ccd)
                        """
                        fitted_d = self.confidence(err_d['scmp_source'][0],
                                                   self.scmp_cf, self.sex_cf)
                        """
                        fitted_d = {'ra': 'empty', 'dec': 'empty'}
                        ok = '2'
                        self.plot_err(output_path, err_d, fits_, ok, fitted_d)
                    else:
                        pass  # nothing to do
            else:
                if self.checkerror:
                    logger.debug('wrong detection - too many false detections')
                    # Output folder creation
                    plots_dir = self.prfs_d['plots_dir']
                    output_path = '{}/err/{}/{}/{}'.format(plots_dir,
                                                           self.scmp_cf,
                                                           self.sex_cf,
                                                           err_d['i_pm'][0])
                    create_folder(logger, output_path)

                    fits_ = []  # todo - move to dictionary
                    for idx, catalog_ in enumerate(err_d['catalog']):
                        ccd = self.get_ccd(catalog_)
                        fits_.append(ccd)

                        """
                        fitted_d = self.confidence(err_d['scmp_source'][0],
                                                   self.scmp_cf, self.sex_cf)
                        """
                        fitted_d = {'ra': 'empty', 'dec': 'empty'}
                    ok = '3'
                    self.plot_err(output_path, err_d, fits_, ok, fitted_d)
                else:
                    pass  # nothing to do

        return stats_d

    def pm_filter(self, o_df, pm, confidence_):
        """

        :param o_df:
        :param pm:
        :param confidence_:
        :return:
        """
        pm_ranges = speeds_range(self.prfs_d, confidence_)
        pm_range = pm_ranges[pm]

        if pm_range[0] < float(o_df['PM']) < pm_range[1]:
            return True
        elif self.bypassfilter:
            return True
        else:
            return False

    def plot(self, tmp_d, pm, fitted_d, fits):
        """

        :param tmp_d:
        :param pm:
        :param fitted_d:
        :param fits:
        :return:
        """
        # Set True to plot input and output data
        both = False

        if both:
            mode = 'io'
            # d_ = {'i_alpha': tmp_d['i_alpha'], 'i_delta': tmp_d['i_delta'],
            #       'i_pm_alpha': tmp_d['i_pm_alpha'],
            #       'i_pm_delta': tmp_d['i_pm_delta'],
            #       'o_alpha': tmp_d['o_alpha'], 'o_delta': tmp_d['o_delta'],
            #       'o_pm_alpha': tmp_d['o_pm_alpha'],
            #       'o_pm_delta': tmp_d['o_pm_delta'],
            #       'o_pm_alpha_err': tmp_d['o_pm_alpha_err'],
            #       'o_pm_delta_err': tmp_d['o_pm_delta_err'],
            #       'error_a': tmp_d['error_a'], 'error_b': tmp_d['error_b'],
            #       'epoch': tmp_d['epoch'], 'i_pm': tmp_d['i_pm']}
            # plot = PlotBothConfidence(self.output_path, tmp_d['source'][0],
            #                           pm, mode, fitted_d, d_)
            # if not plot:
            #     raise Exception
        else:
            mode = 'o'
            dict_ = {'sex_cf': tmp_d['sex_cf'], 'scmp_cf': tmp_d['scmp_cf'],
                     'alpha': tmp_d['o_alpha'], 'delta': tmp_d['o_delta'],
                     'error_a': tmp_d['error_a'], 'error_b': tmp_d['error_b'],
                     'i_pm': tmp_d['i_pm'], 'i_pm_alpha': tmp_d['i_pm_alpha'],
                     'i_pm_delta': tmp_d['i_pm_delta'],
                     'o_pm_alpha': tmp_d['o_pm_alpha'],
                     'o_pm_delta': tmp_d['o_pm_delta'],
                     'o_pm_alpha_err': tmp_d['o_pm_alpha_err'],
                     'o_pm_delta_err': tmp_d['o_pm_delta_err'],
                     'epoch': tmp_d['epoch']}
            plot = PlotConfidence(self.output_path, tmp_d['source'][0], pm,
                                  mode, fitted_d, dict_, fits, self.mag)

            if not plot:
                raise Exception

    def plot_err(self, output_dir, err_d, fits_, ok, fitted_d):
        """

        :param output_dir:
        :param err_d:
        :param fits_:
        :param ok:
        :param fitted_d:
        :return:
        """
        PlotError(output_dir, err_d, fits_, self.mag, ok, fitted_d)

    def confidence(self, logger, source, scmp_cf, sex_cf):
        """

        :param source:
        :param scmp_cf:
        :param sex_cf:
        :return:
        """
        catalogs_dir = self.prfs_d['catalogs_dir']
        configuration_dir = '/{}/{}/{}'.format(self.mag, sex_cf, scmp_cf)
        cat_name = '/full_{}_20-21_1.cat'.format(scmp_cf)
        hdu_list = fits.open('{}{}{}'.format(catalogs_dir,
                                             configuration_dir, cat_name))
        db = Table(hdu_list[2].data).to_pandas()

        ra = db.loc[db['SOURCE_NUMBER'] == source, 'ALPHA_J2000'].tolist()
        dec = db.loc[db['SOURCE_NUMBER'] == source, 'DELTA_J2000'].tolist()
        epoch = db.loc[db['SOURCE_NUMBER'] == source, 'EPOCH'].tolist()

        # print('ra {}'.format(ra))
        # print('dec {}'.format(dec))

        # edims = []
        epochs = []
        dimensions = []
        fitted_d = {}

        for dimension in [ra, dec]:
            x = array(epoch)
            epochs.append(x)  # epochs list
            y = array(dimension)
            dimensions.append(y)  # dimensions list
            if dimension == ra:
                sigma = db.loc[db['SOURCE_NUMBER'] == source,
                               'ERRA_WORLD'].tolist()
                dimension_ = 'ra'
                print('sigma ra {}'.format(sigma))
                edim = array([1 / var for var in sigma])
                # edims.append(edim)
                x = sm.add_constant(x)
                model = sm.WLS(y, x, weigths=edim)
                fitted = model.fit()
                print(fitted.summary())
                chi_squared = fitted.rsquared
                fitted_d[dimension_] = chi_squared
            elif dimension == dec:
                sigma = db.loc[db['SOURCE_NUMBER'] == source,
                               'ERRB_WORLD'].tolist()

                print('sigma dec {}'.format(sigma))
                dimension_ = 'dec'
                edim = array([1 / var for var in sigma])
                # edims.append(edim)
                x = sm.add_constant(x)
                model = sm.WLS(y, x, weigths=edim)
                fitted = model.fit()
                print(fitted.summary())
                chi_squared = fitted.rsquared
                fitted_d[dimension_] = chi_squared

        x = array(ra)
        y = array(dec)
        sigma = db.loc[db['SOURCE_NUMBER'] == source, 'ERRA_WORLD'].tolist()
        edim = array([1 / var for var in sigma])
        x = sm.add_constant(x)
        model = sm.WLS(y, x, weigths=edim)
        fitted = model.fit()
        print('new {}'.format(fitted.summary()))
        logger.debug('new chi-squared is {}'.format(fitted.rsquared))

        return fitted_d

    def create_regions(self, source, o_cat, tmp_d):
        """

        :param source:
        :param o_cat:
        :param tmp_d:
        :return:
        """
        # En azul! Valores totales
        # Gets full list of alpha/delta coordinates
        alpha_l = tmp_d['i_alpha'] + tmp_d['o_alpha']
        delta_l = tmp_d['i_delta'] + tmp_d['o_delta']

        d_limits = {'max_alpha': float(max(alpha_l)),
                    'min_alpha': float(min(alpha_l)),
                    'max_delta': float(max(delta_l)),
                    'min_delta': float(min(delta_l))}

        margin = 0.001
        o_df = cut_catalog(o_cat, margin, d_limits)

        # o_df.to_csv('{}.reg'.format(source), index=False,
        #             header=False, sep=" ")
        # Tests reasons
        f_regions_filename = 'f_{}.reg'.format(source)
        o_df.to_csv(f_regions_filename, index=False, sep=" ")

        # En rojo! Input values from Luca's
        alpha_list = Series(tmp_d['i_alpha'], name='ALPHA_J2000')
        delta_list = Series(tmp_d['i_delta'], name='DELTA_J2000')

        i_df = concat([alpha_list, delta_list], axis=1)
        # i_df.to_csv('i_{}.reg'.format(source), index=False,
        #             header=False, sep=" ")
        i_regions_filename = 'i_{}.reg'.format(source)
        i_df.to_csv(i_regions_filename, index=False, sep=" ")

        d_regions_filename = {'f_filename': f_regions_filename,
                              'i_filename': i_regions_filename}

        return d_regions_filename

    def create_regions_files(self, i_df):
        """

        :param i_df:
        :return:
        """
        # Saves a csv file populated by input sources
        # and a regions file of them.
        i_df.to_csv('inputs.csv')
        alpha_df = i_df['alpha_j2000']
        delta_df = i_df['delta_j2000']

        df = concat([alpha_df, delta_df], axis=1)
        input_regs = '{}/full.reg'.format(self.prfs_d['dithers_out'])
        df.to_csv(input_regs, index=False, header=False, sep=" ")

        for dither in range(1, 5, 1):
            i_df_dither = i_df[i_df['dither_values'].isin([dither])]
            alpha_df = i_df_dither['alpha_j2000']
            delta_df = i_df_dither['delta_j2000']

            df = concat([alpha_df, delta_df], axis=1)
            input_regs = '{}/dither_{}.reg'.format(self.prfs_d['dithers_out'],
                                                   dither)
            df.to_csv(input_regs, index=False, header=False, sep=" ")

    def populate_tmp_dict(self, catalog_n, i_pm, i_pm_alpha, i_pm_delta,
                          i_alpha, i_delta, o_df, tmp_d):
        """

        :param catalog_n:
        :param i_pm:
        :param i_pm_alpha:
        :param i_pm_delta:
        :param i_alpha:
        :param i_delta:
        :param o_df:
        :param tmp_d:
        :return:
        """
        tmp_d['sex_cf'].append(self.sex_cf)
        tmp_d['scmp_cf'].append(self.scmp_cf)
        if o_df['SOURCE_NUMBER'].size != 1:
            tmp_d['boolean_l'].append('False')
        else:
            tmp_d['boolean_l'].append('True')
        tmp_d['catalog'].append(catalog_n)
        source = int(o_df['SOURCE_NUMBER'].iloc[0])
        tmp_d['source'].append(source)
        epoch = float(o_df['EPOCH'].iloc[0])
        tmp_d['epoch'].append(epoch)
        tmp_d['i_pm'].append(i_pm)
        tmp_d['i_pm_alpha'].append(i_pm_alpha)
        tmp_d['i_pm_delta'].append(i_pm_delta)
        tmp_d['i_alpha'].append(i_alpha)
        tmp_d['i_delta'].append(i_delta)
        # Extracted proper motions
        pm_alpha = float(o_df['PMALPHA'])
        tmp_d['o_pm_alpha'].append(pm_alpha)
        pm_delta = float(o_df['PMDELTA'])
        tmp_d['o_pm_delta'].append(pm_delta)
        # Errors related to proper motion
        pm_alpha_error = float(o_df['PMALPHAERR'])
        tmp_d['o_pm_alpha_err'].append(pm_alpha_error)
        pm_delta_error = float(o_df['PMDELTAERR'])
        tmp_d['o_pm_delta_err'].append(pm_delta_error)
        o_alpha = float(o_df['ALPHA_J2000'].iloc[0])
        tmp_d['o_alpha'].append(o_alpha)
        o_delta = float(o_df['DELTA_J2000'].iloc[0])
        tmp_d['o_delta'].append(o_delta)
        errora_world = float(o_df['ERRA_WORLD'].iloc[0])
        tmp_d['error_a'].append(errora_world)
        errorb_world = float(o_df['ERRB_WORLD'].iloc[0])
        tmp_d['error_b'].append(errorb_world)

        return tmp_d

    def populate_err_dict(self, catalog_n, i_pm, i_pm_alpha, i_alpha,
                          i_pm_delta, i_delta, o_pm, o_pm_alpha, o_pm_delta,
                          cat_df, err_d, scmp_source, o_alpha, o_delta):
        """

        :param catalog_n:
        :param i_pm:
        :param i_pm_alpha:
        :param i_alpha:
        :param i_pm_delta:
        :param i_delta:
        :param o_pm:
        :param o_pm_alpha:
        :param o_pm_delta:
        :param cat_df:
        :param err_d:
        :param scmp_source:
        :param o_alpha:
        :param o_delta:
        :return:
        """
        err_d['sex_cf'].append(self.sex_cf)
        err_d['scmp_cf'].append(self.scmp_cf)
        err_d['catalog'].append(catalog_n)
        source = int(cat_df['source'].iloc[0])
        err_d['source'].append(source)
        err_d['scmp_source'].append(scmp_source)
        err_d['i_pm'].append(i_pm)
        err_d['i_pm_alpha'].append(i_pm_alpha)
        err_d['i_pm_delta'].append(i_pm_delta)
        err_d['i_alpha'].append(i_alpha)
        err_d['i_delta'].append(i_delta)
        if o_pm is False:
            err_d['o_pm'].append(False)
        else:
            err_d['o_pm'].append(o_pm)
        if o_pm_alpha is False:
            err_d['o_pm_alpha'].append(False)
        else:
            err_d['o_pm_alpha'].append(o_pm_alpha)
        if o_pm_delta is False:
            err_d['o_pm_delta'].append(False)
        else:
            err_d['o_pm_delta'].append(o_pm_delta)
        if o_alpha is False:
            err_d['o_alpha'].append(False)
        else:
            err_d['o_alpha'].append(o_alpha)
        if o_delta is False:
            err_d['o_delta'].append(False)
        else:
            err_d['o_delta'].append(o_delta)

        return err_d

    def get_ccd(self, catalog_n):
        """

        :param catalog_n:
        :return: ccd
        """
        cats = [['x0_y0', 1, 1], ['x0_y0', 2, 2], ['x0_y0', 3, 3],
                ['x0_y0', 4, 4], ['x0_y1', 1, 5], ['x0_y1', 2, 6],
                ['x0_y1', 3, 7], ['x0_y1', 4, 8], ['x0_y2', 1, 9],
                ['x0_y2', 2, 10], ['x0_y2', 3, 11], ['x0_y2', 4, 12],
                ['x1_y0', 1, 13], ['x1_y0', 2, 14], ['x1_y0', 3, 15],
                ['x1_y0', 4, 16], ['x1_y1', 1, 17], ['x1_y1', 2, 18],
                ['x1_y1', 3, 19], ['x1_y1', 4, 20], ['x1_y2', 1, 21],
                ['x1_y2', 2, 22], ['x1_y2', 3, 23], ['x1_y2', 4, 24],
                ['x2_y0', 1, 25], ['x2_y0', 2, 26], ['x2_y0', 3, 27],
                ['x2_y0', 4, 28], ['x2_y1', 1, 29], ['x2_y1', 2, 30],
                ['x2_y1', 3, 31], ['x2_y1', 4, 32], ['x2_y2', 1, 33],
                ['x2_y2', 2, 34], ['x2_y2', 3, 35], ['x2_y2', 4, 36]]

        ccd = ''
        dither = ''
        for cat_ in cats:
            if cat_[2] == catalog_n:
                ccd = cat_[0]
                dither = cat_[1]

        fits_loc = '{}/{}/CCDs'.format(self.prfs_d['fits_dir'], self.mag)
        fits_name = 'mag_{}_CCD_{}_d{}.fits'.format(self.mag, ccd, dither)
        fits_file = '{}/{}'.format(fits_loc, fits_name)

        return fits_file
