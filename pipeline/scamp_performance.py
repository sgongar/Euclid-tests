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
from numpy import array
from pandas import concat, read_csv, DataFrame, Series
import statsmodels.api as sm

from misc import all_same, extract_settings
from misc import create_folder, speeds_range
from plots import PlotConfidence, PlotError
from regions import Create_regions


def redo_tmp_d():
    """ Creates a dictionary

    :return: tmp_d
    """
    tmp_d = {'sex_cf': [], 'scmp_cf': [], 'boolean_l': [], 'catalog': [],
             'CCD': [], 'source': [], 'epoch': [], 'i_pm': [],
             'i_pm_alpha': [], 'i_pm_delta': [], 'i_alpha': [],
             'i_delta': [], 'o_pm_alpha': [], 'o_pm_delta': [],
             'o_pm_alpha_err': [], 'o_pm_delta_err': [], 'o_alpha': [],
             'o_delta': [], 'error_a': [], 'error_b': []}

    return tmp_d


def redo_err_d():
    """ Creates a dictionary

    :return: err_d
    """
    err_d = {'sex_cf': [], 'scmp_cf': [], 'catalog': [],
             'CCD': [], 'source': [], 'scmp_source': [], 'i_pm': [],
             'i_pm_alpha': [], 'i_pm_delta': [], 'o_pm': [], 'o_pm_alpha': [],
             'o_pm_delta': [], 'i_alpha': [], 'i_delta': [], 'o_alpha': [],
             'o_delta': []}

    return err_d


def check_star(catalog_n, i_df, o_alpha, o_delta):
    """

    :param catalog_n:
    :param i_df:0
    :param o_alpha:
    :param o_delta:
    :return:
    """
    tolerance = 0.0005  # was 0.0001

    i_df = i_df[i_df['catalog'].isin([catalog_n])]
    i_df = i_df[i_df['alpha_j2000'] + tolerance > o_alpha]
    i_df = i_df[o_alpha > i_df['alpha_j2000'] - tolerance]
    i_df = i_df[i_df['delta_j2000'] + tolerance > o_delta]
    i_df = i_df[o_delta > i_df['delta_j2000'] - tolerance]

    return i_df


def cut_catalog(o_cat, margin, limits):
    """

    :param o_cat:
    :param margin:
    :param limits:
    :return:
    """
    o_df = o_cat[o_cat['ALPHA_J2000'] + margin > limits['max_alpha']]
    o_df = o_df[limits['min_alpha'] > o_df['ALPHA_J2000'] - margin]
    o_df = o_df[o_df['DELTA_J2000'] + margin > limits['max_delta']]
    o_df = o_df[limits['min_delta'] > o_df['DELTA_J2000'] - margin]

    return o_df


def check_cat_order(cat_list):
    """

    :param cat_list:
    :return:
    """
    tmp_order = []

    for idx_cat in range(0, len(cat_list), 1):
        cat = cat_list[idx_cat]
        dither = cat - (cat / 4) * 4
        if dither == 0:
            dither = 4
        tmp_order.append(dither)

    if sorted(tmp_order) == tmp_order:
        return True
    else:
        return False


def check_source(catalog_n, o_cat, i_alpha, i_delta):
    """

    :param catalog_n:
    :param o_cat:
    :param i_alpha:
    :param i_delta:
    :return:
    """
    tolerance = 0.0005  # 1.8 arcsecond

    o_df = o_cat[o_cat['CATALOG_NUMBER'].isin([catalog_n])]
    o_df = o_df[o_df['ALPHA_J2000'] + tolerance > i_alpha]
    o_df = o_df[i_alpha > o_df['ALPHA_J2000'] - tolerance]
    o_df = o_df[o_df['DELTA_J2000'] + tolerance > i_delta]
    o_df = o_df[i_delta > o_df['DELTA_J2000'] - tolerance]

    return o_df


def create_dict(scmp_cf, sex_cf, confidence_):
    """

    :param scmp_cf:
    :param sex_cf:
    :param confidence_:
    :return:
    """
    stats_keys = ['total', 'right', 'false', 'f_dr', 'f_pur', 'f_com']

    stats_d = {'PM': [0.001, 0.003, 0.01, 0.03, 0.1, 0.3,
                      1, 3, 10, 30]}

    scamp_parameters = scmp_cf.split('_')
    sex_parameters = sex_cf.split('_')

    stats_d['crossid'] = []
    stats_d['pixscale'] = []
    stats_d['posangle'] = []
    stats_d['position'] = []
    stats_d['deblending'] = []
    stats_d['threshold'] = []
    stats_d['mincount'] = []
    stats_d['area'] = []
    stats_d['confidence'] = []

    for value_ in range(len(stats_d['PM'])):
        stats_d['crossid'].append(scamp_parameters[0])
        stats_d['pixscale'].append(scamp_parameters[1])
        stats_d['posangle'].append(scamp_parameters[2])
        stats_d['position'].append(scamp_parameters[3])
        stats_d['deblending'].append(sex_parameters[0])
        stats_d['threshold'].append(sex_parameters[1])
        stats_d['mincount'].append(sex_parameters[3])
        stats_d['area'].append(sex_parameters[4])
        # Confidence
        stats_d['confidence'].append(confidence_)

    for key_ in stats_keys:
        stats_d[key_] = []
        for value_ in range(len(stats_d['PM'])):
            stats_d[key_].append(0)

    # out dictionary
    out_keys = ['alpha_j2000', 'delta_j2000',
                'catalog', 'PM', 'source', 'CCD', 'dither']
    out_d = {}

    for key_ in out_keys:
        out_d[key_] = []

    return stats_d, out_d


class ScampPerformanceSSOs:

    def __init__(self, logger, mag, scmp_cf, sex_cf, confidence_):
        """

        """
        self.bypassfilter = True
        self.checkerror = True
        self.prfs_d = extract_settings()

        self.mag = mag
        self.scmp_cf = scmp_cf
        self.sex_cf = sex_cf
        self.confidence_ = confidence_

        self.check(logger)

    def check(self, logger):
        """

        :param logger:
        :return:
        """
        # For now any file is saved
        save = False

        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {} and {}'.format(self.scmp_cf,
                                                                 self.sex_cf))

        input_d = {}
        for d in range(1, 5, 1):
            cat_loc = '{}/{}/Catalogs'.format(self.prfs_d['fits_dir'],
                                              self.mag)
            cat_name = '{}/Cat_20-21_d{}'.format(cat_loc, d)
            input_d[d] = '{}.dat'.format(cat_name)
        input_d = Create_regions(input_d).check_luca(self.mag, False, True)

        # Creates a DataFrame from an input dictionary
        input_l = []
        for key_ in input_d.keys():
            input_l.append(input_d[key_])

        i_df = concat(input_l, axis=0)
        # Look for < 3 coincidences
        i_df = concat(g for _, g in i_df.groupby('source')
                      if len(g) >= 3)
        i_df = i_df.reset_index(drop=True)
        i_df.to_csv('inputs.csv')

        if save:
            self.create_regions_files(i_df)

        # Open particular file!
        filt_n = 'filt_{}_{}_2.csv'.format(self.scmp_cf, self.mag)
        filter_o_n = '{}/{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                             self.mag, self.sex_cf,
                                             self.scmp_cf, filt_n)

        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        stats_d, out_d = create_dict(self.scmp_cf, self.sex_cf,
                                     self.confidence_)

        # Gets unique sources from input data
        unique_sources = list(set(i_df['source'].tolist()))

        unique_sources = [24, 274, 284, 314]
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
                # source_ = row.source
                # ccd_ = row.CCD
                # dither_ = row.dither_values
                catalog_n = row.catalog
                i_pm = row.pm_values
                i_pm_alpha = row.pm_alpha  # this is wrong! fixme
                i_pm_delta = row.pm_delta
                i_alpha = row.alpha_j2000
                i_delta = row.delta_j2000

                # Checks if there is a source closed to input one
                # o_cat contains data from output (filtered) catalog
                o_df = check_source(catalog_n, o_cat, i_alpha, i_delta)

                # If there is one saves data from input data
                if o_df.empty is not True and o_df['PM'].size == 1:
                    pm_mask = self.pm_filter(o_df, i_pm, self.confidence_)
                    if pm_mask:
                        if self.checkerror:
                            scmp_source = float(o_df['SOURCE_NUMBER'].iloc[0])
                            o_alpha = float(o_df['ALPHA_J2000'].iloc[0])
                            o_delta = float(o_df['DELTA_J2000'].iloc[0])
                            o_pm = float(o_df['PM'].iloc[0])
                            o_pm_alpha = float(o_df['PMALPHA'].iloc[0])
                            o_pm_delta = float(o_df['PMDELTA'].iloc[0])
                            print('o_alpha {}'.format(o_alpha))
                            print('o_delta {}'.format(o_delta))
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

            # Total number
            idx = stats_d['PM'].index(i_pm)
            stats_d['total'][idx] += 1

            order_mask = check_cat_order(tmp_d['catalog'])

            if not order_mask:
                raise Exception

            # plot errors!
            if len(err_d['scmp_source']) is not 0 and err_d['scmp_source'][0] is not False:
                flag_detection, sources_number = all_same(err_d['scmp_source'])

                if flag_detection and sources_number >= 3:
                    if self.checkerror:
                        print('yes')
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

                        fitted_d = self.confidence(err_d['scmp_source'][0],
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
                        print('2')
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

                        fitted_d = self.confidence(err_d['scmp_source'][0],
                                                   self.scmp_cf, self.sex_cf)

                        ok = '2'
                        self.plot_err(output_path, err_d, fits_, ok, fitted_d)
                    else:
                        pass  # nothing to do
            else:
                if self.checkerror:
                    print('3')
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

    def confidence(self, source, scmp_cf, sex_cf):
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

        edims = []
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
                edim = array([1 / var for var in sigma])
                edims.append(edim)
                x = sm.add_constant(x)
                model = sm.WLS(y, x, weigths=edim)
                fitted = model.fit()
                chi_squared = fitted.rsquared
                fitted_d[dimension_] = chi_squared
            elif dimension == dec:
                sigma = db.loc[db['SOURCE_NUMBER'] == source,
                               'ERRB_WORLD'].tolist()
                dimension_ = 'dec'
                edim = array([1 / var for var in sigma])
                edims.append(edim)
                x = sm.add_constant(x)
                model = sm.WLS(y, x, weigths=edim)
                fitted = model.fit()
                chi_squared = fitted.rsquared
                fitted_d[dimension_] = chi_squared

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


class ScampPerformanceStars:

    def __init__(self):
        """

        """
        self.bypassfilter = True
        self.prfs_d = extract_settings()

    def check(self, logger, mag, scmp_cf, sex_cf):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_cf:
        :return:
        """
        # Creates an input dictionary with all input sources
        logger.debug('checking performance for {} and {}'.format(scmp_cf,
                                                                 sex_cf))
        input_d = {}
        for d in range(1, 5, 1):
            cat_name = '{}/Cat_20-21_d{}'.format(self.prfs_d['input_ref'], d)
            input_d[d] = '{}.dat'.format(cat_name)
        input_d = Create_regions(input_d).check_luca(True, True)

        # Creates a DataFrame from an input dictionary
        input_l = []
        for key_ in input_d.keys():
            input_l.append(input_d[key_])

        i_df = concat(input_l, axis=0)
        # Look for < 3 coincidences
        i_df = concat(g for _, g in i_df.groupby('source')
                      if len(g) >= 3)
        i_df = i_df.reset_index(drop=True)

        # Open particular file!
        filt_n = 'filt_{}_{}_3.csv'.format(scmp_cf, mag)
        filter_o_n = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'],
                                          sex_cf, scmp_cf, filt_n)

        # Cross with filtered data - Opens datafile
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        # Gets unique sources from input data
        o_unique_sources = list(set(o_cat['SOURCE_NUMBER'].tolist()))

        # Creates variables
        right = 0
        not_right = 0
        total = 0
        stats_d = {'source_l': [], 'catalog_l': [], 'alpha_l': [],
                   'delta_l': [], 'pm_l': [], 'mag': [], 'error_a': [],
                   'error_b': []}

        for source_ in o_unique_sources:
            flag = False
            tmp_d = {'source': [], 'catalog': [], 'i_pm': [],
                     'i_alpha': [], 'i_delta': [], 'mag': [],
                     'error_a': [], 'error_b': []}

            o_df = o_cat[o_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                source = row.SOURCE_NUMBER
                catalog_n = row.CATALOG_NUMBER
                i_alpha = row.ALPHA_J2000
                i_delta = row.DELTA_J2000
                pm = row.PM
                mag = row.MAG
                error_a = row.ERRA_WORLD
                error_b = row.ERRB_WORLD

                out_df = check_star(catalog_n, i_df, i_alpha, i_delta)

                if out_df.empty:
                    # it's a star
                    right += 1
                    flag = True
                    tmp_d['source'].append(source)
                    tmp_d['catalog'].append(catalog_n)
                    tmp_d['i_alpha'].append(i_alpha)
                    tmp_d['i_delta'].append(i_delta)
                    tmp_d['i_pm'].append(pm)
                    tmp_d['mag'].append(mag)
                    tmp_d['error_a'].append(error_a)
                    tmp_d['error_b'].append(error_b)
                else:
                    # it's a SSO
                    not_right += 1
                    flag = False

                total += 1

            if flag:
                order_mask = check_cat_order(tmp_d['catalog'])
                flag_detection, sources_number = all_same(tmp_d['source'])

                if order_mask and sources_number >= 3:
                    for tmp_source_ in tmp_d['source']:
                        stats_d['source_l'].append(tmp_source_)
                    for catalog_ in tmp_d['catalog']:
                        stats_d['catalog_l'].append(catalog_)
                    for i_alpha_ in tmp_d['i_alpha']:
                        stats_d['alpha_l'].append(i_alpha_)
                    for i_delta_ in tmp_d['i_delta']:
                        stats_d['delta_l'].append(i_delta_)
                    for pm_ in tmp_d['i_pm']:
                        stats_d['pm_l'].append(pm_)
                    for mag_ in tmp_d['mag']:
                        stats_d['mag'].append(mag_)
                    for erra_ in tmp_d['error_a']:
                        stats_d['erra_world'].append(erra_)
                    for errb_ in tmp_d['error_b']:
                        stats_d['errb_world'].append(errb_)

        stats_df = DataFrame(stats_d)
        stats_df.to_csv('test_{}_{}.csv'.format(sex_cf, scmp_cf))

        return stats_d
