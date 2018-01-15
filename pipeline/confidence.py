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


def f(b, x):
    """

    :param b: is a vector of the parameters.
    :param x: is an array of the current x values. is in the same format as
    the x passed to Data or RealData.
    :return: an array in the same format as y passed to RealData.
    """

    return b[0]*x + b[1]


class WLSConfidence:

    def __init__(self, source):
        self.source = source

        # Extracts fitting values
        self.compute()

    def get_source_data(self):
        catalogs_dir = self.prfs_d['catalogs_dir']
        configuration_dir = '/20-21/{}/{}'.format(self.sex_cf, self.scmp_cf)
        cat_name = '/full_{}_20-21_1.cat'.format(self.scmp_cf)
        hdu_list = fits.open('{}{}{}'.format(catalogs_dir,
                                             configuration_dir, cat_name))
        db = Table(hdu_list[2].data).to_pandas()

        ra = db.loc[db['SOURCE_NUMBER'] == self.source, 'ALPHA_J2000'].tolist()
        dec = db.loc[db['SOURCE_NUMBER'] == self.source, 'DELTA_J2000'].tolist()
        # epoch = db.loc[db['SOURCE_NUMBER'] == source, 'EPOCH'].tolist()

    def compute(self):
        """

        :return:
        """


        x = array(ra)
        y = array(dec)
        sigma = db.loc[db['SOURCE_NUMBER'] == source, 'ERRA_WORLD'].tolist()
        x = sm.add_constant(x)
        edim = array([1 / var for var in sigma])
        model = sm.WLS(y, x, weigths=edim)
        fitted = model.fit()

    return fitted.rsquared


def odr_confidence(self, star, source, tmp_d):
    """

    :param star:
    :param source:
    :param tmp_d:
    :return:
    """
    catalogs_dir = self.prfs_d['catalogs_dir']
    configuration_dir = '/20-21/{}/{}'.format(self.sex_cf, self.scmp_cf)
    cat_name = '/full_{}_20-21_1.cat'.format(self.scmp_cf)
    hdu_list = fits.open('{}{}{}'.format(catalogs_dir,
                                         configuration_dir, cat_name))
    db = Table(hdu_list[2].data).to_pandas()

    ra = db.loc[db['SOURCE_NUMBER'] == source, 'ALPHA_J2000'].tolist()
    x = array(ra)
    dec = db.loc[db['SOURCE_NUMBER'] == source, 'DELTA_J2000'].tolist()
    y = array(dec)
    # epoch = db.loc[db['SOURCE_NUMBER'] == source, 'EPOCH'].tolist()

    err_ra = db.loc[db['SOURCE_NUMBER'] == source, 'ERRA_WORLD'].tolist()
    sx = array([1 / var for var in err_ra])
    err_dec = db.loc[db['SOURCE_NUMBER'] == source, 'ERRB_WORLD'].tolist()
    sy = array([1 / var for var in err_dec])

    linreg = linregress(x, y)

    linear = Model(f)
    mydata = RealData(x, y, sx=sx, sy=sy)
    tmp_linreg = list(linreg)  # temporary list for fitting output
    myodr = ODR(mydata, linear,
                beta0=array([tmp_linreg[0], tmp_linreg[1]]),
                delta0=err_ra, job=0, maxit=100)

    myoutput = myodr.run()  # run Orthogonal Distance Regression

    # print fit parameters and 1-sigma estimates
    popt = myoutput.beta
    perr = myoutput.sd_beta

    # prepare confidence level curves
    nstd = 5.  # to draw 5-sigma intervals
    # print('interval {}'.format(nstd * perr))
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr

    x_fit = linspace(min(x), max(x), 100)
    fit = f(popt, x_fit)
    fit_up = f(popt_up, x_fit)
    fit_dw = f(popt_dw, x_fit)

    cov_ra = float(myoutput.cov_beta.item((0, 0)))
    cov_dec = float(myoutput.cov_beta.item((1, 1)))
    cov_radec = float(myoutput.cov_beta.item((0, 1)))
    fit_odr = float(cov_radec / sqrt((cov_ra * cov_dec)))
    print('fit_odr {}'.format(fit_odr))

    try:
        dependant_var = float(myoutput.cov_beta.item((1, 1)))
        residual_var = myoutput.res_var
        coff_det_1 = 1 - (residual_var / dependant_var)
        print('coff_det_1 {}'.format(coff_det_1))
    except ZeroDivisionError:
        coff_det_1 = 0

    try:
        print('1 {}'.format(std(ra)))
        std_ra = std(ra) ** 2
        print('2 {}'.format(std_ra))
        std_dec = std(dec) ** 2
        cov_radec = float(myoutput.cov_beta.item((0, 1)))
        coff_det_2 = float(cov_radec ** 2 / (std_ra * std_dec))
        print('coff_det_2 {}'.format(coff_det_2))
    except ZeroDivisionError:
        coff_det_2 = 0

    a, b = myoutput.beta
    sa, sb = myoutput.sd_beta

    # Prints useful information about correlation
    # print(myoutput.pprint())

    predicted_dec = []
    for ra_ in ra:
        slope = float(list(linreg)[0])
        intercept = float(list(linreg)[1])
        predicted_dec.append((slope * ra_) + intercept)

    # print('dec {}'.format(dec))
    # print('predicted_dec {}'.format(predicted_dec))

    xp = linspace(min(ra), max(ra), 1000)
    yp = a * xp + b

    fits_d = {'x_odr': xp, 'y_odr': yp, 'fit_odr': fit_odr}
    plots_dir = self.prfs_d['plots_dir']

    # Moves output to different directory
    if len(tmp_d['i_pm']) == 0:
        output_path = '{}/fit/{}/{}/star/{}'.format(plots_dir, self.scmp_cf,
                                                    self.sex_cf,
                                                    tmp_d['o_pm_norm'][0])
        create_folder(self.logger, output_path)
    elif len(tmp_d['i_pm']) != 0:
        output_path = '{}/fit/{}/{}/SSO/{}'.format(plots_dir, self.scmp_cf,
                                                   self.sex_cf,
                                                   tmp_d['i_pm'][0])
        create_folder(self.logger, output_path)
    else:
        raise Exception
    try:
        plot = PlotFitting(star, tmp_d, output_path, fits_d)
    except UnboundLocalError:
        pass
    return myoutput
