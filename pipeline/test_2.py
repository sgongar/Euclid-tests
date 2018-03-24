#!/usr/bin/python
# -*- coding: utf-8 -*-


from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
from numpy import arange, mean, std, array, polyfit, sum, sqrt, linspace, size
from pandas import read_csv
from scipy.stats import t, linregress
from scipy.odr import Model, ODR, RealData


def scatterfit(x, y, a=None, b=None):
    """
   Compute the mean deviation of the data about the linear model given if A,B
   (y=ax+b) provided as arguments. Otherwise, compute the mean deviation about
   the best-fit line.

   x,y assumed to be Numpy arrays. a,b scalars.
   Returns the float sd with the mean deviation.

   Author: Rodrigo Nemmen
    """

    if a == None:
        # Performs linear regression
        a, b, r, p, err = linregress(x, y)

    # Std. deviation of an individual measurement (Bevington, eq. 6.15)
    N = size(x)
    sd = 1. / (N - 2.) * sum((y - a * x - b) ** 2)
    sd = sqrt(sd)

    return sd


def predband(xd, yd, a, b, conf=0.40, x=None):
    """
    Calculates the prediction band of the linear regression model at the desired confidence
    level.
    Clarification of the difference between confidence and prediction bands:
    "The 2sigma confidence interval is 95% sure to contain the best-fit regression line.
    This is not the same as saying it will contain 95% of the data points. The prediction bands are
    further from the best-fit line than the confidence bands, a lot further if you have many data
    points. The 95% prediction interval is the area in which you expect 95% of all data points to fall."
    (from http://graphpad.com/curvefit/linear_regression.htm)
    Arguments:
    - conf: desired confidence level, by default 0.95 (2 sigma)
    - xd,yd: data arrays
    - a,b: linear fit parameters as in y=ax+b
    - x: (optional) array with x values to calculate the confidence band. If none is provided, will
     by default generate 100 points in the original x-range of the data.

    Usage:
    >>> lpb,upb,x=nemmen.predband(all.kp,all.lg,a,b,conf=0.95)
    calculates the prediction bands for the given input arrays
    >>> pylab.fill_between(x, lpb, upb, alpha=0.3, facecolor='gray')
    plots a shaded area containing the prediction band
    Returns:
    Sequence (lpb,upb,x) with the arrays holding the lower and upper confidence bands
    corresponding to the [input] x array.
    References:
    1. http://www.JerryDallal.com/LHSP/slr.htm, Introduction to Simple Linear Regression, Gerard
    E. Dallal, Ph.D.
    Rodrigo Nemmen
    v1 Dec. 2011
    v2 Jun. 2012: corrected bug in dy.
    """
    alpha = 1. - conf  # significance
    n = xd.size  # data sample size
    if x == None:
        x = linspace(xd.min(), xd.max(), 100)
    # Predicted values (best-fit model)
    y = a * x + b
    # Auxiliary definitions
    sd = scatterfit(xd, yd, a, b)  # Scatter of data about the model
    sxd = sum((xd - xd.mean()) ** 2)
    sx = (x - xd.mean()) ** 2  # array
    # Quantile of Student's t distribution for p=1-alpha/2
    q = t.ppf(1. - alpha / 2., n - 2)
    # Prediction band
    dy = q * sd * sqrt(1. + 1. / n + sx / sxd)
    upb = y + dy  # Upper prediction band
    lpb = y - dy  # Lower prediction band
    return lpb, upb, x


def f(b, x):
    """
    Linear function y = m*x + b
    B is a vector of the parameters.
    x is an array of the current x values.
    x is in the same format as the x passed to Data or RealData.

    Return an array in the same format as y passed to Data or RealData.
    """
    return b[0] * x + b[1]


def b_image_plot():

    data_dir = '/home/sgongora/Documents/CarpetaCompartida/full_stats_ssos'

    b_image = []
    b_image_d = {}
    errb_image = []
    mag_iso = []
    mag_iso_d = {}
    magerr_iso = []

    # reads b_image
    for mag_ in ['20-21', '21-22', '22-23', '23-24', '24-25', '25-26', '26-27']:
        b_image_d[mag_] = []
        mag_iso_d[mag_] = []

        b_image_df = read_csv(
            '{}/f_{}_median_b_image_3.csv'.format(data_dir, mag_),
            index_col=0)
        print(b_image_df.columns[:-3])
        for column_ in b_image_df.columns[:-3]:
            list_to_append = b_image_df[column_].tolist()
            list_to_append = [x for x in list_to_append if type(x) != float]
            list_to_append = [x for x in list_to_append if 'mean' not in x]
            list_to_append = [float(x) for x in list_to_append]
            for element_ in list_to_append:
                b_image.append(element_)
                b_image_d[mag_].append(element_)

        mag_iso_df = read_csv(
            '{}/f_{}_median_mag_iso_3.csv'.format(data_dir, mag_),
            index_col=0)
        for column_ in mag_iso_df.columns[:-3]:
            list_to_append = mag_iso_df[column_].tolist()
            list_to_append = [x for x in list_to_append if type(x) != float]
            list_to_append = [x for x in list_to_append if 'mean' not in x]
            list_to_append = [float(x) for x in list_to_append]
            for element_ in list_to_append:
                mag_iso.append(element_)
                mag_iso_d[mag_].append(element_)

        errb_image_df = read_csv(
            '{}/f_{}_median_errb_image_3.csv'.format(data_dir,
                                                     mag_),
            index_col=0)
        for column_ in errb_image_df.columns[:-3]:
            list_to_append = errb_image_df[column_].tolist()
            list_to_append = [x for x in list_to_append if type(x) != float]
            list_to_append = [x for x in list_to_append if 'mean' not in x]
            list_to_append = [float(x) for x in list_to_append]
            for element_ in list_to_append:
                errb_image.append(element_)

        magerr_iso_df = read_csv(
            '{}/f_{}_median_magerr_iso_3.csv'.format(data_dir,
                                                     mag_),
            index_col=0)
        for column_ in magerr_iso_df.columns[:-3]:
            list_to_append = magerr_iso_df[column_].tolist()
            list_to_append = [x for x in list_to_append if type(x) != float]
            list_to_append = [x for x in list_to_append if 'mean' not in x]
            list_to_append = [float(x) for x in list_to_append]
            for element_ in list_to_append:
                magerr_iso.append(element_)

    test_fit = polyfit(mag_iso, b_image, 2)

    linear = Model(f)
    mydata = RealData(mag_iso, b_image, sx=magerr_iso, sy=errb_image)

    myodr = ODR(mydata, linear, beta0=[test_fit[0], test_fit[1]])

    myoutput = myodr.run()
    # myoutput.pprint()

    lower, upper, x = predband(array(mag_iso), array(b_image),
                               myoutput.beta[0], myoutput.beta[1])
    fit = myoutput.beta

    upper_fit = polyfit(x, upper, 1)
    lower_fit = polyfit(x, lower, 1)

    # reads mag_iso
    plot_size = [16.53, 11.69]
    plot_dpi = 100

    with PdfPages('mag_iso_vs_b_image_2.pdf') as pdf:
        fig = pyplot.figure(figsize=plot_size, dpi=plot_dpi)
        ax_1 = fig.add_subplot(1, 1, 1)
        ax_1.plot(mag_iso, b_image, 'bs')
        colors = ['gs', 'rs', 'cs', 'ms', 'ys', 'ks', 'ws']
        for idx_color, mag_ in enumerate(['20-21', '21-22', '22-23',
                                          '23-24', '24-25', '25-26',
                                          '26-27']):
            ax_1.plot(mag_iso_d[mag_], b_image_d[mag_],
                      colors[idx_color])
        ax_1.plot([19, 27], [19 * myoutput.beta[0] + myoutput.beta[1],
                             27 * myoutput.beta[0] + myoutput.beta[1]])
        ax_1.plot([19, 27], [19 * upper_fit[0] + upper_fit[1],
                             27 * upper_fit[0] + upper_fit[1]])
        ax_1.plot([19, 27], [19 * lower_fit[0] + lower_fit[1],
                             27 * lower_fit[0] + lower_fit[1]])
        ax_1.set_ylabel('b_image')
        ax_1.set_xlabel('mag_iso')
        ax_1.grid(True)

        pdf.savefig()

    """
    fig = pyplot.figure(figsize=plot_size, dpi=plot_dpi)
    ax_1 = fig.add_subplot(1, 1, 1)
    ax_1.plot(mag_iso, b_image, 'bs')
    ax_1.plot([19, 27], [19 * myoutput.beta[0] + myoutput.beta[1],
                         27 * myoutput.beta[0] + myoutput.beta[1]])
    ax_1.plot([19, 27], [19 * upper_fit[0] + upper_fit[1],
                         27 * upper_fit[0] + upper_fit[1]])
    ax_1.plot([19, 27], [19 * lower_fit[0] + lower_fit[1],
                         27 * lower_fit[0] + lower_fit[1]])
    ax_1.set_ylabel('b_image')
    ax_1.set_xlabel('mag_iso')
    ax_1.grid(True)

    pyplot.show()
    """

    return fit, lower_fit, upper_fit, b_image, mag_iso

def a_image_plot():

    data_dir = '/home/sgongora/Documents/CarpetaCompartida/full_stats_ssos'

    a_image = []
    erra_image = []
    mag_iso = []
    magerr_iso = []

    # reads b_image
    for mag_ in ['20-21', '21-22', '22-23']:
        a_image_df = read_csv(
            '{}/f_{}_median_a_image_3.csv'.format(data_dir, mag_),
            index_col=0)
        for column_ in a_image_df.columns:
            list_to_append = a_image_df[column_].tolist()
            list_to_append = [x for x in list_to_append if type(x) != float]
            list_to_append = [x for x in list_to_append if 'mean' not in x]
            list_to_append = [float(x) for x in list_to_append]
            for element_ in list_to_append:
                a_image.append(element_)

        mag_iso_df = read_csv(
            '{}/f_{}_median_mag_iso_3.csv'.format(data_dir, mag_),
            index_col=0)
        for column_ in mag_iso_df.columns:
            list_to_append = mag_iso_df[column_].tolist()
            list_to_append = [x for x in list_to_append if type(x) != float]
            list_to_append = [x for x in list_to_append if 'mean' not in x]
            list_to_append = [float(x) for x in list_to_append]
            for element_ in list_to_append:
                mag_iso.append(element_)

        erra_image_df = read_csv(
            '{}/f_{}_median_erra_image_3.csv'.format(data_dir,
                                                     mag_),
            index_col=0)
        for column_ in erra_image_df.columns:
            list_to_append = erra_image_df[column_].tolist()
            list_to_append = [x for x in list_to_append if type(x) != float]
            list_to_append = [x for x in list_to_append if 'mean' not in x]
            list_to_append = [float(x) for x in list_to_append]
            for element_ in list_to_append:
                erra_image.append(element_)

        magerr_iso_df = read_csv(
            '{}/f_{}_median_magerr_iso_3.csv'.format(data_dir,
                                                     mag_),
            index_col=0)
        for column_ in magerr_iso_df.columns:
            list_to_append = magerr_iso_df[column_].tolist()
            list_to_append = [x for x in list_to_append if type(x) != float]
            list_to_append = [x for x in list_to_append if 'mean' not in x]
            list_to_append = [float(x) for x in list_to_append]
            for element_ in list_to_append:
                magerr_iso.append(element_)

    test_fit = polyfit(mag_iso, a_image, 1)

    linear = Model(f)
    mydata = RealData(mag_iso, a_image, sx=magerr_iso, sy=erra_image)

    myodr = ODR(mydata, linear, beta0=[test_fit[0], test_fit[1]])

    myoutput = myodr.run()
    myoutput.pprint()

    # reads mag_iso
    plot_size = [16.53, 11.69]
    plot_dpi = 100

    with PdfPages('mag_iso_vs_b_image.pdf') as pdf:
        fig = pyplot.figure(figsize=plot_size, dpi=plot_dpi)
        ax_1 = fig.add_subplot(1, 1, 1)
        ax_1.plot(mag_iso, a_image, 'bs')
        ax_1.plot([19, 27], [19 * myoutput.beta[0] + myoutput.beta[1],
                             27 * myoutput.beta[0] + myoutput.beta[1]])
        ax_1.set_ylabel('a_image')
        ax_1.set_xlabel('mag_iso')
        ax_1.grid(True)

        pdf.savefig()


if __name__ == "__main__":
    fit, lower_fit, upper_fit, b_image, mag_iso = b_image_plot()

    print(fit)
    print(lower_fit)
    print(upper_fit)

    idx_in = 0
    idx_out = 0
    total = len(b_image)

    for idx, mag_iso_ in enumerate(mag_iso):
        upper_test = (upper_fit[0] * float(mag_iso_)) + upper_fit[1]
        lower_test = (lower_fit[0] * float(mag_iso_)) + lower_fit[1]

        if float(lower_test) < float(b_image[idx]) < float(upper_test):
            idx_in += 1
        else:
            idx_out += 1

    print('total: {}'.format(total))
    print('idx_in: {}'.format(idx_in))
    print('idx_out: {}'.format(idx_out))

    # print(b_image_test, upper_test, lower_test)
    # a_image_plot()