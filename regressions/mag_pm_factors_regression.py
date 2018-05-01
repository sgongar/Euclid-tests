# !/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: First version. Plotting functions for completeness and purity factors
       are working. Creates two pdf files, one for each factor.

Todo:
    * Create a logger instance for this particular script.
    * Reshape scales in order to improve legibility. (v 0.2)
    * Add new factor values of new s-pline fit. (v 0.2)
    * Unit testing!

*GNU Terry Pratchett*
"""
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import array, meshgrid


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def f(x_, y_, b0, b1, intercept):
    """

    :param x_:
    :param y_:
    :param b0:
    :param b1:
    :param intercept:
    :return:
    """

    return (b0 * x_) + (b1 * y_) + intercept


analysis = ['purity', 'completeness']
coefficients = {'purity': {'b0': -0.117821, 'b1': -0.019716,
                           'intercept': 3.298143},
                'completeness': {'b0': -0.02146, 'b1': 0.15938,
                                 'intercept': 1.12410}}
factors = {'purity': [0.77, 0.23, 0.13, 0.81, 1.0, 1.0, 1.0, 1.0, 0.67, 0.19,
                      0.25, 0.91, 1.0, 1.0, 1.0, 1.0, 0.57, 0.2, 0.24, 1.0,
                      1.0, 1.0, 1.0, 1.0, 0.72, 0.25, 0.19, 0.62, 1.0, 1.0,
                      1.0, 1.0, 1.0, 0.5, 0.18, 0.23, 1.0, 1.0, 1.0, 1.0, 1.0,
                      0.33, 0.08, 0.04, 0.4, 1.0, 1.0, 1.0],
           'completeness': [0.81, 0.91, 0.86, 0.76, 0.68, 0.77, 0.78, 0.68,
                            0.28, 0.0, 0.82, 0.81, 0.95, 0.88, 0.81, 0.83,
                            0.89, 0.71, 0.29, 0.0, 0.72, 0.79, 0.95, 0.73,
                            0.81, 0.91, 0.75, 0.8, 0.24, 0.0, 0.62, 0.68, 0.8,
                            0.57, 0.68, 0.63, 0.25, 0.59, 0.15, 0.0, 0.32,
                            0.58, 0.7, 0.71, 0.55, 0.53, 0.24, 0.55, 0.0, 0.0,
                            0.07, 0.17, 0.32, 0.19, 0.24, 0.14, 0.09, 0.33,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0]}

mags = [20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5,
        21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5,
        22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5,
        23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5,
        24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5,
        25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5,
        26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5]

pms = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001,
       0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003,
       0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003, 0.01,
       0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003, 0.01, 0.03,
       0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003, 0.01, 0.03, 0.1,
       0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3,
       1.0, 3.0, 10.0, 30.0]


for analysis_ in analysis:

    with PdfPages('{}.pdf'.format(analysis_)) as pdf:
        print('running analysis for {} factor'.format(analysis_))

        fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        mags_tmp = []
        pms_tmp = []
        if analysis_ == 'completeness':
            mags_tmp = [x for x in mags if x != 26.5]
            pms_tmp = [x for x in pms if x != 10.0]
            pms_tmp = [x for x in pms_tmp if x != 30.0]

        elif analysis_ == 'purity':
            mags_tmp = mags
            pms_tmp = [x for x in pms if x != 10.0]
            pms_tmp = [x for x in pms_tmp if x != 30.0]

        x = array(mags_tmp)
        y = array(pms_tmp)

        X, Y = meshgrid(x, y)
        Z = f(X, Y, coefficients[analysis_]['b0'],
              coefficients[analysis_]['b1'],
              coefficients[analysis_]['intercept'])

        CS = plt.contourf(X, Y, Z, 20, cmap='RdGy')
        # print(CS.levels)  Todo - Reshaped levels
        # print(CS.levels[::2])
        plt.clabel(CS, inline=True, fontsize=12, colors='b')

        plt.colorbar()
        plt.axis(aspect='image')
        plt.grid(True)

        ax.set_title('analysis type: {}'.format(analysis_))
        ax.set_xlabel('Magnitude')
        ax.set_ylabel('Proper motion')

        pdf.savefig()
