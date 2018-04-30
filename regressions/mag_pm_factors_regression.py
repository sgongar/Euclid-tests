# !/usr/bin/python
# -*- coding: utf-8 -*-

"""

Todo:
-

Versions:
- 0.1

Todo:
    * Improve log messages

*GNU Terry Pratchett*
"""
import matplotlib.pyplot as plt
from numpy import array, meshgrid


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


# def f(x_, y_):
#     """ coefficient number from R regression
#
#     :param x_:
#     :param y_:
#     :return:
#     """
#
#     return (-0.117821 * x_) + (-0.019716 * y_) + 3.298143
#
#
# # todo get real mag number!
# mags = [20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5,
#         21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5,
#         22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5,
#         23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5,
#         24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5,
#         25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5,
#         26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5]
#
# pms = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001,
#        0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003,
#        0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003, 0.01,
#        0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003, 0.01, 0.03,
#        0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003, 0.01, 0.03, 0.1,
#        0.3, 1.0, 3.0, 10.0, 30.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3,
#        1.0, 3.0, 10.0, 30.0]
#
# # just to check it!
# f_com = [0.81, 0.91, 0.86, 0.76, 0.68, 0.77, 0.78, 0.68, 0.28, 0.0, 0.82,
#          0.81, 0.95, 0.88, 0.81, 0.83, 0.89, 0.71, 0.29, 0.0, 0.72, 0.79,
#          0.95, 0.73, 0.81, 0.91, 0.75, 0.8, 0.24, 0.0, 0.62, 0.68, 0.8,
#          0.57, 0.68, 0.63, 0.25, 0.59, 0.15, 0.0, 0.32, 0.58, 0.7, 0.71,
#          0.55, 0.53, 0.24, 0.55, 0.0, 0.0, 0.07, 0.17, 0.32, 0.19, 0.24,
#          0.14, 0.09, 0.33, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#          0.0, 0.0, 0.0]
#
# x = array(mags)
# y = array(pms)
#
# X, Y = meshgrid(x, y)
# Z = f(X, Y)
# print(type(Z), Z)
#
# contours = plt.contour(X, Y, Z, 20, cmap='RdGy')
# plt.clabel(contours, inline=True, fontsize=12)
# plt.colorbar()
# plt.show()


def f(x_, y_):
    """ coefficient number from R regression

    :param x_:
    :param y_:
    :return:
    """

    return (-0.13987 * x_) + (-0.01809 * y_) + 3.83190


# todo get real mag number!
mags = [20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5, 20.5,
        21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5, 21.5,
        22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5, 22.5,
        23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5, 23.5,
        24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5, 24.5,
        25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5, 25.5,
        26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5, 26.5]
pms = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0,
       0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0,
       0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0,
       0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0,
       0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0,
       0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0,
       0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
f_com = [0.81, 0.91, 0.86, 0.76, 0.68, 0.77, 0.78, 0.68,
         0.82, 0.81, 0.95, 0.88, 0.81, 0.83, 0.89, 0.71,
         0.72, 0.79, 0.95, 0.73, 0.81, 0.91, 0.75, 0.8,
         0.62, 0.68, 0.8, 0.57, 0.68, 0.63, 0.25, 0.59,
         0.32, 0.58, 0.7, 0.71, 0.55, 0.53, 0.24, 0.55,
         0.07, 0.17, 0.32, 0.19, 0.24, 0.14, 0.09, 0.33,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
x = array(mags)
y = array(pms)

X, Y = meshgrid(x, y)
Z = f(X, Y)
print(type(Z), Z)

contours = plt.contourf(X, Y, Z, 20, cmap='RdGy')

plt.colorbar()
plt.axis(aspect='image')
plt.show()

