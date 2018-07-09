#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
For now intermediate versions are stored in physical memory, this will not
be longer needed in the future but in this moment could be useful for
pipeline development reasons.


Versions:
- 0.1: Initial release.

Todo:
    * Improve documentation

*GNU Terry Pratchett*
"""
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from numpy import arange
from pandas import read_csv

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"

# 75 % sources

"""
fig = pyplot.figure()
ax = Axes3D(fig)

x_vals = []
y_vals = []
z_vals = []
a_image_list = arange(0, 10, 0.15)
b_image_list = arange(0, 6, 0.15)
for a_image in a_image_list:
    for b_image in b_image_list:
        pm = 3.7 + (-4.54174735 * b_image) + (2.83112412 * a_image)
        x_vals.append(b_image)
        y_vals.append(a_image)
        z_vals.append(pm)
ax.scatter(x_vals, y_vals, z_vals, c='b')

x_vals = []
y_vals = []
z_vals = []
a_image_list = arange(0, 10, 0.15)
b_image_list = arange(0, 6, 0.15)
for a_image in a_image_list:
    for b_image in b_image_list:
        pm = 2.37524438132 + (-4.54174735 * b_image) + (2.83112412 * a_image)
        x_vals.append(b_image)
        y_vals.append(a_image)
        z_vals.append(pm)
ax.scatter(x_vals, y_vals, z_vals, c='r')

x_vals = []
y_vals = []
z_vals = []
a_image_list = arange(0, 10, 0.15)
b_image_list = arange(0, 6, 0.15)
for a_image in a_image_list:
    for b_image in b_image_list:
        pm = 0.5 + (-4.54174735 * b_image) + (2.83112412 * a_image)
        x_vals.append(b_image)
        y_vals.append(a_image)
        z_vals.append(pm)
ax.scatter(x_vals, y_vals, z_vals, c='b')
"""

margin = 0.5
total_ = 0
in_ = 0
out_ = 0
out_x_vals = []
out_y_vals = []
out_z_vals = []
catalogue = read_csv('catalogue.csv', index_col=0)
for i, row in enumerate(catalogue.itertuples(), 1):
    pm = 2.37524438132 + (-4.54174735 * row.B_IMAGE) + (2.83112412 * row.A_IMAGE)

    upper_value = (2.37524438132 * (1 + margin)) + (-4.54174735 * row.B_IMAGE) + (2.83112412 * row.A_IMAGE)
    lower_value = (2.37524438132 * (1 - margin)) + (-4.54174735 * row.B_IMAGE) + (2.83112412 * row.A_IMAGE)

    total_ += 1
    if lower_value < row.PM < upper_value:
        # print('in {}'.format(row.PM))
        in_ += 1
    else:
        # print('out {}'.format(row.PM))
        out_ += 1

    out_x_vals.append(row.B_IMAGE)
    out_y_vals.append(row.A_IMAGE)
    out_z_vals.append(row.PM)

print('total: {}'.format(total_))
print('in: {}'.format(in_))
print('out: {}'.format(out_))

"""
ax.scatter(out_x_vals, out_y_vals, out_z_vals, c='g')

pyplot.show()
"""
