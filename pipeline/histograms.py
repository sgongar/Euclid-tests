#!/usr/bin/python
# -*- coding: utf-8 -*-

""" #TODO

Versions:
* 0.1 - 

In order to improve #TODO
* c = catalog
* n = name
* d = dictionary
* loc = location
* cmp = comparation

Todo:
    * Improve log messages
    * Improve docstring
    * Create versions history
    * License??
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pandas import read_csv

# from misc import get_ticks, pm_compute, pm_filter, get_limits
from misc import extract_settings, get_fits

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class DrawHistograms:

    def __init__(self):
        """

        """
        prfs_d = extract_settings()
        self.read_stats(prfs_d)

    def read_stats(self, prfs_d):
        """ Reads stats.

        """
        # Creates dict.
        stats_d = {}
        # Gets fits files from directory.
        fits_files = get_fits(unique=False)

        for fits_n in fits_files:
            csv_n = '{}/sources_{}.csv'.format(prfs_d['tmp_out'],
                                               fits_n[-13:-5])
            test = read_csv(csv_n, index_col=0)
            print test

    def draw_histograms(self):
        """

        """
        fig, ((ax0, ax1),
              (ax2, ax3),
              (ax4, ax5)) = plt.subplots(ncols=2, nrows=3,
                                         figsize=(8.27, 11.69), dpi=100)


#     # input analysis values
#     ax0.hist(pm_in_values, bins=20,
#              facecolor='gray', label='test')
#     ax0.set_title("Proper motion distribution")
#     ax0.set_xlim(x0_limits)
#     x0_ticks.insert(0, 20)  # TODO workaround about automatic range
#     ax0.set_xticks(x0_ticks)
#     ax0.grid(True)

#     ax1.hist(mag_in_values, bins=20,
#              facecolor='gray')
#     ax1.set_title("Magnitude distribution")
#     ax1.set_xlim(x1_limits)
#     ax1.set_xticks(x1_ticks)

#     ax1.grid(True)

#     # output analysis values
#     print "proper motion out values", len(pm_out_values)
#     ax2.hist(pm_out_values, bins=20,
#              facecolor='gray', label='test')
#     ax2.set_title("Proper motion distribution")
#     ax2.set_xlim(x0_limits)
#     ax2.set_xticks(x0_ticks)
#     ax2.grid(True)

#     print "magnitude out values", len(mag_out_values)
#     ax3.hist(mag_out_values, bins=20,
#              facecolor='gray')
#     ax3.set_title("Magnitude distribution")
#     ax3.set_xlim(x1_limits)
#     ax3.set_xticks(x1_ticks)

#     ax3.grid(True)

#     print "flags distribution", len(flags_values)
#     ax4.hist(flags_values, bins=20,
#              facecolor='gray')
#     ax4.set_title("Flags distribution")
#     # ax4.set_xlim(x1_limits)
#     # ax4.set_xticks(x1_ticks)

#     ax4.grid(True)

#     print "flags distribution", len(flags_dict.keys())
#     labels = flags_dict.keys()
#     sizes = [len(flags_dict[flags_dict.keys()[x]]) for x in range(len(flags_dict.keys()))]
#     ax5.pie(sizes, labels=labels, autopct='%1.1f%%',
#             shadow=True, startangle=90)
#     ax5.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

#     ax5.grid(True)

#     fig.tight_layout()

#     with PdfPages('test.pdf') as pdf:
#         pdf.savefig()


if __name__ == '__main__':
    DrawHistograms()
