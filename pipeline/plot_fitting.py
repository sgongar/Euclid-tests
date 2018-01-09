#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance

Versions:
- 0.1 Initial release

Todo:
    *

"""
from datetime import datetime, timedelta
from decimal import Decimal
from math import modf

from astropy import units
from astropy.coordinates import Angle

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from pyds9 import DS9

from misc import significant_l, extract_settings


class PlotFitting:

    def __init__(self, plot_d, output_path):
        """

        :param plot_d:
        """
        self.prfs_d = extract_settings()
        self.plot_d = plot_d
        self.output_path = output_path

        # Page size
        self.plot_size = [16.53, 11.69]
        self.plot_dpi = 100

        self.plot()

    def plot(self):
        """

        :return:
        """
        pdf_name = '{}/{}.pdf'.format(self.output_path, self.err_d['source'][0])

        with PdfPages(pdf_name) as pdf:
            # ALPHA PARAMETERS
            fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
            ax_1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)

    # def __init__(self, output_path, err_d, fits_, mag, ok, fitted_d):
    #     """
    #
    #     :param output_path:
    #     :param err_d:
    #     """
    #     self.prfs_d = extract_settings()
    #     self.output_path = output_path
    #     self.err_d = err_d
    #     self.fits_files = fits_
    #     self.ccds = self.gets_ccd_names(fits_)
    #     self.mag = mag
    #     self.ok = ok
    #     self.fitted_d = fitted_d
    #
    #     # Page size
    #     self.plot_size = [16.53, 11.69]
    #     self.plot_dpi = 100
    #
    #     self.plot()
    #
    # def gets_ccd_names(self, fits_):
    #     """
    #
    #     :return:
    #     """
    #     ccds = []
    #     for ccd in fits_:
    #         ccds.append(ccd[-13:-5])
    #
    #     return ccds
    #
    # def plot(self):
    #     """
    #
    #     :return:
    #     """
    #     pdf_name = '{}/{}.pdf'.format(self.output_path, self.err_d['source'][0])
    #
    #     with PdfPages(pdf_name) as pdf:
    #         # ALPHA PARAMETERS
    #         fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
    #         ax_1 = plt.subplot2grid((1, 5), (0, 0), colspan=4)
    #
    #         source_str = 'source {}'.format(self.err_d['source'][0])
    #         pm_str = 'input_pm {}'.format(self.err_d['i_pm'][0])
    #         # ok_str = 'ok {}'.format(self.ok)
    #         ax_1.set_title('{}\n{}'.format(source_str, pm_str))
    #
    #         # alpha coordinates
    #         i_alpha_tmpstmp = []
    #         i_alpha_seconds = []
    #         i_alpha_arcseconds = []
    #         tmp_hour = []
    #         tmp_minute = []
    #         for alpha_ in self.err_d['i_alpha']:
    #             a = Angle(alpha_, units.degree)
    #             dms = a.dms
    #             degree = int(dms[0])
    #             tmp_hour.append(degree)
    #             minute = int(dms[1])
    #             tmp_minute.append(minute)
    #             second = float("{0:.6f}".format(dms[2]))
    #             arcsecond = float("{0:.6f}".format(a.arcsecond))
    #             i_alpha_tmpstmp.append('{}:{}:{}'.format(degree, minute,
    #                                                      second))
    #             i_alpha_seconds.append('{}'.format(second))
    #             i_alpha_arcseconds.append(arcsecond)
    #
    #         o_alpha_tmpstmp = []
    #         o_alpha_seconds = []
    #         tmp_hour = []
    #         tmp_minute = []
    #         for alpha_ in self.err_d['o_alpha']:
    #             if alpha_ is not False:
    #                 a = Angle(alpha_, units.degree)
    #                 dms = a.dms
    #                 degree = int(dms[0])
    #                 tmp_hour.append(degree)
    #                 minute = int(dms[1])
    #                 tmp_minute.append(minute)
    #                 second = float("{0:.6f}".format(dms[2]))
    #                 o_alpha_tmpstmp.append('{}:{}:{}'.format(degree, minute,
    #                                                          second))
    #                 o_alpha_seconds.append('{}'.format(second))
    #
    #         alpha_seconds = i_alpha_seconds + o_alpha_seconds
    #
    #         # delta coordinates
    #         i_delta_tmpstmp = []
    #         i_delta_seconds = []
    #         tmp_degree = []
    #         tmp_minute = []
    #         for delta_ in self.err_d['i_delta']:
    #             d = Angle(delta_, units.degree)
    #             dms = d.dms
    #             degree = int(dms[0])
    #             tmp_degree.append(degree)
    #             minute = int(dms[1])
    #             tmp_minute.append(minute)
    #             second = float("{0:.6f}".format(dms[2]))
    #             i_delta_tmpstmp.append(
    #                 '{}:{}.{}'.format(degree, minute, second))
    #             i_delta_seconds.append('{}'.format(second))
    #
    #         o_delta_tmpstmp = []
    #         o_delta_seconds = []
    #         tmp_degree = []
    #         tmp_minute = []
    #         for delta_ in self.err_d['o_delta']:
    #             if delta_ is not False:
    #                 d = Angle(delta_, units.degree)
    #                 dms = d.dms
    #                 degree = int(dms[0])
    #                 tmp_degree.append(degree)
    #                 minute = int(dms[1])
    #                 tmp_minute.append(minute)
    #                 second = float("{0:.6f}".format(dms[2]))
    #                 o_delta_tmpstmp.append('{}:{}.{}'.format(degree, minute,
    #                                                          second))
    #                 o_delta_seconds.append('{}'.format(second))
    #
    #         delta_seconds = i_delta_seconds + o_delta_seconds
    #
    #         i_alpha_seconds = [float(i) for i in i_alpha_seconds]  # needed?
    #         i_alpha_seconds = [float("{0:.6f}".format(i)) for i in
    #                            i_alpha_seconds]
    #         i_delta_seconds = [float(i) for i in i_delta_seconds]  # needed?
    #         i_delta_seconds = [float("{0:.6f}".format(i)) for i in
    #                            i_delta_seconds]
    #
    #         o_alpha_seconds = [float(i) for i in o_alpha_seconds]  # needed?
    #         o_alpha_seconds = [float("{0:.6f}".format(i)) for i in
    #                            o_alpha_seconds]
    #         o_delta_seconds = [float(i) for i in o_delta_seconds]  # needed?
    #         o_delta_seconds = [float("{0:.6f}".format(i)) for i in
    #                            o_delta_seconds]
    #
    #         alpha_seconds = [float(i) for i in alpha_seconds]  # needed?
    #         alpha_seconds = [float("{0:.6f}".format(i)) for i in alpha_seconds]
    #         delta_seconds = [float(i) for i in delta_seconds]  # needed?
    #         delta_seconds = [float("{0:.6f}".format(i)) for i in delta_seconds]
    #
    #         # Plots data
    #         blue_colors = cm.Blues(np.linspace(0, 1, len(i_alpha_seconds) + 1))
    #         for idx in range(0, len(i_alpha_seconds), 1):
    #             ax_1.scatter(i_alpha_seconds[idx],
    #                          i_delta_seconds[idx],
    #                          c=blue_colors[idx + 1], s=32)
    #         """
    #         ax_1.plot(i_alpha_seconds, i_delta_seconds, 'bs', markersize=6,
    #                   label='catalog position')
    #         """
    #         # print('o_alpha {}'.format(o_alpha_seconds))
    #         # print('o_delta {}'.format(o_delta_seconds))
    #         # print('i_alpha {}'.format(i_alpha_seconds))
    #         # print('i_delta {}'.format(i_delta_seconds))
    #         """
    #         ax_1.plot(o_alpha_seconds, o_delta_seconds, 'rs', markersize=6,
    #                   label='extracted position')
    #         """
    #         red_colors = cm.Reds(np.linspace(0, 1, len(o_alpha_seconds) + 1))
    #         for idx in range(0, len(o_alpha_seconds), 1):
    #             ax_1.scatter(o_alpha_seconds[idx],
    #                          o_delta_seconds[idx],
    #                          c=red_colors[idx + 1], s=32)
    #
    #         try:
    #             x_ticks = create_y_ticks(alpha_seconds)
    #             y_ticks = create_y_ticks(delta_seconds)
    #         except ZeroDivisionError:
    #             print('alpha: {}'.format(self.err_d['i_alpha']))
    #             print('delta: {}'.format(self.err_d['i_delta']))
    #             print('source {}'.format(self.err_d['source'][0]))
    #             raise Exception
    #
    #         # x-ticks assignation
    #         ax_1.set_xticks(x_ticks['major_t'], minor=False)
    #         ax_1.set_xticklabels(x_ticks['major_t'])
    #         ax_1.set_xticks(x_ticks['minor_t'], minor=True)
    #
    #         # y-ticks assignation
    #         ax_1.set_yticks(y_ticks['major_t'], minor=False)
    #         ax_1.set_yticklabels(y_ticks['major_t'])
    #         ax_1.set_yticks(y_ticks['minor_t'], minor=True)
    #
    #         # Formats grids
    #         ax_1.grid(b=True, which='major', linestyle='-', linewidth=2)
    #         ax_1.grid(b=True, which='minor', linestyle='--', linewidth=1)
    #
    #         # x-axis
    #         x_label_ra = 'Right ascension \n'
    #         major_s = x_ticks['major_s']
    #         x_label_major_step = 'major step size {}"\n'.format(major_s)
    #         minor_s = x_ticks['minor_s']
    #         x_label_minor_step = 'minor step size {}"'.format(minor_s)
    #         ax_1.set_xlabel('{}{}{}'.format(x_label_ra, x_label_major_step,
    #                                         x_label_minor_step))
    #         # x-axis
    #         y_label_ra = 'Declination \n'
    #         major_s = y_ticks['major_s']
    #         y_label_major_step = 'major step size {}"\n'.format(major_s)
    #         minor_s = y_ticks['minor_s']
    #         y_label_minor_step = 'minor step size {}"'.format(minor_s)
    #         ax_1.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
    #                                         y_label_minor_step))
    #
    #         # Plots a legend
    #         plt.legend(loc=0, ncol=2, borderaxespad=0.)
    #         plt.setp(ax_1.get_xticklabels(), visible=True)
    #         plt.setp(ax_1.get_yticklabels(), visible=True)
    #         plt.draw()
    #
    #         ax_2 = plt.subplot2grid((1, 5), (0, 4))
    #
    #         i_pm = float("{0:.6f}".format(self.err_d['i_pm'][1]))
    #         i_pm_alpha = float("{0:.6f}".format(self.err_d['i_pm_alpha'][1]))
    #         i_pm_delta = float("{0:.6f}".format(self.err_d['i_pm_delta'][1]))
    #         o_pm = float("{0:.6f}".format(self.err_d['o_pm'][1]))
    #         o_pm_alpha = float("{0:.6f}".format(self.err_d['o_pm_alpha'][1]))
    #         o_pm_delta = float("{0:.6f}".format(self.err_d['o_pm_delta'][1]))
    #         if type(self.fitted_d['ra']) is str:
    #             chi_ra = 'empty'
    #         else:
    #             chi_ra = float("{0:.6f}".format(self.fitted_d['ra']))
    #         if type(self.fitted_d['dec']) is str:
    #             chi_dec = 'empty'
    #         else:
    #             chi_dec = float("{0:.6f}".format(self.fitted_d['dec']))
    #
    #         table_ = [['', 'catalog values'],
    #                   ['cat_pm', i_pm],
    #                   ['cat_pm_alpha', i_pm_alpha],
    #                   ['cat_pm_delta', i_pm_delta],
    #                   ['', 'extracted values'],
    #                   ['ext_pm', o_pm],
    #                   ['ext_pm_alpha', o_pm_alpha],
    #                   ['ext_pm_delta', o_pm_delta],
    #                   ['', 'chi squared'],
    #                   ['chi_squared_ra', chi_ra],
    #                   ['chi_squared_dec', chi_dec]]
    #
    #         # ax_2.axis('tight')
    #         ax_2.axis('off')
    #         # ax_2.axis('on')
    #         data_table = ax_2.table(cellText=table_, colLabels=None,
    #                                 loc='center')
    #         data_table.set_fontsize(14)
    #
    #         # fig.tight_layout()
    #
    #         # Saves the current figure in pdf file
    #         pdf.savefig()  # saves current figure
    #         plt.close(fig)
    #
    #         #
    #         # IMAGES
    #         #
    #
    #         # Gets limits
    #         alpha_center = 0
    #         delta_center = 0
    #         if len(self.err_d['i_alpha']) == 2:
    #             alpha_sum = self.err_d['i_alpha'][0] + self.err_d['i_alpha'][1]
    #             alpha_center = alpha_sum / 2
    #             delta_sum = self.err_d['i_delta'][0] + self.err_d['i_delta'][1]
    #             delta_center = delta_sum / 2
    #         elif len(self.err_d['i_alpha']) == 3:
    #             alpha_center = self.err_d['i_alpha'][1]
    #             delta_center = self.err_d['i_delta'][1]
    #         elif len(self.err_d['i_alpha']) == 4:
    #             alpha_sum = self.err_d['i_alpha'][1] + self.err_d['i_alpha'][2]
    #             alpha_center = alpha_sum / 2
    #             delta_sum = self.err_d['i_delta'][1] + self.err_d['i_delta'][2]
    #             delta_center = delta_sum / 2
    #
    #         if alpha_center is 0 or delta_center is 0:
    #             print('ERROR')
    #             print("self.err_d['i_alpha'] {}".format(self.err_d['i_alpha']))
    #             print("self.err_d['i_delta'] {}".format(self.err_d['i_delta']))
    #
    #             raise Exception
    #
    #         for idx in range(0, len(self.err_d['i_alpha']), 1):
    #             # Get regions for desired fits file
    #             regions_file = get_regions(self.fits_files[idx],
    #                                        self.err_d['sex_cf'][0],
    #                                        self.prfs_d, self.mag)
    #
    #             # Creates image
    #             dither = self.fits_files[idx][-6:-5]
    #             i_regs = '{}/dither_{}.reg'.format(self.prfs_d['dithers_out'],
    #                                                dither)
    #             if idx == 0:
    #                 p_alpha = 0
    #                 p_delta = 0
    #             else:
    #                 p_alpha = self.err_d['i_alpha'][idx - 1]
    #                 p_delta = self.err_d['i_delta'][idx - 1]
    #
    #             # Zoomed image
    #             zoom = True
    #             img = image(i_regs, self.fits_files[idx], regions_file,
    #                         alpha_center, delta_center, idx,
    #                         float(self.err_d['i_pm'][0]), p_alpha, p_delta,
    #                         zoom)
    #
    #             img = img[1:-50, :]  # Removes ds9's legend
    #
    #             fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
    #             ax_1 = fig.add_subplot(1, 1, 1)
    #             ax_1.set_title('Dither {} - Enlarged view '.format(dither))
    #
    #             plt.imshow(img)
    #
    #             labels = [item.get_text() for item in ax_1.get_xticklabels()]
    #
    #             empty_string_labels = [''] * len(labels)
    #             ax_1.set_xticklabels(empty_string_labels)
    #
    #             labels = [item.get_text() for item in ax_1.get_yticklabels()]
    #
    #             empty_string_labels = [''] * len(labels)
    #             ax_1.set_yticklabels(empty_string_labels)
    #
    #             ax_1.set_xlabel('right ascension')
    #             ax_1.set_ylabel('declination')
    #
    #             pdf.savefig()
    #             plt.close(fig)
    #
    #             # Wide image
    #             zoom = False
    #             img = image(i_regs, self.fits_files[idx], regions_file,
    #                         alpha_center, delta_center, idx,
    #                         float(self.err_d['i_pm'][0]), p_alpha, p_delta,
    #                         zoom)
    #
    #             img = img[1:-50, :]  # Removes ds9's legend
    #
    #             fig = plt.figure(figsize=self.plot_size, dpi=self.plot_dpi)
    #             ax_1 = fig.add_subplot(1, 1, 1)
    #             ax_1.set_title('Dither {} - Wide view '.format(dither))
    #
    #             plt.imshow(img)
    #
    #             labels = [item.get_text() for item in ax_1.get_xticklabels()]
    #
    #             empty_string_labels = [''] * len(labels)
    #             ax_1.set_xticklabels(empty_string_labels)
    #
    #             labels = [item.get_text() for item in ax_1.get_yticklabels()]
    #
    #             empty_string_labels = [''] * len(labels)
    #             ax_1.set_yticklabels(empty_string_labels)
    #
    #             ax_1.set_xlabel('right ascension')
    #             ax_1.set_ylabel('declination')
    #
    #             pdf.savefig()
    #             plt.close(fig)
    #
    #
