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
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from misc import significant_l
import numpy as np


__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def create_labels(datetimes):
    """

    :param datetimes:
    :return: labels
    """
    labels = []
    for datetime_ in datetimes:
        label_ = '{}.{}'.format(datetime_.second,
                                datetime_.microsecond)
        labels.append(label_)

    return labels


def round_number(number):
    """
    TODO Improve description


    :param number:
    :return:
    """
    number_t = significant_l(number)

    # Rounds the float value of difference
    first_digit = '%.2E' % Decimal(number - number_t)

    if first_digit[0] == '-':
        first_digit = int(first_digit[1])
    else:
        first_digit = int(first_digit[0])

    # Redondea al alza el ultimo digio
    if first_digit > 5:
        last_digit = str(number_t)[-1]
        last_digit = int(last_digit)
        last_digit += 1
        number_t = list(str(number_t))
        number_t[-1] = last_digit
        number_t = [str(i) for i in number_t]
        number_t = ''.join(number_t)
        number_t = float(number_t)
    else:
        pass  # do nothing

    return number_t


def create_x_ticks(epoch_seconds):
    """

    :param epoch_seconds:
    :return: x_ticks
    """
    x_difference = float(max(epoch_seconds)) - float(min(epoch_seconds))
    x_major_stp = (x_difference / 3)
    x_minor_stp = (x_difference / 6)

    # X-SCALE
    # Major steps
    x_major_stps = np.arange((epoch_seconds[0] - x_major_stp * 1),
                             (epoch_seconds[0] + x_major_stp * 4), x_major_stp)

    # Minor steps
    x_minor_stps = np.arange((epoch_seconds[0] - x_minor_stp * 2),
                             (epoch_seconds[0] + x_minor_stp * 8),
                             x_minor_stp)

    x_ticks = {'major_t': x_major_stps, 'minor_t': x_minor_stps,
               'major_s': x_major_stp, 'minor_s': x_minor_stp}

    return x_ticks


def create_alpha_y_ticks(alpha_seconds, hour, minute):
    """

    :param alpha_seconds:
    :param hour:
    :param minute:
    :return:
    """
    divisions = 4  #

    # Gets the major step between ticks thought the difference
    # between the maximum and the minimum value of alpha
    difference = float(max(alpha_seconds)) - float(min(alpha_seconds))
    major_stp = (difference / divisions)
    #
    major_stp = round_number(major_stp)

    # In order to avoid too crowned lines divisions will be less if
    # steps are little
    if float(major_stp) == 0.0001:
        divisions = 2
    elif float(major_stp) < 0.0001:
        divisions = 2
    else:
        divisions = 4
    minor_stp = (float(major_stp) / divisions)

    # Gets maximum decimal position of major step
    decimals = int(str(major_stp)[::-1].find('.'))

    # todo, document!
    major_margin = major_stp * (divisions / 2)
    initial_major_stp = round(min(alpha_seconds), decimals) - major_margin
    final_major_stp = round(max(alpha_seconds), decimals) + major_margin
    major_stps = np.arange(initial_major_stp, final_major_stp, major_stp)
    # todo, document!
    minor_margin = minor_stp * divisions
    initial_minor_stp = round(min(alpha_seconds), decimals) - minor_margin
    final_minor_stp = round(max(alpha_seconds), decimals) + minor_margin
    minor_stps = np.arange(initial_minor_stp, final_minor_stp, minor_stp)

    # Formats major steps for be readable by matplotlib axis
    alpha_major_stps = []
    for idx_stp, stp_ in enumerate(major_stps):
        if float('{:.6f}'.format(stp_)) > 60.000000:
            minute_ = minute + 1
            # todo decimal aproximation should b automatically, not hardcoded
            # from decimal import Decimal
            # d = Decimal(str(major_stp))
            # u_decimals_pos = d.as_tuple().exponent
            # decimals_pos = u_decimals_pos + (u_decimals_pos * -2)
            stp_ = '{:.6f}'.format(stp_)
            second_ = float(stp_) - 60.0
            second_ = '{:.6f}'.format(second_)
            alpha_major_stps.append('{}:{}:{}'.format(hour, minute_, second_))
        elif float('{:.6f}'.format(stp_)) == 60.000000:
            minute_ = minute + 1
            # todo decimal aproximation should b automatically, not hardcoded
            # from decimal import Decimal
            # d = Decimal(str(major_stp))
            # u_decimals_pos = d.as_tuple().exponent
            # decimals_pos = u_decimals_pos + (u_decimals_pos * -2)
            stp_ = '{:.6f}'.format(stp_)
            second_ = float(stp_) - 60.0
            second_ = '{:.6f}'.format(second_)
            alpha_major_stps.append('{}:{}:{}'.format(hour, minute_, second_))
        elif float('{:.6f}'.format(stp_)) < 0.000000:
            # todo implement backward number
            pass
        elif float('{:.6f}'.format(stp_)) < 60.000000:
            minute_ = minute
            second_ = stp_
            alpha_major_stps.append('{}:{}:{}'.format(hour, minute_, second_))
    # Formats minor steps
    alpha_minor_stps = []
    for idx_stp, stp_ in enumerate(minor_stps):
        # According step value... todo
        if stp_ > 60.0:
            minute_ = minute + 1
            stp_ = '{:.6f}'.format(stp_)
            second_ = float(stp_) - 60.0
            second_ = '{:.6f}'.format(second_)
            alpha_minor_stps.append('{}:{}:{}'.format(hour, minute_, second_))
        elif float('{:.1f}'.format(stp_)) == 60.0:
            minute_ = minute + 1
            second_ = 0.0
            second_ = '{:.6f}'.format(second_)
            alpha_minor_stps.append('{}:{}:{}'.format(hour, minute_, second_))
        else:
            minute_ = minute
            second_ = '{:.6f}'.format(stp_)
            alpha_minor_stps.append('{}:{}:{}'.format(hour, minute_, second_))

    # Creates a list of datatime objects
    # Sometimes, seconds are greater than 59, this values should be
    # filter in a better way
    major_t = [datetime.strptime(t, "%H:%M:%S.%f") for t in alpha_major_stps]
    minor_t = [datetime.strptime(t, "%H:%M:%S.%f") for t in alpha_minor_stps]
    y_ticks = {'major_t': major_t, 'minor_t': minor_t,
               'major_s': major_stp, 'minor_s': minor_stp}

    return y_ticks


def create_delta_y_ticks(delta_seconds):
    """
    step = stp

    :param delta_seconds:
    :return: y_ticks
    """
    divisions = 4  #

    # Gets the major step between ticks thought the difference
    # between the maximum and the minimum value of alpha
    difference = float(max(delta_seconds)) - float(min(delta_seconds))
    major_stp = (difference / divisions)
    #
    major_stp = float(round_number(major_stp))
    minor_stp = (major_stp / 4)

    # Gets maximum decimal position of major step
    decimals = int(str(major_stp)[::-1].find('.'))

    # Major step list starts two times before and end two times after
    # known values
    major_stps = np.arange(round(min(delta_seconds), decimals) - major_stp * 2,
                           round(max(delta_seconds), decimals) + major_stp * 2,
                           major_stp)
    # Minor step list starts eight times before and fend eight times
    # after know values
    minor_stps = np.arange(round(min(delta_seconds), decimals) - minor_stp * 8,
                           round(max(delta_seconds), decimals) + minor_stp * 8,
                           minor_stp)

    y_ticks = {'major_t': major_stps, 'minor_t': minor_stps,
               'major_s': major_stp, 'minor_s': minor_stp}

    return y_ticks


class PlotConfidence:

    def __init__(self, output_path, source, pm, mode, fitted_d, tmp_d):
        """

        :param output_path:
        :param source:
        :param pm:
        :param mode:
        :param fitted_d:
        :param tmp_d: a dictionary ['alpha', 'delta',
                                    'error_a', 'error_b', 'epoch']
        """
        self.output_path = output_path
        self.source = source
        self.pm = pm
        self.mode = mode  # Can be 'i' or 'o'
        self.fitted_d = fitted_d
        self.tmp_d = tmp_d

        plot_size = [16.53, 11.69]
        plot_dpi = 100

        self.plot(plot_size, plot_dpi)

    def plot(self, plot_size, plot_dpi):
        """

        :param plot_size:
        :param plot_dpi:
        :return:
        """
        pdf_name = '{}/{}_{}_{}.pdf'.format(self.output_path, self.source,
                                            self.pm, self.mode)
        with PdfPages(pdf_name) as pdf:
            # ALPHA PARAMETERS
            fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax = fig.add_subplot(1, 1, 1)
            chi_squared = float("{0:.6f}".format(self.fitted_d['ra']))
            chi_squared_title = 'chi_squared: {}'.format(chi_squared)

            pm_alpha_cat = float("{0:.6f}".format(self.tmp_d['i_pm_alpha'][0]))
            pm_alpha_cat_str = 'catalog {}'.format(pm_alpha_cat)

            pm_alpha_ext = float("{0:.6f}".format(self.tmp_d['o_pm_alpha'][0]))
            pm_alpha_ext_str = 'extracted {}'.format(pm_alpha_ext)

            o_pm_alpha = float(self.tmp_d['o_pm_alpha'][0])
            o_pm_alpha_err = float(self.tmp_d['o_pm_alpha_err'][0])
            o_pm_alpha_sn = o_pm_alpha / o_pm_alpha_err
            o_pm_alpha_sn = float("{0:.6f}".format(o_pm_alpha_sn))
            pm_alpha_sn_str = 'SN {}'.format(o_pm_alpha_sn)
            pm_alpha_ext_err = self.tmp_d['o_pm_alpha_err'][0]
            pm_alpha_ext_err = float("{0:.6f}".format(pm_alpha_ext_err))
            pm_alpha_ext_err_str = 'pm_error {}'.format(pm_alpha_ext_err)

            pm_alpha_comp = '{} / {}'.format(pm_alpha_ext_str,
                                             pm_alpha_cat_str)
            pm_error_comp = '{} / {}'.format(pm_alpha_sn_str,
                                             pm_alpha_ext_err_str)

            ax.set_title('{}\n{}\n{}'.format(chi_squared_title, pm_alpha_comp,
                                             pm_error_comp))

            alpha_tmpstmp = []
            alpha_seconds = []
            tmp_hour = []
            tmp_minute = []
            for alpha_ in self.tmp_d['alpha']:
                a = Angle(alpha_, units.degree)
                dms = a.dms
                degree = int(dms[0])
                tmp_hour.append(degree)
                minute = int(dms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(dms[2]))
                alpha_tmpstmp.append('{}:{}:{}'.format(degree, minute, second))
                alpha_seconds.append('{}'.format(second))

            # X-AXIS Shared between input and output catalogs
            epoch_seconds = []
            for epoch_ in self.tmp_d['epoch']:
                frac, whole = modf(epoch_)
                seconds_ = frac * 365.25 * 24 * 60
                epoch_seconds.append(seconds_)

            # Creates x ticks (major and minor ones)
            x_ticks = create_x_ticks(epoch_seconds)

            # Y-AXIS
            myformat = mdates.DateFormatter('%S.%f')
            ax.yaxis.set_major_formatter(myformat)

            # format elements in alpha_seconds list to floats
            alpha_seconds = [float(i) for i in alpha_seconds]  # needed?
            alpha_seconds = [float("{0:.6f}".format(i)) for i in
                             alpha_seconds]

            # Check if all hours/minutes are same
            if len(list(set(tmp_hour))) != 1:
                raise Exception
            if len(list(set(tmp_minute))) != 1:
                raise Exception

            # Creates y ticks (major and minor ones)
            # y_ticks = create_alpha_y_ticks(alpha_seconds, tmp_hour[0],
            #                                tmp_minute[0])
            # y_labels = create_labels(y_ticks['minor_t'])

            y_ticks = create_delta_y_ticks(alpha_seconds)

            # x-ticks assignation
            ax.set_xticks(x_ticks['major_t'], minor=False)
            ax.set_xticks(x_ticks['minor_t'], minor=True)
            # x-ticks labels

            # y-ticks assignation
            ax.set_yticks(y_ticks['major_t'], minor=False)
            ax.set_yticks(y_ticks['minor_t'], minor=True)
            # y-ticks labels
            empty_string_labels = [''] * len(y_ticks['major_t'])
            ax.set_yticklabels(empty_string_labels, minor=False)
            ax.set_yticklabels(y_ticks['minor_t'], minor=True)

            # Formats grids
            ax.grid(b=True, which='major', linestyle='-', linewidth=2)
            ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

            # Annotations
            for idx_alpha in range(0, len(alpha_tmpstmp), 1):
                x = epoch_seconds[idx_alpha]  # x position
                y = alpha_seconds[idx_alpha]  # y position
                # variables for position and error representation
                alpha_position = alpha_tmpstmp[idx_alpha]
                alpha_str = '   alpha {}'.format(alpha_position)
                error = self.tmp_d['error_b'][idx_alpha] * 3600  # error-y
                error_fmt = float("{0:.6f}".format(error))  # error rounded
                error_str = '   error {}'.format(error_fmt)
                # annotation
                ax.annotate('{}\n{}'.format(alpha_str, error_str),
                            xy=(x, y), textcoords='data', fontsize=13)

            for idx_alpha in range(0, len(alpha_tmpstmp), 1):
                x = epoch_seconds[idx_alpha]  # x position
                y = alpha_seconds[idx_alpha]  # y position
                error = float(self.tmp_d['error_b'][idx_alpha] * 3600)
                ax.errorbar(x, y, yerr=error, fmt='o',
                            ecolor='g', capthick=2, elinewidth=4)

            # Axis labels creation
            # x-axis
            x_label = 'EPOCH (seconds)'
            ax.set_xlabel(x_label)
            # y-axis
            y_label_ra = 'Right ascension (")\n'
            major_s = y_ticks['major_s']
            y_label_major_step = 'major step size {}"\n'.format(major_s)
            minor_s = y_ticks['minor_s']
            y_label_minor_step = 'minor step size {}"'.format(minor_s)
            ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                          y_label_minor_step))

            # Plots data
            ax.plot(epoch_seconds, alpha_seconds, 'bs', markersize=6,
                    label='extracted position')
            # In order to avoid scientific notation plot should be redrawn
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)

            # Plots a legend
            plt.legend(loc=0, ncol=2, borderaxespad=0.)
            plt.draw()

            # Saves the current figure in pdf file
            pdf.savefig()  # saves current figure
            # plt.clf()  # clear current figure
            plt.close(fig)

            #
            # DELTA PARAMETERS
            # Another figure is created in order to plot delta output.
            # X-axis values are shared between both figures but Y-axis for
            # declination's output is showed in seconds (floats) instead
            # datetime objects.
            #
            fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax = fig.add_subplot(1, 1, 1)
            chi_squared = float("{0:.6f}".format(self.fitted_d['dec']))
            chi_squared_title = 'chi_squared: {}'.format(chi_squared)

            pm_delta_cat = float("{0:.6f}".format(self.tmp_d['i_pm_delta'][0]))
            pm_delta_cat_str = 'catalog {}'.format(pm_delta_cat)

            pm_delta_ext = float("{0:.6f}".format(self.tmp_d['o_pm_delta'][0]))
            pm_delta_ext_str = 'extracted {}'.format(pm_delta_ext)

            o_pm_delta = float(self.tmp_d['o_pm_delta'][0])
            o_pm_delta_err = float(self.tmp_d['o_pm_delta_err'][0])
            o_pm_delta_sn = o_pm_delta / o_pm_delta_err
            o_pm_delta_sn = float("{0:.6f}".format(o_pm_delta_sn))
            pm_delta_sn_str = 'SN {}'.format(o_pm_delta_sn)
            pm_delta_ext_err = self.tmp_d['o_pm_delta_err'][0]
            pm_delta_ext_err = float("{0:.6f}".format(pm_delta_ext_err))
            pm_delta_ext_err_str = 'pm_error {}'.format(pm_delta_ext_err)

            pm_delta_comp = '{} / {}'.format(pm_delta_ext_str,
                                             pm_delta_cat_str)
            pm_error_comp = '{} / {}'.format(pm_delta_sn_str,
                                             pm_delta_ext_err_str)

            ax.set_title('{}\n{}\n{}'.format(chi_squared_title, pm_delta_comp,
                                             pm_error_comp))

            delta_tmpstmp = []
            delta_seconds = []
            tmp_degree = []
            tmp_minute = []
            for delta_ in self.tmp_d['delta']:
                d = Angle(delta_, units.degree)
                dms = d.dms
                degree = int(dms[0])
                tmp_degree.append(degree)
                minute = int(dms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(dms[2]))
                delta_tmpstmp.append('{}:{}.{}'.format(degree, minute, second))
                delta_seconds.append('{}'.format(second))

            ax.set_xlabel('EPOCH')

            # Y-SCALE
            # format elements in alpha_seconds list to floats
            delta_seconds = [float(i) for i in delta_seconds]
            delta_seconds = [float("{0:.6f}".format(i)) for i in delta_seconds]

            # Check if all hours/minutes are same
            if len(list(set(tmp_degree))) != 1:
                raise Exception
            if len(list(set(tmp_minute))) != 1:
                raise Exception

            y_ticks = create_delta_y_ticks(delta_seconds)

            # x-ticks assignation
            ax.set_xticks(x_ticks['major_t'], minor=False)
            ax.set_xticks(x_ticks['minor_t'], minor=True)
            # x-ticks labels

            # y-ticks assignation
            ax.set_yticks(y_ticks['major_t'], minor=False)
            ax.set_yticks(y_ticks['minor_t'], minor=True)
            # y-ticks labels
            empty_string_labels = [''] * len(y_ticks['major_t'])
            ax.set_yticklabels(empty_string_labels, minor=False)
            # Converts floats into strings
            # y_minor_ticks = [str(i) for i in y_minor_ticks]
            ax.set_yticklabels(y_ticks['minor_t'], minor=True)

            # Formats grids
            ax.grid(b=True, which='major', linestyle='-', linewidth=2)
            ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

            for idx_delta in range(0, len(delta_tmpstmp), 1):
                x = epoch_seconds[idx_delta]  # x position
                y = delta_seconds[idx_delta]  # y position
                # variables for position and error representation
                delta_position = delta_tmpstmp[idx_delta]
                delta_str = '   delta {}'.format(delta_position)
                error = self.tmp_d['error_b'][idx_delta] * 3600  # error-y
                error_fmt = float("{0:.6f}".format(error))  # error rounded
                error_str = '   error {}'.format(error_fmt)
                # annotation
                ax.annotate('{}\n{}'.format(delta_str, error_str),
                            xy=(x, y), textcoords='data', fontsize=13)

            for idx_delta in range(0, len(delta_tmpstmp), 1):
                x = epoch_seconds[idx_delta]  # x position
                y = delta_seconds[idx_delta]  # y position
                error = float(self.tmp_d['error_b'][idx_delta] * 3600)
                ax.errorbar(x, y, yerr=error, fmt='o',
                            ecolor='g', capthick=2, elinewidth=4)

            # Label creation
            y_label_ra = 'Declination (")\n'
            major_s = y_ticks['major_s']
            y_label_major_step = 'major step size {}"\n'.format(major_s)
            minor_s = y_ticks['minor_s']
            y_label_minor_step = 'minor step size {}"'.format(minor_s)
            ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                          y_label_minor_step))

            ax.plot(epoch_seconds, delta_seconds, 'bs', markersize=6,
                    label='extracted position')

            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            plt.legend(loc=0, ncol=2, borderaxespad=0.)
            plt.draw()

            pdf.savefig()  # saves current figure
            # plt.clf()  # clear current figure
            plt.close(fig)

        return True


class PlotBothConfidence:

    def __init__(self, output_path, source, pm, mode, fitted_d, tmp_d):
        """

        :param output_path:
        :param source:
        :param pm:
        :param mode:
        :param fitted_d:
        :param tmp_d: a dictionary ['i_alpha', 'i_delta', 'o_alpha', 'o_delta',
                                    'error_a', 'error_b', 'epoch']
        """
        self.output_path = output_path
        self.source = source
        self.pm = pm
        self.mode = mode  # Can be 'i' or 'o'
        self.fitted_d = fitted_d
        self.tmp_d = tmp_d

        plot_size = [16.53, 11.69]
        plot_dpi = 100

        self.plot(plot_size, plot_dpi)

    def plot(self, plot_size, plot_dpi):
        """

        :param plot_size:
        :param plot_dpi:
        :return:
        """
        pdf_name = '{}/{}_{}_{}.pdf'.format(self.output_path, self.source,
                                            self.pm, self.mode)
        with PdfPages(pdf_name) as pdf:
            # ALPHA PARAMETERS
            fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax = fig.add_subplot(1, 1, 1)
            chi_squared = float("{0:.6f}".format(self.fitted_d['ra']))
            chi_squared_title = 'chi_squared: {}'.format(chi_squared)

            pm_alpha_cat = float("{0:.6f}".format(self.tmp_d['i_pm_alpha'][0]))
            pm_alpha_cat_str = 'catalog {}'.format(pm_alpha_cat)

            pm_alpha_ext = float("{0:.6f}".format(self.tmp_d['o_pm_alpha'][0]))
            pm_alpha_ext_str = 'extracted {}'.format(pm_alpha_ext)

            o_pm_alpha = float(self.tmp_d['o_pm_alpha'][0])
            o_pm_alpha_err = float(self.tmp_d['o_pm_alpha_err'][0])
            o_pm_alpha_sn = o_pm_alpha / o_pm_alpha_err
            o_pm_alpha_sn = float("{0:.6f}".format(o_pm_alpha_sn))
            pm_alpha_sn_str = 'SN {}'.format(o_pm_alpha_sn)
            pm_alpha_ext_err = self.tmp_d['o_pm_alpha_err'][0]
            pm_alpha_ext_err = float("{0:.6f}".format(pm_alpha_ext_err))
            pm_alpha_ext_err_str = 'pm_error {}'.format(pm_alpha_ext_err)

            pm_alpha_comp = '{} / {}'.format(pm_alpha_ext_str,
                                             pm_alpha_cat_str)
            pm_error_comp = '{} / {}'.format(pm_alpha_sn_str,
                                             pm_alpha_ext_err_str)

            ax.set_title('{}\n{}\n{}'.format(chi_squared_title, pm_alpha_comp,
                                             pm_error_comp))

            # input_parameters
            plot_flag = True
            i_a_tmpstmp = []
            i_alpha_seconds = []
            i_tmp_hour = []
            i_tmp_minute = []
            for alpha_ in self.tmp_d['i_alpha']:  # fixme change to dms
                a = Angle(alpha_, units.degree)
                hms = a.hms
                hour = int(hms[0])
                i_tmp_hour.append(hour)
                minute = int(hms[1])
                i_tmp_minute.append(minute)
                second = float("{0:.6f}".format(hms[2]))
                i_a_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
                i_alpha_seconds.append('{}'.format(second))

            i_datetimes = []
            for idx, i_datetime_ in enumerate(i_a_tmpstmp):
                i_datetime_tmp = datetime.strptime(i_datetime_, "%H:%M:%S.%f")
                i_datetimes.append(i_datetime_tmp)

            # output parameters
            o_a_tmpstmp = []
            o_alpha_seconds = []
            o_tmp_hour = []
            o_tmp_minute = []
            for alpha_ in self.tmp_d['o_alpha']:
                a = Angle(alpha_, units.degree)
                hms = a.hms
                hour = int(hms[0])
                o_tmp_hour.append(hour)
                minute = int(hms[1])
                o_tmp_minute.append(minute)
                second = float("{0:.6f}".format(hms[2]))
                o_a_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
                o_alpha_seconds.append('{}'.format(second))

            o_datetimes = []
            for idx, o_datetime_ in enumerate(o_a_tmpstmp):
                o_datetime_tmp = datetime.strptime(o_datetime_, "%H:%M:%S.%f")
                o_datetimes.append(o_datetime_tmp)

            # X-AXIS Shared between input and output catalogs
            epoch_seconds = []
            for epoch_ in self.tmp_d['epoch']:
                frac, whole = modf(epoch_)
                seconds_ = frac * 365.25 * 24 * 60
                epoch_seconds.append(seconds_)

            # Creates x ticks (major and minor ones)
            x_ticks = create_x_ticks(epoch_seconds)

            # todo reformat!
            # Y-AXIS
            myformat = mdates.DateFormatter('%S.%f')
            ax.yaxis.set_major_formatter(myformat)

            # format elements in alpha_seconds list to floats
            i_alpha_seconds = [float(i) for i in i_alpha_seconds]  # needed?
            i_alpha_seconds = [float("{0:.6f}".format(i)) for i in
                               i_alpha_seconds]
            o_alpha_seconds = [float(i) for i in o_alpha_seconds]  # needed?
            o_alpha_seconds = [float("{0:.6f}".format(i)) for i in
                               o_alpha_seconds]

            # Check if all hours/minutes are same
            if len(list(set(i_tmp_hour))) != 1:
                plot_flag = False
            if len(list(set(o_tmp_hour))) != 1:
                plot_flag = False
            if len(list(set(i_tmp_minute))) != 1:
                plot_flag = False
            if len(list(set(o_tmp_minute))) != 1:
                plot_flag = False

            if plot_flag:
                alpha_seconds = i_alpha_seconds + o_alpha_seconds
                hour = i_tmp_hour[0]
                minute = i_tmp_minute[0]

                # Creates y ticks (major and minor ones)
                y_ticks = create_alpha_y_ticks(alpha_seconds, hour, minute)

                y_labels = create_labels(y_ticks['minor_t'])

                # x-ticks assignation
                ax.set_xticks(x_ticks['major_t'], minor=False)
                ax.set_xticks(x_ticks['minor_t'], minor=True)
                # x-ticks labels

                # y-ticks assignation
                ax.set_yticks(y_ticks['major_t'], minor=False)
                ax.set_yticks(y_ticks['minor_t'], minor=True)
                # y-ticks labels
                empty_string_labels = [''] * len(y_ticks['major_t'])
                ax.set_yticklabels(empty_string_labels, minor=False)
                ax.set_yticklabels(y_labels, minor=True)

                # Formats grids
                ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                # Annotations
                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    # Format alpha
                    hour = datetime_.hour
                    minute = datetime_.minute
                    second = datetime_.second
                    msecond = datetime_.microsecond
                    alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
                                                              second, msecond)

                    # Annotate position and error associated
                    ax.annotate('{}'.format(alpha_str),
                                xy=(epoch_seconds[idx_datetime_], datetime_),
                                textcoords='data', fontsize=13)

                # Annotations
                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    # Format alpha
                    hour = datetime_.hour
                    minute = datetime_.minute
                    second = datetime_.second
                    msecond = datetime_.microsecond
                    alpha_str = '   alpha {}:{}:{}:{}'.format(hour, minute,
                                                              second, msecond)
                    # Format error
                    error_seconds = self.tmp_d['error_a'][idx_datetime_] * 3600
                    error = float("{0:.6f}".format(error_seconds))
                    error_str = '   error  {}"'.format(error)
                    # Annotate position and error associated
                    ax.annotate('{}\n{}'.format(alpha_str, error_str),
                                xy=(epoch_seconds[idx_datetime_], datetime_),
                                textcoords='data', fontsize=13)
                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    error_seconds = self.tmp_d['error_a'][idx_datetime_] * 3600
                    errorbar_size = timedelta(0, error_seconds)
                    ax.errorbar(epoch_seconds[idx_datetime_], datetime_,
                                yerr=errorbar_size,
                                ecolor='g', capthick=2, elinewidth=4)

                # Axis labels creation
                # x-axis
                x_label = 'EPOCH (seconds)'
                ax.set_xlabel(x_label)
                # y-axis
                y_label_ra = 'Right ascension (")\n'
                major_s = y_ticks['major_s']
                y_label_major_step = 'major step size {}"\n'.format(major_s)
                minor_s = y_ticks['minor_s']
                y_label_minor_step = 'minor step size {}"'.format(minor_s)
                ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                              y_label_minor_step))

                # Plots data
                ax.plot(epoch_seconds, i_datetimes, 'bs', markersize=6,
                        label='catalog position')
                ax.plot(epoch_seconds, o_datetimes, 'rs', markersize=6,
                        label='extracted position')
                ax.plot(epoch_seconds, i_datetimes, linewidth=1)
                # In order to avoid scientific notation plot should be redrawn
                ax = plt.gca()
                ax.get_xaxis().get_major_formatter().set_useOffset(False)

                # Plots a legend
                plt.legend(loc=0, ncol=2, borderaxespad=0.)
                plt.draw()

                # Saves the current figure in pdf file
                pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure
            plt.close(fig)

            #
            # DELTA PARAMETERS
            # Another figure is created in order to plot delta output.
            # X-axis values are shared between both figures but Y-axis for
            # declination's output is showed in seconds (floats) instead
            # datetime objects.
            #
            fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax = fig.add_subplot(1, 1, 1)
            chi_squared = float("{0:.6f}".format(self.fitted_d['dec']))
            chi_squared_title = 'chi_squared: {}'.format(chi_squared)

            pm_delta_cat = float("{0:.6f}".format(self.tmp_d['i_pm_delta'][0]))
            pm_delta_cat_str = 'catalog {}'.format(pm_delta_cat)

            pm_delta_ext = float("{0:.6f}".format(self.tmp_d['o_pm_delta'][0]))
            pm_delta_ext_str = 'extracted {}'.format(pm_delta_ext)

            o_pm_delta = float(self.tmp_d['o_pm_delta'][0])
            o_pm_delta_err = float(self.tmp_d['o_pm_delta_err'][0])
            o_pm_delta_sn = o_pm_delta / o_pm_delta_err
            o_pm_delta_sn = float("{0:.6f}".format(o_pm_delta_sn))
            pm_delta_sn_str = 'SN {}'.format(o_pm_delta_sn)
            pm_delta_ext_err = self.tmp_d['o_pm_delta_err'][0]
            pm_delta_ext_err = float("{0:.6f}".format(pm_delta_ext_err))
            pm_delta_ext_err_str = 'pm_error {}'.format(pm_delta_ext_err)

            pm_delta_comp = '{} / {}'.format(pm_delta_ext_str,
                                             pm_delta_cat_str)
            pm_error_comp = '{} / {}'.format(pm_delta_sn_str,
                                             pm_delta_ext_err_str)

            ax.set_title('{}\n{}\n{}'.format(chi_squared_title, pm_delta_comp,
                                             pm_error_comp))

            plot_flag = True
            i_d_tmpstmp = []
            i_delta_seconds = []
            i_tmp_degree = []
            i_tmp_minute = []
            for delta_ in self.tmp_d['i_delta']:
                d = Angle(delta_, units.degree)
                dms = d.dms
                degree = int(dms[0])
                i_tmp_degree.append(degree)
                minute = int(dms[1])
                i_tmp_minute.append(minute)
                second = float("{0:.6f}".format(dms[2]))
                i_d_tmpstmp.append('{}:{}.{}'.format(degree, minute, second))
                i_delta_seconds.append('{}'.format(second))

            o_d_tmpstmp = []
            o_delta_seconds = []
            o_tmp_degree = []
            o_tmp_minute = []
            for delta_ in self.tmp_d['o_delta']:
                d = Angle(delta_, units.degree)
                dms = d.dms
                degree = int(dms[0])
                o_tmp_degree.append(degree)
                minute = int(dms[1])
                o_tmp_minute.append(minute)
                second = float("{0:.6f}".format(dms[2]))
                o_d_tmpstmp.append('{}:{}.{}'.format(degree, minute, second))
                o_delta_seconds.append('{}'.format(second))

            ax.set_xlabel('EPOCH')

            # Y-SCALE
            # format elements in alpha_seconds list to floats
            i_delta_seconds_ = []
            for idx_, i_delta_second_ in enumerate(i_delta_seconds):
                i_delta_tmp = float("{0:.6f}".format(float(i_delta_second_)))
                i_delta_seconds_.append(i_delta_tmp)
            i_delta_seconds = i_delta_seconds_  # it's really needed?

            o_delta_seconds_ = []
            for idx_, o_delta_second_ in enumerate(o_delta_seconds):
                o_delta_tmp = float("{0:.6f}".format(float(o_delta_second_)))
                o_delta_seconds_.append(o_delta_tmp)
            o_delta_seconds = o_delta_seconds_  # it's really needed?

            delta_seconds = i_delta_seconds + o_delta_seconds

            # Check if all hours/minutes are same
            if len(list(set(i_tmp_degree))) != 1:
                plot_flag = False
            if len(list(set(i_tmp_minute))) != 1:
                plot_flag = False
            if len(list(set(o_tmp_degree))) != 1:
                plot_flag = False
            if len(list(set(o_tmp_minute))) != 1:
                plot_flag = False

            if plot_flag:
                y_ticks = create_delta_y_ticks(delta_seconds)

                # x-ticks assignation
                ax.set_xticks(x_ticks['major_t'], minor=False)
                ax.set_xticks(x_ticks['minor_t'], minor=True)
                # x-ticks labels

                # y-ticks assignation
                ax.set_yticks(y_ticks['major_t'], minor=False)
                ax.set_yticks(y_ticks['minor_t'], minor=True)
                # y-ticks labels
                empty_string_labels = [''] * len(y_ticks['major_t'])
                ax.set_yticklabels(empty_string_labels, minor=False)
                # Converts floats into strings
                # y_minor_ticks = [str(i) for i in y_minor_ticks]
                ax.set_yticklabels(y_ticks['minor_t'], minor=True)

                # Formats grids
                ax.grid(b=True, which='major', linestyle='-', linewidth=2)
                ax.grid(b=True, which='minor', linestyle='--', linewidth=1)

                for idx_datetime_, datetime_ in enumerate(i_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = i_delta_seconds[idx_datetime_]  # y position
                    # variables for position and error representation
                    delta_position = i_d_tmpstmp[idx_datetime_]
                    delta_str = '   delta {}'.format(delta_position)
                    # annotation
                    ax.annotate('{}'.format(delta_str),
                                xy=(x, y), textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    # variables for position and error representation
                    delta_position = o_d_tmpstmp[idx_datetime_]
                    delta_str = '   delta {}'.format(delta_position)
                    error = self.tmp_d['error_b'][idx_datetime_] * 3600
                    error_fmt = float("{0:.6f}".format(error))  # error rounded
                    error_str = '   error {}'.format(error_fmt)
                    # annotation
                    ax.annotate('{}\n{}'.format(delta_str, error_str),
                                xy=(x, y), textcoords='data', fontsize=13)

                for idx_datetime_, datetime_ in enumerate(o_datetimes):
                    x = epoch_seconds[idx_datetime_]  # x position
                    y = o_delta_seconds[idx_datetime_]  # y position
                    error = float(self.tmp_d['error_b'][idx_datetime_] * 3600)
                    ax.errorbar(x, y, yerr=error,
                                ecolor='g', capthick=2, elinewidth=4)

                # Label creation
                y_label_ra = 'Declination (")\n'
                major_s = y_ticks['major_s']
                y_label_major_step = 'major step size {}"\n'.format(major_s)
                minor_s = y_ticks['minor_s']
                y_label_minor_step = 'minor step size {}"'.format(minor_s)
                ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                              y_label_minor_step))

                # Plots data
                ax.plot(epoch_seconds, i_delta_seconds, 'bs', markersize=6,
                        label='extracted position')
                ax.plot(epoch_seconds, o_delta_seconds, 'rs', markersize=6,
                        label='extracted position')
                ax.plot(epoch_seconds, i_delta_seconds, linewidth=1)
                # In order to avoid scientific notation plot should be redrawn
                ax = plt.gca()
                ax.get_xaxis().get_major_formatter().set_useOffset(False)
                plt.legend(loc=0, ncol=2, borderaxespad=0.)
                plt.draw()

                pdf.savefig()  # saves current figure

            plt.clf()  # clear current figure
            plt.close(fig)
