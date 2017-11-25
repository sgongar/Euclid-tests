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

from astropy import units as u
from astropy.coordinates import Angle
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from numpy import arange, array, median


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


def significant_l(number):
    """

    :param number:
    :return:
    """
    len_ = len(str(number))  # make sure there is enough precision
    a = ('%.' + str(len_) + 'E') % Decimal(number)
    significate_d = a.split(".")[0]
    times = a.split("E")[1]
    result = int(significate_d) * (10 ** int(times))

    return result


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
    # divisions = 4  # TODO

    x_difference = float(max(epoch_seconds)) - float(min(epoch_seconds))
    x_major_stp = (x_difference / 3)
    x_minor_stp = (x_difference / 6)

    # X-SCALE
    # Major steps
    x_major_stps = arange((epoch_seconds[0] - x_major_stp * 1),
                          (epoch_seconds[0] + x_major_stp * 4), x_major_stp)
    # Minor steps
    x_minor_stps = arange((epoch_seconds[0] - x_minor_stp * 2),
                          (epoch_seconds[0] - x_minor_stp * 8),
                          x_minor_stp)
    x_ticks = {'major_t': x_major_stps, 'minor_t': x_minor_stps,
               'major_s': x_major_stp, 'minor_s': x_minor_stp}

    return x_ticks


def create_y_ticks(alpha_seconds, hour, minute):
    """
    step = stp

    :param alpha_seconds:
    :param hour:
    :param minute:
    :return: y_ticks
    """
    divisions = 4  #

    # Gets the major step between ticks thought the difference
    # between the maximum and the minimum value of alpha
    difference = float(max(alpha_seconds)) - float(min(alpha_seconds))
    major_stp = (difference / divisions)
    #
    major_stp = round_number(major_stp)
    minor_stp = (major_stp / 4)

    # Gets maximum decimal position of major step
    decimals = int(str(major_stp)[::-1].find('.'))

    # Major step list starts two times before and end two times after
    # known values
    major_stps = arange(round(min(alpha_seconds), decimals) - major_stp * 2,
                        round(max(alpha_seconds), decimals) + major_stp * 2,
                        major_stp)
    # Minor step list starts eight times before and fend eight times
    # after know values
    minor_stps = arange(round(min(alpha_seconds), decimals) - minor_stp * 8,
                        round(max(alpha_seconds), decimals) + minor_stp * 8,
                        minor_stp)

    # Formats major steps for be readable by matplotlib axis
    alpha_major_stps = []
    for idx_stp, stp_ in enumerate(major_stps):
        alpha_major_stps.append('{}:{}:{}'.format(hour, minute, stp_))
    # Formats minor steps
    alpha_minor_stps = []
    for idx_stp, stp_ in enumerate(minor_stps):
        alpha_minor_stps.append('{}:{}:{}'.format(hour, minute, stp_))

    # Creates a list of datatime objects
    # Sometimes, seconds are greater than 59, this values should be
    # filter in a better way
    y_ticks = {'major_t': [datetime.strptime(t,
                                             "%H:%M:%S.%f") for t in alpha_major_stps],
               'minor_t': [datetime.strptime(t,
                                             "%H:%M:%S.%f") for t in alpha_minor_stps],
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
            ax.set_title('chi_squared: {}'.format(self.fitted_d['ra']))

            alpha_tmpstmp = []
            alpha_seconds = []
            tmp_hour = []
            tmp_minute = []
            for alpha_ in self.tmp_d['alpha']:
                a = Angle(alpha_, u.degree)
                hms = a.hms
                hour = int(hms[0])
                tmp_hour.append(hour)
                minute = int(hms[1])
                tmp_minute.append(minute)
                second = float("{0:.6f}".format(hms[2]))
                alpha_tmpstmp.append('{}:{}:{}'.format(hour, minute, second))
                alpha_seconds.append('{}'.format(second))

            datetimes = [datetime.strptime(t, "%H:%M:%S.%f") for t in alpha_tmpstmp]

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
            y_ticks = create_y_ticks(alpha_seconds, tmp_hour[0], tmp_minute[0])
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
            for idx_datetime_, datetime_ in enumerate(datetimes):
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

            for idx_datetime_, datetime_ in enumerate(datetimes):
                error_seconds = self.tmp_d['error_a'][idx_datetime_] * 3600
                errorbar_size = timedelta(0, error_seconds)
                ax.errorbar(epoch_seconds[idx_datetime_], datetime_,
                            yerr=errorbar_size,
                            fmt='o', ecolor='g', capthick=2, elinewidth=4)

            # Axis labels creation
            # x-axis
            x_label = 'EPOCH (seconds)'
            ax.set_xlabel(x_label)
            # y-axis
            y_label_ra = 'Right ascension (")\n'
            y_label_major_step = 'major step size {}"\n'.format(y_ticks['major_s'])
            y_label_minor_step = 'minor step size {}"'.format(y_ticks['minor_s'])
            ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                          y_label_minor_step))
            # Plots data
            ax.plot(epoch_seconds, datetimes, 'bs', markersize=6,
                    label='extracted position')
            # In order to avoid scientific notation plot should be redrawn
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)

            # Plots a legend
            plt.legend(loc=0, ncol=2, borderaxespad=0.)
            plt.draw()

            # Saves the current figure in pdf file
            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure

            #
            #
            # DELTA PARAMETERS
            #
            #
            fig = plt.figure(figsize=plot_size, dpi=plot_dpi)
            ax = fig.add_subplot(1, 1, 1)
            ax.set_title('chi_squared: {}'.format(self.fitted_d['dec']))

            delta_tmpstmp = []
            delta_seconds = []
            for delta_ in self.tmp_d['delta']:
                d = Angle(delta_, u.degree)
                dms = d.dms
                degree = int(dms[0])
                minute = int(dms[1])
                second = float("{0:.6f}".format(dms[2]))
                delta_tmpstmp.append('{}:{}.{}'.format(degree, minute, second))
                delta_seconds.append('{}'.format(float("{0:.6f}".format(dms[2]))))

            ax.set_xlabel('EPOCH')

            # Y-SCALE
            # format elements in alpha_seconds list to floats
            delta_seconds = [float(i) for i in delta_seconds]
            delta_seconds = [float("{0:.6f}".format(i)) for i in delta_seconds]

            y_ticks = create_y_ticks(delta_seconds, int(dms[0]), int(dms[1]))

            print(len(y_ticks))

            # x-ticks assignation
            ax.set_xticks(x_ticks['major_t'], minor=False)
            ax.set_xticks(x_ticks['minor_t'], minor=True)
            # x-ticks labels
            # TODO
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

            for idx_datetime_, datetime_ in enumerate(datetimes):
                x = epoch_seconds[idx_datetime_]  # x position
                y = delta_seconds[idx_datetime_]  # y position
                # variables for position and error representation
                delta_position = delta_tmpstmp[idx_datetime_]
                delta_str = '   delta {}'.format(delta_position)
                error = self.tmp_d['error_b'][idx_datetime_] * 3600  # error on y
                error_fmt = float("{0:.6f}".format(error))  # error rounded
                error_str = '   error {}'.format(error_fmt)
                # annotation
                ax.annotate('{}\n{}'.format(delta_str, error_str),
                            xy=(x, y), textcoords='data', fontsize=13)

            for idx_datetime_, datetime_ in enumerate(datetimes):
                x = epoch_seconds[idx_datetime_]  # x position
                y = delta_seconds[idx_datetime_]  # y position
                error = float(self.tmp_d['error_b'][idx_datetime_] * 3600)
                ax.errorbar(x, y, yerr=error, fmt='o',
                            ecolor='g', capthick=2, elinewidth=4)

            # Label creation
            y_label_ra = 'Declination (")\n'
            y_label_major_step = 'major step size {}"\n'.format(y_ticks['major_s'])
            y_label_minor_step = 'minor step size {}"'.format(y_ticks['minor_s'])
            ax.set_ylabel('{}{}{}'.format(y_label_ra, y_label_major_step,
                                          y_label_minor_step))

            ax.plot(epoch_seconds, delta_seconds, 'bs', markersize=6,
                    label='extracted position')

            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            plt.legend(loc=0, ncol=2, borderaxespad=0.)
            plt.draw()

            pdf.savefig()  # saves current figure
            plt.clf()  # clear current figure
