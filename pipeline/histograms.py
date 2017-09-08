#!/usr/bin/python
# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.table import Table
import numpy as np
from misc import extract_settings
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from misc import get_ticks, pm_compute, pm_filter, get_limits
from pandas import concat


"""
def database_extract(prefs_dict):
    """

    """
    full_catalogue = fits.open(prefs_dict['results_dir'] + '/full_1.cat')
    full_table = full_catalogue[2].data
    data_table = Table(full_table).to_pandas()

    merged_catalogue = fits.open(prefs_dict['results_dir'] + '/merged_1.cat')
    merged_table = merged_catalogue[2].data
    data_merged = Table(merged_table).to_pandas()

    # removing "zeros"
    db = data_table.loc[~data_table['CATALOG_NUMBER'].isin([0])]

    # computing pm
    db = pm_compute(merged_table, db)

    db = concat(g for _, g in db.groupby("SOURCE_NUMBER")
                if len(g) >= int(prefs_dict['detections']))

    db = pm_filter(db, prefs_dict['pm_low'], prefs_dict['pm_up'],
                   prefs_dict['pm_sn'])

    # pm_list = list(set(list(db['PM'])))
    mag_list = []
    pm_list = []

    g = db.groupby('SOURCE_NUMBER')
    d = dict(iter(g))

    flags_g = data_merged.groupby('FLAGS_SCAMP')
    flags_d = dict(iter(flags_g))

    flags_list = []

    for i in range(len(d.keys()) - 1):
        mag = float(data_merged['MAG'].iloc[int(d.keys()[i]) - 1])
        flag_64 = float(data_merged['FLAGS_SCAMP'].iloc[int(d.keys()[i]) - 1])
        flag_0 = float(data_merged['FLAGS_SCAMP'].iloc[int(d.keys()[i]) - 1])

        check_mag = mag != 99.0
        check_flag_64 = flag_64 == 64.0
        check_flag_0 = flag_0 == 0.0

        if check_mag and check_flag_64 or check_mag and check_flag_0:
            if db.loc[db['SOURCE_NUMBER'] == int(d.keys()[i]) - 1, 'PM'].size > 1:
                mag_list.append(data_merged['MAG'].iloc[int(d.keys()[i]) - 1])
                pm_list.append(db.loc[db['SOURCE_NUMBER'] == int(d.keys()[i]) - 1, 'PM'].iloc[0])

        flags_list.append(float(data_merged['FLAGS_SCAMP'].iloc[int(d.keys()[i])]))

    return pm_list, mag_list, flags_list, flags_d


def get_input_catalogue(prefs_dict):
    """

    """
    catalogue_file = np.genfromtxt(prefs_dict['input_cats'] + '/Cat_d1.dat')
    list_magnitudes = catalogue_file[:, 2]
    list_pm = catalogue_file[:, 3]

    pm_values = []
    mag_values = []

    for i in range(len(list_pm)):
        if list_pm[i] != 0.:
            pm_values.append(list_pm[i] - 1000)
            mag_values.append(list_magnitudes[i])

    return pm_values, mag_values


def draw_histograms(prefs_dict, pm_in_values, mag_in_values):
    """
    todo share x axis!
    """

    # values
    pm_out_values, mag_out_values, flags_values, flags_dict = database_extract(prefs_dict)

    # shared limits and axis for proper motion values
    limits = get_limits(pm_in_values, pm_out_values)
    x0_limits = [limits[0] - 10, limits[1] + 10]
    x0_ticks = get_ticks(limits[0], limits[1], 20)

    # shared limits and axis for magnitudes values
    limits = get_limits(mag_in_values, mag_out_values)
    x1_limits = [limits[0] - 1, limits[1] + 1]
    x1_ticks = get_ticks(limits[0], limits[1], 1)

    # subplots
    fig, ((ax0, ax1), (ax2, ax3), (ax4, ax5)) = plt.subplots(ncols=2, nrows=3,
                                                             figsize=(8.27,
                                                                      11.69),
                                                             dpi=100)

    # input analysis values
    ax0.hist(pm_in_values, bins=20,
             facecolor='gray', label='test')
    ax0.set_title("Proper motion distribution")
    ax0.set_xlim(x0_limits)
    x0_ticks.insert(0, 20)  # TODO workaround about automatic range
    ax0.set_xticks(x0_ticks)
    ax0.grid(True)

    ax1.hist(mag_in_values, bins=20,
             facecolor='gray')
    ax1.set_title("Magnitude distribution")
    ax1.set_xlim(x1_limits)
    ax1.set_xticks(x1_ticks)

    ax1.grid(True)

    # output analysis values
    print "proper motion out values", len(pm_out_values)
    ax2.hist(pm_out_values, bins=20,
             facecolor='gray', label='test')
    ax2.set_title("Proper motion distribution")
    ax2.set_xlim(x0_limits)
    ax2.set_xticks(x0_ticks)
    ax2.grid(True)

    print "magnitude out values", len(mag_out_values)
    ax3.hist(mag_out_values, bins=20,
             facecolor='gray')
    ax3.set_title("Magnitude distribution")
    ax3.set_xlim(x1_limits)
    ax3.set_xticks(x1_ticks)

    ax3.grid(True)

    print "flags distribution", len(flags_values)
    ax4.hist(flags_values, bins=20,
             facecolor='gray')
    ax4.set_title("Flags distribution")
    # ax4.set_xlim(x1_limits)
    # ax4.set_xticks(x1_ticks)

    ax4.grid(True)

    print "flags distribution", len(flags_dict.keys())
    labels = flags_dict.keys()
    sizes = [len(flags_dict[flags_dict.keys()[x]]) for x in range(len(flags_dict.keys()))]
    ax5.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax5.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    ax5.grid(True)

    fig.tight_layout()

    with PdfPages('test.pdf') as pdf:
        pdf.savefig()


def extract_data(prefs_dict):
    pm_values, mag_values = get_input_catalogue(prefs_dict)
    draw_histograms(prefs_dict, pm_values, mag_values)


if __name__ == '__main__':
    prefs_dict = extract_settings()
    extract_data(prefs_dict)
"""
