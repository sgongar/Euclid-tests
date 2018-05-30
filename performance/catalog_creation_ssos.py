#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release. Split from ssos_catalog_creation.py

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    *

*GNU Terry Pratchett*

"""
from math import cos, sin

from pandas import concat, DataFrame, read_csv, Series

from images_management_elvis import get_borders
from misc import extract_settings_elvis

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def create_empty_catalog_dict():
    """

    :return: cat_d
    """
    cat_d = {'IDX': [], 'SOURCE': [], 'DITHER': [], 'ALPHA_J2000': [],
             'DELTA_J2000': [], 'PM': [], 'MAG': [], 'PA': []}

    return cat_d


def extract_ssos_df():
    """

    :return:
    """
    columns = ['ALPHA_J2000', 'DELTA_J2000', 'PM', 'PA',
               'MAG', 'MAG', 'MAG', 'MAG']
    # TODO actual catalogues is only valid for dither -> 1
    # TODO create catalogues for dithers 2, 3 and 4
    cat_ssos = read_csv('{}/ssos_cat.txt'.format(prfs_dict['references']),
                        delim_whitespace=True, header=None, names=columns)

    ssos_source = range(0, cat_ssos['ALPHA_J2000'].size, 1)
    cat_ssos['SOURCE'] = ssos_source
    ssos_df = cat_ssos[['SOURCE', 'ALPHA_J2000', 'DELTA_J2000',
                        'PM', 'PA', 'MAG']]

    return ssos_df


def propagate_dithers():
    """

    :return:
    """
    ssos_d = create_empty_catalog_dict()
    ssos_df = extract_ssos_df()
    unique_sources = ssos_df['SOURCE']
    idx = 0
    # Move over sources and dithers
    for idx_source, source_ in enumerate(unique_sources):
        source_df = ssos_df[ssos_df['SOURCE'].isin([source_])]
        # Dither 1
        ssos_d['IDX'].append(idx)
        idx_source = source_df['SOURCE'].iloc[0]
        ssos_d['SOURCE'].append(idx_source)
        dither_source = 1
        ssos_d['DITHER'].append(dither_source)
        alpha_source = source_df['ALPHA_J2000'].iloc[0]
        ssos_d['ALPHA_J2000'].append(alpha_source)
        delta_source = source_df['DELTA_J2000'].iloc[0]
        ssos_d['DELTA_J2000'].append(delta_source)
        pm_source = source_df['PM'].iloc[0]
        ssos_d['PM'].append(pm_source)
        pa_source = source_df['PA'].iloc[0]
        ssos_d['PA'].append(pa_source)
        mag_source = source_df['MAG'].iloc[0]
        ssos_d['MAG'].append(mag_source)

        for dither_source in range(2, 5, 1):
            # dither_time is equal to fraction of hour
            dither_time = 0.27861
            idx += 1
            alpha_increment_per_hour = cos(float(pa_source)) * float(pm_source)
            alpha_increment_per_dither = alpha_increment_per_hour / dither_time
            alpha_source = alpha_source + (alpha_increment_per_dither / 3600)
            delta_increment_per_hour = sin(float(pa_source)) * float(pm_source)
            delta_increment_per_dither = delta_increment_per_hour / dither_time
            delta_source = delta_source + (delta_increment_per_dither / 3600)

            ssos_d['IDX'].append(idx)
            ssos_d['SOURCE'].append(idx_source)
            ssos_d['DITHER'].append(dither_source)
            ssos_d['ALPHA_J2000'].append(alpha_source)
            ssos_d['DELTA_J2000'].append(delta_source)
            ssos_d['PM'].append(pm_source)
            ssos_d['PA'].append(pa_source)
            ssos_d['MAG'].append(mag_source)

        idx += 1

    for key_ in ssos_d.keys():
        print(key_, len(ssos_d[key_]))

    sso_cat = DataFrame(ssos_d)

    return sso_cat


def filter_by_position(sso_df):
    """

    :param sso_df:
    :return:
    """
    save = True
    borders_d = get_borders()
    right_sources = []

    unique_sources = list(set(sso_df['IDX']))

    for idx_source_, source_ in enumerate(unique_sources):
        source_df = sso_df[sso_df['IDX'].isin([source_])]
        for idx_dither_, row in enumerate(source_df.itertuples(), 1):
            alpha = row.ALPHA_J2000
            delta = row.DELTA_J2000
            for ccd_ in borders_d[idx_dither_].keys():
                borders = borders_d[idx_dither_][ccd_]
                alpha_comp = borders['below_ra'] < alpha < borders['above_ra']
                delta_comp = borders['below_dec'] < delta < borders['above_dec']
                comp = alpha_comp and delta_comp

                if comp:
                    right_sources.append(row.IDX)

    if save:
        sso_df.to_csv('cat_ssos.csv')

    sso_clean_df = sso_df[sso_df['IDX'].isin(right_sources)]

    if save:
        sso_clean_df.to_csv('cat_clean_ssos.csv')


def create_regions():
    """

    :return:
    """
    for dither_ in range(1, 5, 1):
        catalog = read_csv('cat_ssos.csv', index_col=0)
        catalog = catalog[catalog['DITHER'].isin([dither_])]
        alpha_list = Series(catalog['ALPHA_J2000'].tolist(),
                            name='ALPHA_J2000')
        delta_list = Series(catalog['DELTA_J2000'].tolist(),
                            name='DELTA_J2000')

        positions_table = concat([alpha_list, delta_list], axis=1)
        positions_table.to_csv('cat_ssos_{}.reg'.format(dither_),
                               index=False, header=False, sep=" ")

    for dither_ in range(1, 5, 1):
        clean_catalog = read_csv('cat_clean_ssos.csv', index_col=0)
        clean_catalog = clean_catalog[clean_catalog['DITHER'].isin([dither_])]
        alpha_list = Series(clean_catalog['ALPHA_J2000'].tolist(),
                            name='ALPHA_J2000')
        delta_list = Series(clean_catalog['DELTA_J2000'].tolist(),
                            name='DELTA_J2000')

        positions_table = concat([alpha_list, delta_list], axis=1)
        positions_table.to_csv('cat_clean_ssos_{}.reg'.format(dither_),
                               index=False, header=False, sep=" ")


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    sso_df = propagate_dithers()
    filter_by_position(sso_df)  # Rejects source if is out of the bounds

    create_regions()
