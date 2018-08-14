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

from astropy.units import degree
from astropy.coordinates import SkyCoord
from numpy import pi
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
    cat_d = {'IDX': [], 'SOURCE': [], 'DITHER': [], 'RA': [],
             'DEC': [], 'VEL': [], 'ABMAG': [], 'THETA': []}

    return cat_d


def extract_ssos_df():
    """

    :return:
    """
    # columns = ['ALPHA_J2000', 'DELTA_J2000', 'PM', 'PA',
    #            'MAG', 'MAG', 'MAG', 'MAG']
    # TODO actual catalogues is only valid for dither -> 1
    # TODO create catalogues for dithers 2, 3 and 4
    # cat_ssos = read_csv('{}/ssos_cat.txt'.format(prfs_dict['references']),
    #                     delim_whitespace=True, header=None, names=columns)

    cat_ssos = read_csv('{}/ssos_cat.txt'.format(prfs_dict['references']),
                        delim_whitespace=True)
    print(cat_ssos.columns)
    ssos_source = range(0, cat_ssos['RA'].size, 1)
    cat_ssos['SOURCE'] = ssos_source
    ssos_df = cat_ssos
    # ssos_df = cat_ssos[['SOURCE', 'ALPHA_J2000', 'DELTA_J2000',
    #                    'PM', 'PA', 'MAG']]

    # Hasta aqui bien
    alpha_list = Series(ssos_df['RA'].tolist(), name='ALPHA_J2000')
    delta_list = Series(ssos_df['DEC'].tolist(), name='DELTA_J2000')

    positions_table = concat([alpha_list, delta_list], axis=1)
    positions_table.to_csv('test_1.reg', index=False, header=False, sep=" ")

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
        pm_source = source_df['VEL'].iloc[0]
        ssos_d['VEL'].append(pm_source)
        pa_source = source_df['THETA'].iloc[0]
        ssos_d['THETA'].append(pa_source)
        mag_source = source_df['ABMAG'].iloc[0]
        ssos_d['ABMAG'].append(mag_source)

        dither_time = (565.0 / 2) / 3600.0
        alpha_source = source_df['RA'].iloc[0]
        alpha_increment_per_hour = cos(
            float(pa_source * pi / 180.0)) * float(pm_source)
        alpha_increment_per_dither = alpha_increment_per_hour * dither_time
        alpha_source = alpha_source + (alpha_increment_per_dither / 3600.0)
        ssos_d['RA'].append(alpha_source)
        delta_source = source_df['DEC'].iloc[0]
        delta_increment_per_hour = sin(
            float(pa_source * pi / 180.0)) * float(pm_source)
        delta_increment_per_dither = delta_increment_per_hour * dither_time
        delta_source = delta_source + (delta_increment_per_dither / 3600.0)
        ssos_d['DEC'].append(delta_source)

        for dither_source in range(2, 5, 1):
            # dither_time is equal to fraction of hour
            dither_time = 1003.0/3600.0

            idx += 1
            alpha_increment_per_hour = cos(float(pa_source*pi/180.0)) * float(pm_source)
            alpha_increment_per_dither = alpha_increment_per_hour * dither_time
            alpha_source = alpha_source + (alpha_increment_per_dither / 3600.0)
            delta_increment_per_hour = sin(float(pa_source*pi/180.0)) * float(pm_source)
            delta_increment_per_dither = delta_increment_per_hour * dither_time
            delta_source = delta_source + (delta_increment_per_dither / 3600.0)

            ssos_d['IDX'].append(idx)
            ssos_d['SOURCE'].append(idx_source)
            ssos_d['DITHER'].append(dither_source)
            ssos_d['RA'].append(alpha_source)
            ssos_d['DEC'].append(delta_source)
            ssos_d['VEL'].append(pm_source)
            ssos_d['THETA'].append(pa_source)
            ssos_d['ABMAG'].append(mag_source)

        idx += 1

    sso_cat = DataFrame(ssos_d)

    return sso_cat


# This is wrong! New SSOs catalog is coming!
def filter_by_position(sso_df):
    """

    :param sso_df:
    :return:
    """
    save = True
    borders_d = get_borders()

    right_sources = []

    unique_sources = list(set(sso_df['SOURCE']))

    for idx_source_, source_ in enumerate(unique_sources):
        source_df = sso_df[sso_df['SOURCE'].isin([source_])]
        for idx_dither_, row in enumerate(source_df.itertuples(), 1):
            alpha = row.RA
            delta = row.DEC
            source_coords = SkyCoord(ra=alpha * degree, dec=delta * degree,
                                     equinox='J2021.5')
            source_coords_ecliptic = source_coords.barycentrictrueecliptic
            lon_e = float(source_coords_ecliptic.lon.degree)
            lat_e = float(source_coords_ecliptic.lat.degree)
            for ccd_ in borders_d[idx_dither_].keys():
                borders = borders_d[idx_dither_][ccd_]
                lon_comp = borders['below_lon'] < lon_e < borders['above_lon']
                lat_comp = borders['below_lat'] < lat_e < borders['above_lat']
                comp = lon_comp and lat_comp

                if comp:
                    right_sources.append(row.IDX)

    if save:
        sso_df.to_csv('catalogues_detected/cat_ssos.csv')
    for dither_ in range(1, 5, 1):
        catalog = sso_df[sso_df['DITHER'].isin([dither_])]
        catalog.to_csv('catalogues_detected/cat_ssos_{}.csv'.format(dither_))

    # Removes non visible sources
    sso_clean_df = sso_df[sso_df['IDX'].isin(right_sources)]
    if save:
        sso_clean_df.to_csv('catalogues_detected/cat_clean_ssos.csv')
    for dither_ in range(1, 5, 1):
        clean_catalog = sso_clean_df[sso_clean_df['DITHER'].isin([dither_])]
        clean_catalog.to_csv('catalogues_detected/cat_clean_ssos_{}.csv'.format(dither_))


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    sso_df = propagate_dithers()
    filter_by_position(sso_df)  # Rejects source if is out of the bounds
