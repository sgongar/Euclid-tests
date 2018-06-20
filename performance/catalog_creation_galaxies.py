#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a catalog populated of galaxies from sextracted catalogs
from single CCDs images.

Versions:
- 0.1: Initial release. Split from stars_catalog_creation.py

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    *

*GNU Terry Pratchett*

"""
from math import hypot
from multiprocessing import Process
from sys import stdout

from astropy.io import fits
from astropy.table import Table
from numpy import median
from pandas import concat, DataFrame, read_csv

from misc import get_cat, get_cats, extract_settings_elvis

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def check_distance(o_df, alpha, delta):
    """

    :param o_df:
    :param alpha:
    :param delta:
    :return:
    """
    distance_l = []
    for ix, row in o_df.iterrows():
        distance = hypot(row.ALPHA_J2000 - alpha, row.DELTA_J2000 - delta)
        distance_l.append(distance)

    index = distance_l.index(min(distance_l))

    return index


def check_source(o_cat, i_alpha, i_delta):
    """

    :param o_cat:
    :param i_alpha:
    :param i_delta:
    :return:
    """
    prfs_d = extract_settings_elvis()

    o_cat = o_cat[o_cat['ALPHA_J2000'] + prfs_d['tolerance'] > i_alpha]
    o_cat = o_cat[i_alpha > o_cat['ALPHA_J2000'] - prfs_d['tolerance']]
    o_cat = o_cat[o_cat['DELTA_J2000'] + prfs_d['tolerance'] > i_delta]
    o_cat = o_cat[i_delta > o_cat['DELTA_J2000'] - prfs_d['tolerance']]

    return o_cat


def extract_cats_d():
    """

    :return:
    """
    cats_d = {}
    for dither in range(1, 5, 1):
        cats_d[dither] = {}
        cats = get_cats(dither)
        for cat_name in cats:
            hdu_list = fits.open('{}/{}'.format(prfs_dict['fits_dir'],
                                                cat_name))
            cat_data = Table(hdu_list[2].data)
            cat_df = cat_data.to_pandas()  # Converts to Pandas format
            cat_number = get_cat(cat_name)  # Gets cat's number from cat's name
            cats_d[dither][cat_name] = cat_df

            cat_list = [cat_number] * cat_df['NUMBER'].size
            cats_d[dither][cat_name]['CATALOG_NUMBER'] = cat_list

    return cats_d


def create_full_cats(cats_d):
    """

    :param cats_d:
    :return:
    """
    full_d = {}

    for dither in range(1, 5, 1):
        dither_l = []
        for key_ in cats_d[dither].keys():
            dither_l.append(cats_d[dither][key_])
        full_d[dither] = concat(dither_l, ignore_index=True)
        full_idx = range(0, full_d[dither]['NUMBER'].size, 1)
        full_d[dither]['IDX'] = full_idx

    return full_d


def extract_inputs_d():
    """

    :return:
    """
    inputs_d = {}

    cat_galaxies_loc = prfs_dict['references']
    cat_galaxies = fits.open('{}/cat_galaxies.fits'.format(cat_galaxies_loc))
    galaxies_data = Table(cat_galaxies[1].data)
    galaxies_df = galaxies_data.to_pandas()
    galaxies_idx = range(0, 143766, 1)  # hardcoded - todo!
    galaxies_df['IDX'] = galaxies_idx
    inputs_d['galaxies'] = galaxies_df

    return inputs_d


def create_empty_catalog_dict():
    """

    :return: cat_d
    """
    cat_d = {'NUMBER': [], 'IDX': [], 'CATALOG_NUMBER': [], 'X_WORLD': [],
             'Y_WORLD': [], 'MAG_AUTO': [], 'MAGERR_AUTO': [],
             'ERRA_WORLD': [], 'ERRB_WORLD': [], 'ERRTHETA_WORLD': []}

    return cat_d


def create_catalog():
    """

    :return:
    """
    cats_d = extract_cats_d()  # extracts dataframes from catalogues
    full_d = create_full_cats(cats_d)  # creates dataframe from CCDs catalogues
    inputs_d = extract_inputs_d()
    cats = {}
    save = True

    unique_sources = inputs_d['galaxies']['IDX']
    total_galaxies = inputs_d['galaxies']['IDX'].size

    sub_list_size = total_galaxies / 18

    sub_list_l = []
    for idx_sub_list in range(0, 18, 1):
        if idx_sub_list != (18 - 1):
            idx_down = sub_list_size * idx_sub_list
            idx_up = sub_list_size * (idx_sub_list + 1)
            sub_list_l.append(unique_sources[idx_down:idx_up])
        else:
            idx_down = sub_list_size * idx_sub_list
            sub_list_l.append(unique_sources[idx_down:])

    areas_j = []
    for idx_l in range(0, 18, 1):
        areas_p = Process(target=create_galaxies_catalog_thread,
                          args=(idx_l, sub_list_l[idx_l], inputs_d, full_d))
        areas_j.append(areas_p)
        areas_p.start()

    active_areas = list([job.is_alive() for job in areas_j])
    while True in active_areas:
        active_areas = list([job.is_alive() for job in areas_j])
        pass

    # Merges areas
    # Merges catalogs
    galaxies_list = []
    for idx_csv in range(0, 18, 1):
        galaxies_ = read_csv('tmp_galaxies/galaxies_{}.csv'.format(idx_csv),
                             index_col=0)
        galaxies_list.append(galaxies_)

    galaxies_df = concat(galaxies_list)
    cats['galaxies'] = galaxies_df

    if save:
        galaxies_df.to_csv('tmp_galaxies/galaxies.csv')

    return cats


def create_galaxies_catalog_thread(idx_l, sub_list, inputs_d, full_d):
    """

    :param idx_l:
    :param sub_list:
    :param inputs_d:
    :param full_d:
    :return:
    """
    save = True

    cat_d = create_empty_catalog_dict()
    total_thread = len(sub_list)
    stdout.write('total galaxies {} of thread {}\n'.format(total_thread,
                                                           idx_l))
    for idx, galaxy in enumerate(sub_list):
        galaxies_df = inputs_d['galaxies']
        source_df = galaxies_df[galaxies_df['IDX'].isin([galaxy])]
        alpha = source_df['ra'].iloc[0]
        delta = source_df['dec'].iloc[0]

        dither_n = 0
        source_d = create_empty_catalog_dict()
        for dither in range(1, 5, 1):
            o_df = check_source(full_d[dither], alpha, delta)
            if o_df.empty is not True:
                dither_n += 1
                # Returns the index of the closest found source
                index = check_distance(o_df, alpha, delta)
                o_df = o_df.iloc[[index]]
                for key_ in source_d.keys():
                    source_d[key_].append(o_df[key_].iloc[0])

        if len(source_d['NUMBER']) != 0:
            for key_ in source_d.keys():
                median_value = median(source_d[key_])
                cat_d[key_].append(median_value)

    cat_df = DataFrame(cat_d)
    if save:
        cat_df.to_csv('tmp_galaxies/galaxies_{}.csv'.format(idx_l))


def write_galaxies_catalog(catalogs):
    """

    :param catalogs:
    :return:
    """
    galaxies_df_data = catalogs['galaxies']
    # Galaxies catalogue creation
    test_cat_name = '{}/coadd.cat'.format(prfs_dict['references'])
    test_coadd_cat = fits.open(test_cat_name)

    # Source number
    # c1 = fits.Column(name='NUMBER', format='1J', disp='I10',
    #                  array=stars_df_data['NUMBER'])
    # Kron-like elliptical aperture magnitude
    c1 = fits.Column(name='MAG_AUTO', format='1E', unit='mag',
                     disp='F8.4', array=galaxies_df_data['MAG_AUTO'])
    # RMS error for AUTO magnitude
    c2 = fits.Column(name='MAGERR_AUTO', format='1E', unit='mag',
                     disp='F8.4', array=galaxies_df_data['MAGERR_AUTO'])
    # Barycenter position along world x axis
    c3 = fits.Column(name='X_WORLD', format='1D', unit='deg', disp='E18.10',
                     array=galaxies_df_data['X_WORLD'])
    # Barycenter position along world y axis
    c4 = fits.Column(name='Y_WORLD', format='1D', unit='deg', disp='E18.10',
                     array=galaxies_df_data['Y_WORLD'])
    # World RMS position error along major axis
    c5 = fits.Column(name='ERRA_WORLD', format='1E', unit='deg',
                     disp='G12.7', array=galaxies_df_data['ERRA_WORLD'])
    # World RMS position error along minor axis
    c6 = fits.Column(name='ERRB_WORLD', format='1E', unit='deg',
                     disp='G12.7', array=galaxies_df_data['ERRB_WORLD'])
    # Error ellipse pos.angle(CCW / world - x)
    c7 = fits.Column(name='ERRTHETA_WORLD', format='1E', unit='deg',
                     disp='F6.2', array=galaxies_df_data['ERRTHETA_WORLD'])

    col_defs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7])

    tb_hdu = fits.BinTableHDU.from_columns(col_defs)

    test_coadd_cat[2] = tb_hdu
    test_coadd_cat[2].header['EXTNAME'] = 'LDAC_OBJECTS'

    newcat_name = '{}/galaxies_catalogue.cat'.format(prfs_dict['references'])
    test_coadd_cat.writeto(newcat_name, overwrite=True)


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    cats = create_catalog()
    write_galaxies_catalog(cats)