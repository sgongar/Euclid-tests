#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

La idea es que compare los objetos extraidos en total con los que tengo.

Todo:
    * Improve log messages

"""

from astropy.io import fits
from astropy.table import Table
from pandas import read_csv

from misc import extract_settings_sc3

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def get_cat(cat_n):
    """  TODO get scamp organization in order to reduce the lines number

    :param cat_n:
    :return: cat_file
    """
    cat_part_1 = 'EUC_VIS_SWL-DET-00'

    cats_final = []
    cats_partial = ['', '1-000000-0000000__20170630T011437.3Z_00.00_0',
                    '1-000000-0000000__20170630T011437.3Z_00.00_10',
                    '1-000000-0000000__20170630T011437.3Z_00.00_11',
                    '1-000000-0000000__20170630T011437.3Z_00.00_12',
                    '1-000000-0000000__20170630T011437.3Z_00.00_13',
                    '1-000000-0000000__20170630T011437.3Z_00.00_14',
                    '1-000000-0000000__20170630T011437.3Z_00.00_15',
                    '1-000000-0000000__20170630T011437.3Z_00.00_16',
                    '1-000000-0000000__20170630T011437.3Z_00.00_17',
                    '1-000000-0000000__20170630T011437.3Z_00.00_18',
                    '1-000000-0000000__20170630T011437.3Z_00.00_19',
                    '1-000000-0000000__20170630T011437.3Z_00.00_1',
                    '1-000000-0000000__20170630T011437.3Z_00.00_20',
                    '1-000000-0000000__20170630T011437.3Z_00.00_21',
                    '1-000000-0000000__20170630T011437.3Z_00.00_22',
                    '1-000000-0000000__20170630T011437.3Z_00.00_23',
                    '1-000000-0000000__20170630T011437.3Z_00.00_24',
                    '1-000000-0000000__20170630T011437.3Z_00.00_25',
                    '1-000000-0000000__20170630T011437.3Z_00.00_26',
                    '1-000000-0000000__20170630T011437.3Z_00.00_27',
                    '1-000000-0000000__20170630T011437.3Z_00.00_28',
                    '1-000000-0000000__20170630T011437.3Z_00.00_29',
                    '1-000000-0000000__20170630T011437.3Z_00.00_2',
                    '1-000000-0000000__20170630T011437.3Z_00.00_30',
                    '1-000000-0000000__20170630T011437.3Z_00.00_31',
                    '1-000000-0000000__20170630T011437.3Z_00.00_32',
                    '1-000000-0000000__20170630T011437.3Z_00.00_33',
                    '1-000000-0000000__20170630T011437.3Z_00.00_34',
                    '1-000000-0000000__20170630T011437.3Z_00.00_35',
                    '1-000000-0000000__20170630T011437.3Z_00.00_3',
                    '1-000000-0000000__20170630T011437.3Z_00.00_4',
                    '1-000000-0000000__20170630T011437.3Z_00.00_5',
                    '1-000000-0000000__20170630T011437.3Z_00.00_6',
                    '1-000000-0000000__20170630T011437.3Z_00.00_7',
                    '1-000000-0000000__20170630T011437.3Z_00.00_8',
                    '1-000000-0000000__20170630T011437.3Z_00.00_9',
                    '2-000000-0000000__20170630T011642.0Z_00.00_0',
                    '2-000000-0000000__20170630T011642.0Z_00.00_10',
                    '2-000000-0000000__20170630T011642.0Z_00.00_11',
                    '2-000000-0000000__20170630T011642.0Z_00.00_12',
                    '2-000000-0000000__20170630T011642.0Z_00.00_13',
                    '2-000000-0000000__20170630T011642.0Z_00.00_14',
                    '2-000000-0000000__20170630T011642.0Z_00.00_15',
                    '2-000000-0000000__20170630T011642.0Z_00.00_16',
                    '2-000000-0000000__20170630T011642.0Z_00.00_17',
                    '2-000000-0000000__20170630T011642.0Z_00.00_18',
                    '2-000000-0000000__20170630T011642.0Z_00.00_19',
                    '2-000000-0000000__20170630T011642.0Z_00.00_1',
                    '2-000000-0000000__20170630T011642.0Z_00.00_20',
                    '2-000000-0000000__20170630T011642.0Z_00.00_21',
                    '2-000000-0000000__20170630T011642.0Z_00.00_22',
                    '2-000000-0000000__20170630T011642.0Z_00.00_23',
                    '2-000000-0000000__20170630T011642.0Z_00.00_24',
                    '2-000000-0000000__20170630T011642.0Z_00.00_25',
                    '2-000000-0000000__20170630T011642.0Z_00.00_26',
                    '2-000000-0000000__20170630T011642.0Z_00.00_27',
                    '2-000000-0000000__20170630T011642.0Z_00.00_28',
                    '2-000000-0000000__20170630T011642.0Z_00.00_29',
                    '2-000000-0000000__20170630T011642.0Z_00.00_2',
                    '2-000000-0000000__20170630T011642.0Z_00.00_30',
                    '2-000000-0000000__20170630T011642.0Z_00.00_31',
                    '2-000000-0000000__20170630T011642.0Z_00.00_32',
                    '2-000000-0000000__20170630T011642.0Z_00.00_33',
                    '2-000000-0000000__20170630T011642.0Z_00.00_34',
                    '2-000000-0000000__20170630T011642.0Z_00.00_35',
                    '2-000000-0000000__20170630T011642.0Z_00.00_3',
                    '2-000000-0000000__20170630T011642.0Z_00.00_4',
                    '2-000000-0000000__20170630T011642.0Z_00.00_5',
                    '2-000000-0000000__20170630T011642.0Z_00.00_6',
                    '2-000000-0000000__20170630T011642.0Z_00.00_7',
                    '2-000000-0000000__20170630T011642.0Z_00.00_8',
                    '2-000000-0000000__20170630T011642.0Z_00.00_9',
                    '3-000000-0000000__20170630T011848.6Z_00.00_0',
                    '3-000000-0000000__20170630T011848.6Z_00.00_10',
                    '3-000000-0000000__20170630T011848.6Z_00.00_11',
                    '3-000000-0000000__20170630T011848.6Z_00.00_12',
                    '3-000000-0000000__20170630T011848.6Z_00.00_13',
                    '3-000000-0000000__20170630T011848.6Z_00.00_14',
                    '3-000000-0000000__20170630T011848.6Z_00.00_15',
                    '3-000000-0000000__20170630T011848.6Z_00.00_16',
                    '3-000000-0000000__20170630T011848.6Z_00.00_17',
                    '3-000000-0000000__20170630T011848.6Z_00.00_18',
                    '3-000000-0000000__20170630T011848.6Z_00.00_19',
                    '3-000000-0000000__20170630T011848.6Z_00.00_1',
                    '3-000000-0000000__20170630T011848.6Z_00.00_20',
                    '3-000000-0000000__20170630T011848.6Z_00.00_21',
                    '3-000000-0000000__20170630T011848.6Z_00.00_22',
                    '3-000000-0000000__20170630T011848.6Z_00.00_23',
                    '3-000000-0000000__20170630T011848.6Z_00.00_24',
                    '3-000000-0000000__20170630T011848.6Z_00.00_25',
                    '3-000000-0000000__20170630T011848.6Z_00.00_26',
                    '3-000000-0000000__20170630T011848.6Z_00.00_27',
                    '3-000000-0000000__20170630T011848.6Z_00.00_28',
                    '3-000000-0000000__20170630T011848.6Z_00.00_29',
                    '3-000000-0000000__20170630T011848.6Z_00.00_2',
                    '3-000000-0000000__20170630T011848.6Z_00.00_30',
                    '3-000000-0000000__20170630T011848.6Z_00.00_31',
                    '3-000000-0000000__20170630T011848.6Z_00.00_32',
                    '3-000000-0000000__20170630T011848.6Z_00.00_33',
                    '3-000000-0000000__20170630T011848.6Z_00.00_34',
                    '3-000000-0000000__20170630T011848.6Z_00.00_35',
                    '3-000000-0000000__20170630T011848.6Z_00.00_3',
                    '3-000000-0000000__20170630T011848.6Z_00.00_4',
                    '3-000000-0000000__20170630T011848.6Z_00.00_5',
                    '3-000000-0000000__20170630T011848.6Z_00.00_6',
                    '3-000000-0000000__20170630T011848.6Z_00.00_7',
                    '3-000000-0000000__20170630T011848.6Z_00.00_8',
                    '3-000000-0000000__20170630T011848.6Z_00.00_9',
                    '4-000000-0000000__20170630T012050.1Z_00.00_0',
                    '4-000000-0000000__20170630T012050.1Z_00.00_10',
                    '4-000000-0000000__20170630T012050.1Z_00.00_11',
                    '4-000000-0000000__20170630T012050.1Z_00.00_12',
                    '4-000000-0000000__20170630T012050.1Z_00.00_13',
                    '4-000000-0000000__20170630T012050.1Z_00.00_14',
                    '4-000000-0000000__20170630T012050.1Z_00.00_15',
                    '4-000000-0000000__20170630T012050.1Z_00.00_16',
                    '4-000000-0000000__20170630T012050.1Z_00.00_17',
                    '4-000000-0000000__20170630T012050.1Z_00.00_18',
                    '4-000000-0000000__20170630T012050.1Z_00.00_19',
                    '4-000000-0000000__20170630T012050.1Z_00.00_1',
                    '4-000000-0000000__20170630T012050.1Z_00.00_20',
                    '4-000000-0000000__20170630T012050.1Z_00.00_21',
                    '4-000000-0000000__20170630T012050.1Z_00.00_22',
                    '4-000000-0000000__20170630T012050.1Z_00.00_23',
                    '4-000000-0000000__20170630T012050.1Z_00.00_24',
                    '4-000000-0000000__20170630T012050.1Z_00.00_25',
                    '4-000000-0000000__20170630T012050.1Z_00.00_26',
                    '4-000000-0000000__20170630T012050.1Z_00.00_27',
                    '4-000000-0000000__20170630T012050.1Z_00.00_28',
                    '4-000000-0000000__20170630T012050.1Z_00.00_29',
                    '4-000000-0000000__20170630T012050.1Z_00.00_2',
                    '4-000000-0000000__20170630T012050.1Z_00.00_30',
                    '4-000000-0000000__20170630T012050.1Z_00.00_31',
                    '4-000000-0000000__20170630T012050.1Z_00.00_32',
                    '4-000000-0000000__20170630T012050.1Z_00.00_33',
                    '4-000000-0000000__20170630T012050.1Z_00.00_34',
                    '4-000000-0000000__20170630T012050.1Z_00.00_35',
                    '4-000000-0000000__20170630T012050.1Z_00.00_3',
                    '4-000000-0000000__20170630T012050.1Z_00.00_4',
                    '4-000000-0000000__20170630T012050.1Z_00.00_5',
                    '4-000000-0000000__20170630T012050.1Z_00.00_6',
                    '4-000000-0000000__20170630T012050.1Z_00.00_7',
                    '4-000000-0000000__20170630T012050.1Z_00.00_8',
                    '4-000000-0000000__20170630T012050.1Z_00.00_9']

    for cat_ in cats_partial:
        cats_final.append('{}{}.cat'.format(cat_part_1, cat_))

    cat_file = cats_final[cat_n]

    return cat_file


def check_source(i_alpha, i_delta, e_df):
    """

    :param i_alpha:
    :param i_delta:
    :param e_df:
    :return: e_df
    """
    prfs_d = extract_settings_sc3()

    e_df = e_df[e_df['ALPHA_J2000'] + prfs_d['tolerance'] > i_alpha]
    e_df = e_df[i_alpha > e_df['ALPHA_J2000'] - prfs_d['tolerance']]
    e_df = e_df[e_df['DELTA_J2000'] + prfs_d['tolerance'] > i_delta]
    e_df = e_df[i_delta > e_df['DELTA_J2000'] - prfs_d['tolerance']]

    return e_df


def load_sextractor_cats():
    """

    :return: cat_d
    """
    prfs_d = extract_settings_sc3()

    cats_number = 144  # todo hardcoded!
    cat_d = {}
    for cat_n in range(1, cats_number + 1, 1):
        cat_file = get_cat(cat_n)
        cat_data = fits.open('{}/{}'.format(prfs_d['fits_dir'], cat_file))

        ccd_df = Table(cat_data[2].data)
        print('CCD catalog {} to Pandas'.format(cat_n))
        cat_d[cat_n] = ccd_df.to_pandas()

    return cat_d


def main():
    source_list = []

    # input_catalog = read_csv('sc3_mer_10_starflag1.csv')
    # Load all catalogs
    cat_d = load_sextractor_cats()
    print(cat_d)
    """
    # Get boundaries for all catalogs
    # Look for sources
    # Create two dataframes (stars and galaxies)
    extracted_catalog = fits.open('test_cat.cat')
    full_db = Table(extracted_catalog[2].data)
    full_db = full_db.to_pandas()

    total = input_catalog['X_WORLD'].size
    for i, row in enumerate(input_catalog.itertuples(), 1):
        print('source number: {} - total number: {}'.format(i, total))
        i_alpha = row.X_WORLD
        i_delta = row.Y_WORLD

        e_df = check_source(i_alpha, i_delta, full_db)
        if e_df.empty is not True:
            if e_df['NUMBER'].size == 1:
                source_list.append(e_df['NUMBER'].iloc[0])

    sliced_df = full_db[full_db['NUMBER'].isin(source_list)]
    sliced_df = sliced_df[['X_WORLD', 'Y_WORLD', 'ERRA_WORLD',
                           'ERRB_WORLD', 'MAG_AUTO', 'MAGERR_AUTO']]

    return sliced_df
    """


if __name__ == "__main__":
    main()