#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

La idea es que compare los objetos extraidos en total con los que tengo.

Todo:
    * Improve log messages

"""

from pipeline.misc import extract_settings_sc3

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Gongora-Garcia"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def check_source(i_alpha, i_delta, e_df):
    tolerance = 0.0001389  # 0.5 arcsecond

    e_df = e_df[e_df['ALPHA_J2000'] + tolerance > i_alpha]
    e_df = e_df[i_alpha > e_df['ALPHA_J2000'] - tolerance]
    e_df = e_df[e_df['DELTA_J2000'] + tolerance > i_delta]
    e_df = e_df[i_delta > e_df['DELTA_J2000'] - tolerance]

    return e_df


def main():
    source_list = []

    input_catalog = read_csv('sc3_mer_10_starflag1.csv')
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


def create_data(sliced_df):
    c1 = fits.Column(name='X_WORLD', format='1E', unit='deg', 
                     disp='E18.10', array=sliced_df['X_WORLD'])
    c2 = fits.Column(name='Y_WORLD', format='1E', unit='deg',
                     disp='E18.10', array=sliced_df['Y_WORLD'])
    c3 = fits.Column(name = 'ERRA_WORLD', format = '1E',
                     unit = 'deg', disp = 'G12.7', array=sliced_df['ERRA_WORLD'])
    c4 = fits.Column(name = 'ERRB_WORLD', format = '1E',
                     unit = 'deg', disp = 'G12.7', array=sliced_df['ERRB_WORLD'])
    c5 = fits.Column(name = 'MAG_AUTO', format = '1E',
                     unit = 'mag', disp = 'F8.4', array=sliced_df['MAG_AUTO'])
    c6 = fits.Column(name = 'MAGERR_AUTO', format = '1E',
                     unit = 'mag', disp = 'F8.4', array=sliced_df['MAGERR_AUTO'])

    coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
    data = fits.BinTableHDU.from_columns(coldefs)

    return data


def create_cat(loc, data):
    hdu_list = fits.open(loc)

    # header = hdu_list[2].header
    hdu_list[2] = data
    hdu_list[2].header['EXTNAME'] = 'LDAC_OBJECTS'

    hdu_list.writeto('coadd_final.cat')


if __name__ == "__main__":
    """
    loc = 'test_cat.cat'

    # sliced_df = main()
    # sliced_df.to_csv('test_2.csv')
    sliced_df = read_csv('test_2.csv', index_col=0)
    data = create_data(sliced_df)
    create_cat(loc, data)
    """

    print(extract_settings_sc3)