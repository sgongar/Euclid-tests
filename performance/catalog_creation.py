#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a catalog from sextracted catalogs from single CCDs images.

Versions:
- 0.1: Initial release. Split from check.py
       Recreated for ELViS analysis pipeline.
- 0.2: Single input catalogue and multiple input catalogues added.
- 0.3: Check method for catalogue creation added.
- 0.4: Easiest way implemented.

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    * Get out columns definition. Too long for a single function.
    * Creates units tests.
    * Improve variable nomenclature.
    * Explanations are still not so easy to understand.
    * POO implementation?

*GNU Terry Pratchett*

"""
from math import hypot
from sys import stdout

from astropy.io import fits
from astropy.table import Table
from pandas import concat, DataFrame, read_csv

from misc import get_cat, get_cats, extract_settings_elvis

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.4"
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

    cat_stars = fits.open('{}/cat_stars.fits'.format(prfs_dict['references']))
    stars_data = Table(cat_stars[1].data)
    stars_df = stars_data.to_pandas()
    stars_idx = range(0, 28474, 1)
    stars_df['IDX'] = stars_idx
    inputs_d['stars'] = stars_df

    cat_galaxies = fits.open('{}/cat_galaxies.fits'.format(prfs_dict['references']))
    galaxies_data = Table(cat_galaxies[1].data)
    galaxies_df = galaxies_data.to_pandas()
    galaxies_idx = range(0, 143766, 1)
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
    catalogs = {}
    save = False

    cat_d = create_empty_catalog_dict()
    total_stars = inputs_d['stars']['IDX'].size
    # stdout.write('total stars {} \n'.format(total_stars))
    for idx, star in enumerate(inputs_d['stars']['IDX']):
        # stdout.write('star {} of {} \r'.format(idx, total_stars))
        # stdout.flush()
        stars_df = inputs_d['stars']
        source_df = stars_df[stars_df['IDX'].isin([star])]
        alpha = source_df['RA2000(Gaia)'].iloc[0]
        delta = source_df['DEC2000(Gaia)'].iloc[0]

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

        print(source_d)
        print(' ')
        print(' ')

    cat_df = DataFrame(cat_d)
    catalogs['stars'] = cat_df
    if save:
        cat_df.to_csv('stars.csv')

    cat_d = create_empty_catalog_dict()
    total_galaxies = inputs_d['galaxies']['IDX'].size
    stdout.write('total galaxies {} \n'.format(total_galaxies))
    for idx, galaxy in enumerate(inputs_d['galaxies']['IDX']):
        stdout.write('galaxy {} of {} \r'.format(idx, total_galaxies))
        stdout.flush()
        galaxies_df = inputs_d['galaxies']
        source_df = galaxies_df[galaxies_df['IDX'].isin([galaxy])]
        alpha = source_df['ra'].iloc[0]
        delta = source_df['dec'].iloc[0]

        dither_n = 0
        for dither in range(1, 5, 1):
            o_df = check_source(full_d[dither], alpha, delta)
            if o_df.empty is not True:
                dither_n += 1
                # Returns the index of the closest found source
                index = check_distance(o_df, alpha, delta)
                o_df = o_df.iloc[[index]]
                for key_ in cat_d.keys():
                    cat_d[key_].append(o_df[key_].iloc[0])

    cat_df = DataFrame(cat_d)
    catalogs['galaxies'] = cat_df
    if save:
        cat_df.to_csv('galaxies.csv')

    return catalogs


def write_catalog(catalogs):
    """

    :param catalogs:
    :return:
    """
    stars_df_data = catalogs['stars']
    # Stars catalogue creation
    test_cat_name = '{}/full_coadd.cat'.format(prfs_dict['fits_dir'])
    test_coadd_cat = fits.open(test_cat_name)
    # test_df_data = Table(test_coadd_cat[2].data).to_pandas()

    # Source number
    c1 = fits.Column(name='NUMBER', format='1J', disp='I10',
                     array=stars_df_data['NUMBER'])
    # Object position along x
    c2 = fits.Column(name='X_IMAGE', format='1E', unit='pixel',
                     disp='F11.4', array=stars_df_data['X_IMAGE'])
    # Object position along y
    c3 = fits.Column(name='Y_IMAGE', format='1E', unit='pixel',
                     disp='F11.4', array=stars_df_data['Y_IMAGE'])
    # Variance along x
    c4 = fits.Column(name='X2_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=stars_df_data['X2_IMAGE'])
    # Variance along y
    c5 = fits.Column(name='Y2_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=stars_df_data['Y2_IMAGE'])
    # Covariance between x and y
    c6 = fits.Column(name='XY_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=stars_df_data['XY_IMAGE'])
    # Isophotal area above Analysis threshold
    c7 = fits.Column(name='ISOAREA_IMAGE', format='1J', unit='pixel**2',
                     disp='I9', array=stars_df_data['ISOAREA_IMAGE'])
    # Background at centroid position
    c8 = fits.Column(name='BACKGROUND', format='1E', unit='count',
                     disp='G12.7', array=stars_df_data['BACKGROUND'])
    # Detection threshold above background
    c9 = fits.Column(name='THRESHOLD', format='1E', unit='count',
                     disp='G12.7', array=stars_df_data['THRESHOLD'])
    # Peak flux above background
    c10 = fits.Column(name='FLUX_MAX', format='1E', unit='count',
                      disp='G12.7', array=stars_df_data['FLUX_MAX'])
    # Profile RMS along major axis
    c11 = fits.Column(name='A_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=stars_df_data['A_IMAGE'])
    # Profile RMS along minor axis
    c12 = fits.Column(name='B_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=stars_df_data['B_IMAGE'])
    # Position angle(CCW / x)
    c13 = fits.Column(name='THETA_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=stars_df_data['THETA_IMAGE'])
    # RMS position error along major axis
    c14 = fits.Column(name='ERRA_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=stars_df_data['ERRA_IMAGE'])
    # RMS position error along minor axis
    c15 = fits.Column(name='ERRB_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=stars_df_data['ERRB_IMAGE'])
    # Isophotal flux
    c16 = fits.Column(name='FLUX_ISO', format='1E', unit='count',
                      disp='G12.7', array=stars_df_data['FLUX_ISO'])
    # RMS error for isophotal flux
    c17 = fits.Column(name='FLUXERR_ISO', format='1E', unit='count',
                      disp='G12.7', array=stars_df_data['FLUXERR_ISO'])
    # Isophotal magnitude
    c18 = fits.Column(name='MAG_ISO', format='1E', unit='mag',
                      disp='F8.4', array=stars_df_data['MAG_ISO'])
    # RMS error for isophotal magnitude
    c19 = fits.Column(name='MAGERR_ISO', format='1E', unit='mag',
                      disp='F8.4', array=stars_df_data['MAGERR_ISO'])
    # Flux vector within fixed circular aperture(s)
    c20 = fits.Column(name='FLUX_APER', format='1E', unit='count',
                      disp='G12.7', array=stars_df_data['FLUX_APER'])
    # RMS error vector for aperture flux(es)
    c21 = fits.Column(name='FLUXERR_APER', format='1E', unit='count',
                      disp='G12.7', array=stars_df_data['FLUXERR_APER'])
    # Fixed aperture magnitude vector
    c22 = fits.Column(name='MAG_APER', format='1E', unit='mag',
                      disp='F8.4', array=stars_df_data['MAG_APER'])
    # RMS error vector for fixed aperture mag
    c23 = fits.Column(name='MAGERR_APER', format='1E', unit='mag',
                      disp='F8.4', array=stars_df_data['MAGERR_APER'])
    # Right ascension of barycenter (native)
    c24 = fits.Column(name='ALPHA_SKY', format='1D', unit='deg',
                      disp='F11.7', array=stars_df_data['ALPHA_SKY'])
    # Declination of barycenter (native)
    c25 = fits.Column(name='DELTA_SKY', format='1D', unit='deg',
                      disp='F11.7', array=stars_df_data['DELTA_SKY'])
    # Error ellipse position angle (CCW/x)
    c26 = fits.Column(name='ERRTHETA_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=stars_df_data['ERRTHETA_IMAGE'])
    # Peak surface brightness above background
    c27 = fits.Column(name='MU_MAX', format='1E',
                      unit='mag * arcsec**(-2)', disp='F8.4',
                      array=stars_df_data['MU_MAX'])
    # FWHM assuming a gaussian core
    c28 = fits.Column(name='FWHM_IMAGE', format='1E', unit='pixel',
                      disp='F8.2', array=stars_df_data['FWHM_IMAGE'])
    # S / G classifier output
    c29 = fits.Column(name='CLASS_STAR', format='1E', disp='F6.3',
                      array=stars_df_data['CLASS_STAR'])
    # Fraction-of-light radii
    c30 = fits.Column(name='FLUX_RADIUS', format='1E', unit='pixel',
                      disp='F10.3', array=stars_df_data['FLUX_RADIUS'])
    # A_IMAGE/B_IMAGE
    c31 = fits.Column(name='ELONGATION', format='1E', disp='F8.3',
                      array=stars_df_data['ELONGATION'])
    # 1 - B_IMAGE/A_IMAGE
    c32 = fits.Column(name='ELLIPTICITY', format='1E', disp='F8.3',
                      array=stars_df_data['ELLIPTICITY'])
    # Cxx object ellipse parameter
    c33 = fits.Column(name='CXX_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=stars_df_data['CXX_IMAGE'])
    # Cxy object ellipse parameter
    c34 = fits.Column(name='CXY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=stars_df_data['CXY_IMAGE'])
    # Cyy object ellipse parameter
    c35 = fits.Column(name='CYY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=stars_df_data['CYY_IMAGE'])
    # Cxx error ellipse parameter
    c36 = fits.Column(name='ERRCXX_IMAGE', format='1E',
                      unit='pixel**(-2)', disp='G12.7',
                      array=stars_df_data['ERRCXX_IMAGE'])
    # Cxy error ellipse parameter
    c37 = fits.Column(name='ERRCXY_IMAGE', format='1E',
                      unit='pixel**(-2)', disp='G12.7',
                      array=stars_df_data['ERRCXY_IMAGE'])
    # Cyy error ellipse parameter
    c38 = fits.Column(name='ERRCYY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='G12.7', array=stars_df_data['ERRCYY_IMAGE'])
    # Kron-like elliptical aperture magnitude
    c39 = fits.Column(name='MAG_AUTO', format='1E', unit='mag',
                      disp='F8.4', array=stars_df_data['MAG_AUTO'])
    # Windowed position estimate along x
    c40 = fits.Column(name='XWIN_IMAGE', format='1D', unit='pixel',
                      disp='F11.4', array=stars_df_data['XWIN_IMAGE'])
    # Windowed position estimate along y
    c41 = fits.Column(name='YWIN_IMAGE', format='1D', unit='pixel',
                      disp='F11.4', array=stars_df_data['YWIN_IMAGE'])
    # Flux within a Kron-like elliptical aperture
    c42 = fits.Column(name='FLUX_AUTO', format='1E', unit='count',
                      disp='G12.7', array=stars_df_data['FLUX_AUTO'])
    # RMS error for AUTO flux
    c43 = fits.Column(name='FLUXERR_AUTO', format='1E', unit='count',
                      disp='G12.7', array=stars_df_data['FLUXERR_AUTO'])
    # RMS error for AUTO magnitude
    c44 = fits.Column(name='MAGERR_AUTO', format='1E', unit='mag',
                      disp='F8.4', array=stars_df_data['MAGERR_AUTO'])
    # Gaussian-weighted SNR
    c45 = fits.Column(name='SNR_WIN', format='1E', disp='G10.4',
                      array=stars_df_data['SNR_WIN'])
    # Right ascension of barycenter (J2000)
    c46 = fits.Column(name='ALPHA_J2000', format='1D', unit='deg',
                      disp='F11.7', array=stars_df_data['ALPHA_J2000'])
    # Declination of barycenter (J2000)
    c47 = fits.Column(name='DELTA_J2000', format='1D', unit='deg',
                      disp='F11.7', array=stars_df_data['DELTA_J2000'])
    # Barycenter position along world x axis
    c48 = fits.Column(name='X_WORLD', format='1D', unit='deg', disp='E18.10',
                      array=stars_df_data['X_WORLD'])
    # Barycenter position along world y axis
    c49 = fits.Column(name='Y_WORLD', format='1D', unit='deg', disp='E18.10',
                      array=stars_df_data['Y_WORLD'])
    # Variance of position along X-WORLD (alpha)
    c50 = fits.Column(name='ERRX2_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=stars_df_data['ERRX2_WORLD'])
    # Variance of position along Y-WORLD (delta)
    c51 = fits.Column(name='ERRY2_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=stars_df_data['ERRY2_WORLD'])
    # Covariance of position X-WORLD/Y-WORLD
    c52 = fits.Column(name='ERRXY_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=stars_df_data['ERRXY_WORLD'])
    # Windowed profile RMS along major axis
    c53 = fits.Column(name='AWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=stars_df_data['AWIN_IMAGE'])
    # Windowed profile RMS along minor axis
    c54 = fits.Column(name='BWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=stars_df_data['BWIN_IMAGE'])
    # Windowed position angle (CCW/x)
    c55 = fits.Column(name='THETAWIN_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=stars_df_data['THETAWIN_IMAGE'])
    # RMS windowed position error along major axis
    c56 = fits.Column(name='ERRAWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=stars_df_data['ERRAWIN_IMAGE'])
    # RMS windowed position error along minor axis
    c57 = fits.Column(name='ERRBWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=stars_df_data['ERRBWIN_IMAGE'])
    # Windowed error ellipse position angle (CCW/x)
    c58 = fits.Column(name='ERRTHETAWIN_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=stars_df_data['ERRTHETAWIN_IMAGE'])
    # Extraction flags
    c59 = fits.Column(name='FLAGS', format='1I', disp='I3',
                      array=stars_df_data['FLAGS'])
    # FWHM assuming a gaussian core
    c60 = fits.Column(name='FWHM_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=stars_df_data['FWHM_WORLD'])
    # World RMS position error along major axis
    c61 = fits.Column(name='ERRA_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=stars_df_data['ERRA_WORLD'])
    # World RMS position error along minor axis
    c62 = fits.Column(name='ERRB_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=stars_df_data['ERRB_WORLD'])

    col_defs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                             c13, c14, c15, c16, c17, c18, c19, c20, c21, c22,
                             c23, c24, c25, c26, c27, c28, c29, c30, c31, c32,
                             c33, c34, c35, c36, c37, c38, c39, c40, c41, c42,
                             c43, c44, c45, c46, c47, c48, c49, c50, c51, c52,
                             c53, c54, c55, c56, c57, c58, c59, c60, c61, c62])

    tb_hdu = fits.BinTableHDU.from_columns(col_defs)
    #
    #     test_coadd_cat[2] = tb_hdu
    #     test_coadd_cat[2].header['EXTNAME'] = 'LDAC_OBJECTS'
    #
    #     newcat_name = '{}/stars_catalogue.cat'.format(prfs_dict['fits_dir'])
    #     test_coadd_cat.writeto(newcat_name, overwrite=True)


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    # Create a FPA catalogue
    # single_cat = single_cat(prfs_dict)
    # Checks FPA catalogue against CCD catalogues
    # check_cats(prfs_dict, single_cat)

    # test
    # cat_test()
    catalogs = create_catalog()
    write_catalog(catalogs)





























# def single_cat(prfs_d):
#     """
#
#     :param prfs_d:
#     :return:
#     """
#
#     columns = ['ALPHA_J2000', 'DELTA_J2000', 'pm', 'pa',
#                'mag', 'mag', 'mag', 'mag']
#     input_catalog = read_csv('{}/SSO_Cat.txt'.format(prfs_d['references']),
#                              delim_whitespace=True, header=None, names=columns)
#     input_catalog = input_catalog[['ALPHA_J2000', 'DELTA_J2000', 'pm', 'mag']]
#
#     # Merge output
#     reject_l = []
#

#
# def cat_test():
#     """
#
#     :return:
#     """
#     test_cat_name = '{}/full_coadd.cat'.format(prfs_dict['fits_dir'])
#     test_coadd_cat = fits.open(test_cat_name)
#     test_df_data = Table(test_coadd_cat[2].data).to_pandas()
#     test_df_data['ID'] = range(0, 145417, 1)
#
#     reject_ids = []
#     total_ids = 145417
#     unique_sources = range(0, 145417, 1)
#
#     sub_list_size = total_ids / prfs_dict['cores_number']
#
#     sub_list_l = []
#     for idx_sub_list in range(0, prfs_dict['cores_number'], 1):
#         if idx_sub_list != (prfs_dict['cores_number'] - 1):
#             idx_down = sub_list_size * idx_sub_list
#             idx_up = sub_list_size * (idx_sub_list + 1)
#             sub_list_l.append(unique_sources[idx_down:idx_up])
#         else:
#             idx_down = sub_list_size * idx_sub_list
#             sub_list_l.append(unique_sources[idx_down:])
#
#     areas_j = []
#     for idx_l in range(0, prfs_dict['cores_number'], 1):
#         areas_p = Process(target=cat_test_thread,
#                           args=(idx_l, sub_list_l[idx_l], test_df_data,))
#         areas_j.append(areas_p)
#         areas_p.start()
#
#     active_areas = list([job.is_alive() for job in areas_j])
#     while True in active_areas:
#         active_areas = list([job.is_alive() for job in areas_j])
#         pass
#
#
# def cat_test_thread(idx, sub_list, test_df_data):
#     """
#
#     :param idx:
#     :param sub_list:
#     :param test_df_data:
#     :return:
#     """
#     reject_ids = []
#     for id_ in sub_list:
#         source_df = test_df_data[test_df_data['ID'].isin([id_])]
#         alpha = source_df['ALPHA_J2000'].iloc[0]
#         delta = source_df['DELTA_J2000'].iloc[0]
#         #
#         test_df = check_source(test_df_data, alpha, delta)
#         if int(test_df['NUMBER'].size) != 1:
#             reject_ids.append(id_)
#
#     print(reject_ids)
#     with open('rejected_{}'.format(idx), 'wb') as f:
#         pickle.dump(reject_ids, f)
#
#
# def create_catalog():
#     """
#
#     :return:
#     """
#     rejects_l = []
#

#
#     # Galaxies catalogue creation
#     test_cat_name = '{}/full_coadd.cat'.format(prfs_dict['fits_dir'])
#     test_coadd_cat = fits.open(test_cat_name)
#     test_df_data = Table(test_coadd_cat[2].data).to_pandas()
#     test_df_data['ID'] = range(0, 145408, 1)
#
#     print('Removing duplicates')
#     test_df_data = test_df_data[~test_df_data['NUMBER'].isin(rejects_l)]
#     print('Removing galaxies from catalog')
#     galaxies_df_data = test_df_data[test_df_data['CLASS_STAR'] < 0.95]
#
#     # Source number
#     c1 = fits.Column(name='NUMBER', format='1J', disp='I10',
#                      array=galaxies_df_data['NUMBER'])
#     # Object position along x
#     c2 = fits.Column(name='X_IMAGE', format='1E', unit='pixel',
#                      disp='F11.4', array=galaxies_df_data['X_IMAGE'])
#     # Object position along y
#     c3 = fits.Column(name='Y_IMAGE', format='1E', unit='pixel',
#                      disp='F11.4', array=galaxies_df_data['Y_IMAGE'])
#     # Variance along x
#     c4 = fits.Column(name='X2_IMAGE', format='1D', unit='pixel**2',
#                      disp='E18.10', array=galaxies_df_data['X2_IMAGE'])
#     # Variance along y
#     c5 = fits.Column(name='Y2_IMAGE', format='1D', unit='pixel**2',
#                      disp='E18.10', array=galaxies_df_data['Y2_IMAGE'])
#     # Covariance between x and y
#     c6 = fits.Column(name='XY_IMAGE', format='1D', unit='pixel**2',
#                      disp='E18.10', array=galaxies_df_data['XY_IMAGE'])
#     # Isophotal area above Analysis threshold
#     c7 = fits.Column(name='ISOAREA_IMAGE', format='1J', unit='pixel**2',
#                      disp='I9', array=galaxies_df_data['ISOAREA_IMAGE'])
#     # Background at centroid position
#     c8 = fits.Column(name='BACKGROUND', format='1E', unit='count',
#                      disp='G12.7', array=galaxies_df_data['BACKGROUND'])
#     # Detection threshold above background
#     c9 = fits.Column(name='THRESHOLD', format='1E', unit='count',
#                      disp='G12.7', array=galaxies_df_data['THRESHOLD'])
#     # Peak flux above background
#     c10 = fits.Column(name='FLUX_MAX', format='1E', unit='count',
#                       disp='G12.7', array=galaxies_df_data['FLUX_MAX'])
#     # Profile RMS along major axis
#     c11 = fits.Column(name='A_IMAGE', format='1E', unit='pixel',
#                       disp='F9.3', array=galaxies_df_data['A_IMAGE'])
#     # Profile RMS along minor axis
#     c12 = fits.Column(name='B_IMAGE', format='1E', unit='pixel',
#                       disp='F9.3', array=galaxies_df_data['B_IMAGE'])
#     # Position angle(CCW / x)
#     c13 = fits.Column(name='THETA_IMAGE', format='1E', unit='deg',
#                       disp='F6.2', array=galaxies_df_data['THETA_IMAGE'])
#     # RMS position error along major axis
#     c14 = fits.Column(name='ERRA_IMAGE', format='1E', unit='pixel',
#                       disp='F9.5', array=galaxies_df_data['ERRA_IMAGE'])
#     # RMS position error along minor axis
#     c15 = fits.Column(name='ERRB_IMAGE', format='1E', unit='pixel',
#                       disp='F9.5', array=galaxies_df_data['ERRB_IMAGE'])
#     # Isophotal flux
#     c16 = fits.Column(name='FLUX_ISO', format='1E', unit='count',
#                       disp='G12.7', array=galaxies_df_data['FLUX_ISO'])
#     # RMS error for isophotal flux
#     c17 = fits.Column(name='FLUXERR_ISO', format='1E', unit='count',
#                       disp='G12.7', array=galaxies_df_data['FLUXERR_ISO'])
#     # Isophotal magnitude
#     c18 = fits.Column(name='MAG_ISO', format='1E', unit='mag',
#                       disp='F8.4', array=galaxies_df_data['MAG_ISO'])
#     # RMS error for isophotal magnitude
#     c19 = fits.Column(name='MAGERR_ISO', format='1E', unit='mag',
#                       disp='F8.4', array=galaxies_df_data['MAGERR_ISO'])
#     # Flux vector within fixed circular aperture(s)
#     c20 = fits.Column(name='FLUX_APER', format='1E', unit='count',
#                       disp='G12.7', array=galaxies_df_data['FLUX_APER'])
#     # RMS error vector for aperture flux(es)
#     c21 = fits.Column(name='FLUXERR_APER', format='1E', unit='count',
#                       disp='G12.7', array=galaxies_df_data['FLUXERR_APER'])
#     # Fixed aperture magnitude vector
#     c22 = fits.Column(name='MAG_APER', format='1E', unit='mag',
#                       disp='F8.4', array=galaxies_df_data['MAG_APER'])
#     # RMS error vector for fixed aperture mag
#     c23 = fits.Column(name='MAGERR_APER', format='1E', unit='mag',
#                       disp='F8.4', array=galaxies_df_data['MAGERR_APER'])
#     # Right ascension of barycenter (native)
#     c24 = fits.Column(name='ALPHA_SKY', format='1D', unit='deg',
#                       disp='F11.7', array=galaxies_df_data['ALPHA_SKY'])
#     # Declination of barycenter (native)
#     c25 = fits.Column(name='DELTA_SKY', format='1D', unit='deg',
#                       disp='F11.7', array=galaxies_df_data['DELTA_SKY'])
#     # Error ellipse position angle (CCW/x)
#     c26 = fits.Column(name='ERRTHETA_IMAGE', format='1E', unit='deg',
#                       disp='F6.2', array=galaxies_df_data['ERRTHETA_IMAGE'])
#     # Peak surface brightness above background
#     c27 = fits.Column(name='MU_MAX', format='1E',
#                       unit='mag * arcsec**(-2)', disp='F8.4',
#                       array=galaxies_df_data['MU_MAX'])
#     # FWHM assuming a gaussian core
#     c28 = fits.Column(name='FWHM_IMAGE', format='1E', unit='pixel',
#                       disp='F8.2', array=galaxies_df_data['FWHM_IMAGE'])
#     # S / G classifier output
#     c29 = fits.Column(name='CLASS_STAR', format='1E', disp='F6.3',
#                       array=galaxies_df_data['CLASS_STAR'])
#     # Fraction-of-light radii
#     c30 = fits.Column(name='FLUX_RADIUS', format='1E', unit='pixel',
#                       disp='F10.3', array=galaxies_df_data['FLUX_RADIUS'])
#     # A_IMAGE/B_IMAGE
#     c31 = fits.Column(name='ELONGATION', format='1E', disp='F8.3',
#                       array=galaxies_df_data['ELONGATION'])
#     # 1 - B_IMAGE/A_IMAGE
#     c32 = fits.Column(name='ELLIPTICITY', format='1E', disp='F8.3',
#                       array=galaxies_df_data['ELLIPTICITY'])
#     # Cxx object ellipse parameter
#     c33 = fits.Column(name='CXX_IMAGE', format='1E', unit='pixel**(-2)',
#                       disp='E15.7', array=galaxies_df_data['CXX_IMAGE'])
#     # Cxy object ellipse parameter
#     c34 = fits.Column(name='CXY_IMAGE', format='1E', unit='pixel**(-2)',
#                       disp='E15.7', array=galaxies_df_data['CXY_IMAGE'])
#     # Cyy object ellipse parameter
#     c35 = fits.Column(name='CYY_IMAGE', format='1E', unit='pixel**(-2)',
#                       disp='E15.7', array=galaxies_df_data['CYY_IMAGE'])
#     # Cxx error ellipse parameter
#     c36 = fits.Column(name='ERRCXX_IMAGE', format='1E',
#                       unit='pixel**(-2)', disp='G12.7',
#                       array=galaxies_df_data['ERRCXX_IMAGE'])
#     # Cxy error ellipse parameter
#     c37 = fits.Column(name='ERRCXY_IMAGE', format='1E',
#                       unit='pixel**(-2)', disp='G12.7',
#                       array=galaxies_df_data['ERRCXY_IMAGE'])
#     # Cyy error ellipse parameter
#     c38 = fits.Column(name='ERRCYY_IMAGE', format='1E', unit='pixel**(-2)',
#                       disp='G12.7', array=galaxies_df_data['ERRCYY_IMAGE'])
#     # Kron-like elliptical aperture magnitude
#     c39 = fits.Column(name='MAG_AUTO', format='1E', unit='mag',
#                       disp='F8.4', array=galaxies_df_data['MAG_AUTO'])
#     # Windowed position estimate along x
#     c40 = fits.Column(name='XWIN_IMAGE', format='1D', unit='pixel',
#                       disp='F11.4', array=galaxies_df_data['XWIN_IMAGE'])
#     # Windowed position estimate along y
#     c41 = fits.Column(name='YWIN_IMAGE', format='1D', unit='pixel',
#                       disp='F11.4', array=galaxies_df_data['YWIN_IMAGE'])
#     # Flux within a Kron-like elliptical aperture
#     c42 = fits.Column(name='FLUX_AUTO', format='1E', unit='count',
#                       disp='G12.7', array=galaxies_df_data['FLUX_AUTO'])
#     # RMS error for AUTO flux
#     c43 = fits.Column(name='FLUXERR_AUTO', format='1E', unit='count',
#                       disp='G12.7', array=galaxies_df_data['FLUXERR_AUTO'])
#     # RMS error for AUTO magnitude
#     c44 = fits.Column(name='MAGERR_AUTO', format='1E', unit='mag',
#                       disp='F8.4', array=galaxies_df_data['MAGERR_AUTO'])
#     # Gaussian-weighted SNR
#     c45 = fits.Column(name='SNR_WIN', format='1E', disp='G10.4',
#                       array=galaxies_df_data['SNR_WIN'])
#     # Right ascension of barycenter (J2000)
#     c46 = fits.Column(name='ALPHA_J2000', format='1D', unit='deg',
#                       disp='F11.7', array=galaxies_df_data['ALPHA_J2000'])
#     # Declination of barycenter (J2000)
#     c47 = fits.Column(name='DELTA_J2000', format='1D', unit='deg',
#                       disp='F11.7', array=galaxies_df_data['DELTA_J2000'])
#     # Barycenter position along world x axis
#     c48 = fits.Column(name='X_WORLD', format='1D', unit='deg', disp='E18.10',
#                       array=galaxies_df_data['X_WORLD'])
#     # Barycenter position along world y axis
#     c49 = fits.Column(name='Y_WORLD', format='1D', unit='deg', disp='E18.10',
#                       array=galaxies_df_data['Y_WORLD'])
#     # Variance of position along X-WORLD (alpha)
#     c50 = fits.Column(name='ERRX2_WORLD', format='1D', unit='deg**2',
#                       disp='E18.10', array=galaxies_df_data['ERRX2_WORLD'])
#     # Variance of position along Y-WORLD (delta)
#     c51 = fits.Column(name='ERRY2_WORLD', format='1D', unit='deg**2',
#                       disp='E18.10', array=galaxies_df_data['ERRY2_WORLD'])
#     # Covariance of position X-WORLD/Y-WORLD
#     c52 = fits.Column(name='ERRXY_WORLD', format='1D', unit='deg**2',
#                       disp='E18.10', array=galaxies_df_data['ERRXY_WORLD'])
#     # Windowed profile RMS along major axis
#     c53 = fits.Column(name='AWIN_IMAGE', format='1E', unit='pixel',
#                       disp='F9.3', array=galaxies_df_data['AWIN_IMAGE'])
#     # Windowed profile RMS along minor axis
#     c54 = fits.Column(name='BWIN_IMAGE', format='1E', unit='pixel',
#                       disp='F9.3', array=galaxies_df_data['BWIN_IMAGE'])
#     # Windowed position angle (CCW/x)
#     c55 = fits.Column(name='THETAWIN_IMAGE', format='1E', unit='deg',
#                       disp='F6.2', array=galaxies_df_data['THETAWIN_IMAGE'])
#     # RMS windowed position error along major axis
#     c56 = fits.Column(name='ERRAWIN_IMAGE', format='1E', unit='pixel',
#                       disp='F9.5', array=galaxies_df_data['ERRAWIN_IMAGE'])
#     # RMS windowed position error along minor axis
#     c57 = fits.Column(name='ERRBWIN_IMAGE', format='1E', unit='pixel',
#                       disp='F9.5', array=galaxies_df_data['ERRBWIN_IMAGE'])
#     # Windowed error ellipse position angle (CCW/x)
#     c58 = fits.Column(name='ERRTHETAWIN_IMAGE', format='1E', unit='deg',
#                       disp='F6.2', array=galaxies_df_data['ERRTHETAWIN_IMAGE'])
#     # Extraction flags
#     c59 = fits.Column(name='FLAGS', format='1I', disp='I3',
#                       array=galaxies_df_data['FLAGS'])
#     # FWHM assuming a gaussian core
#     c60 = fits.Column(name='FWHM_WORLD', format='1E', unit='deg',
#                       disp='G12.7', array=galaxies_df_data['FWHM_WORLD'])
#     # World RMS position error along major axis
#     c61 = fits.Column(name='ERRA_WORLD', format='1E', unit='deg',
#                       disp='G12.7', array=galaxies_df_data['ERRA_WORLD'])
#     # World RMS position error along minor axis
#     c62 = fits.Column(name='ERRB_WORLD', format='1E', unit='deg',
#                       disp='G12.7', array=galaxies_df_data['ERRB_WORLD'])
#
#     col_defs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
#                              c13, c14, c15, c16, c17, c18, c19, c20, c21, c22,
#                              c23, c24, c25, c26, c27, c28, c29, c30, c31, c32,
#                              c33, c34, c35, c36, c37, c38, c39, c40, c41, c42,
#                              c43, c44, c45, c46, c47, c48, c49, c50, c51, c52,
#                              c53, c54, c55, c56, c57, c58, c59, c60, c61, c62])
#
#     tb_hdu = fits.BinTableHDU.from_columns(col_defs)
#
#     test_coadd_cat[2] = tb_hdu
#     test_coadd_cat[2].header['EXTNAME'] = 'LDAC_OBJECTS'
#
#     newcat_name = '{}/galaxies_catalogue.cat'.format(prfs_dict['fits_dir'])
#     test_coadd_cat.writeto(newcat_name, overwrite=True)
#
