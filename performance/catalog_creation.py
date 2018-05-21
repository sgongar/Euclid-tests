#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Creates a catalog from sextracted catalogs from single CCDs images.

output -> 1. Catalog is good enough. (?)

Versions:
- 0.1: Initial release. Split from check.py
       Recreated for ELViS analysis pipeline.
- 0.2: Single input catalogue and multiple input catalogues added.
- 0.3: Check method for catalogue creation added.

Information:
- cat: -> hdu_list catalogue
- data: -> Table formatted data
- df: -> dataframe formatted data

Todo:
    * Get out columns definition. Too long for a single function.
    * Creates units tests.
    * Improve variable nomenclature.
    * Explanations are still not so easy to understand.

*GNU Terry Pratchett*

"""
from astropy.io import fits
from astropy.table import Table
from pandas import read_csv

from misc import get_cats, extract_settings_elvis

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.2"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


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


def check_cats(prfs_d, fpa_cat):
    """

    :param prfs_d:
    :param fpa_cat:
    :return:
    """
    # Full catalogue data
    fpa_data = Table(fpa_cat[2].data)
    fpa_df = fpa_data.to_pandas()

    # Merge output
    dither = 1
    cat_files = get_cats(dither)
    cats_d = {}

    for cat_ in cat_files:
        cat_name = '{}/{}'.format(prfs_d['fits_dir'], cat_)
        cat = fits.open(cat_name)
        cat_data = Table(cat[2].data)
        cat_df = cat_data.to_pandas()

        cats_d[cat_] = cat_df

    total_sources = fpa_df['NUMBER'].size  # Gets the total amount of sources
    print('Total sources number {}'.format(total_sources))

    # Loops over the sources trying to find them in CCDs catalogues
    for idx_source, source_ in enumerate(fpa_df['NUMBER']):
        # print('Source {} of {}'.format(idx_source, total_sources))
        source_df = fpa_df[fpa_df['NUMBER'].isin([source_])]

        alpha = source_df['ALPHA_J2000'].iloc[0]
        delta = source_df['DELTA_J2000'].iloc[0]
        # print(alpha, delta)  # Test reasons
        for idx_cat, cat_ in enumerate(cats_d.keys()):
            # print('Checking catalogue number {}'.format(idx_cat))
            test_df = check_source(cats_d[cat_], alpha, delta)
            if test_df.empty is not True:
                # print('source detected')  # Test reasons
                total_sources -= 1
            else:
                pass

    # This final number represents how many sources of created catalogue are
    # present in CCDs catalogues
    print('Final sources number {}'.format(total_sources))


def single_cat(prfs_d):
    """

    :param prfs_d:
    :return:
    """

    columns = ['ALPHA_J2000', 'DELTA_J2000', 'pm', 'pa',
               'mag', 'mag', 'mag', 'mag']
    input_catalog = read_csv('{}/SSO_Cat.txt'.format(prfs_d['references']),
                             delim_whitespace=True, header=None, names=columns)
    input_catalog = input_catalog[['ALPHA_J2000', 'DELTA_J2000', 'pm', 'mag']]

    # Merge output
    reject_l = []

    cat_name = '{}/coadd.cat'.format(prfs_d['fits_dir'])
    coadd_cat = fits.open(cat_name)
    cat_data = Table(coadd_cat[2].data).to_pandas()

    print('Opening catalog -> {}'.format(cat_name))
    print('Total sources to be analysed - > {}'.format(cat_data['NUMBER'].size))
    for idx_source, source_ in enumerate(cat_data['NUMBER']):
        source_df = cat_data[cat_data['NUMBER'].isin([source_])]
        alpha = source_df['ALPHA_J2000'].iloc[0]
        delta = source_df['DELTA_J2000'].iloc[0]

        test_df = check_source(input_catalog, alpha, delta)

        if test_df.empty is not True:
            reject_l.append(source_)

    print('Removing SSOs from catalog')
    cat_data = cat_data[~cat_data['NUMBER'].isin(reject_l)]
    print('Removing galaxies from catalog')
    cat_data = cat_data[cat_data['CLASS_STAR'] > 0.95]

    c1 = fits.Column(name='NUMBER', format='1J', disp='I10',
                     array=cat_data['NUMBER'])
    # Object position along x
    c2 = fits.Column(name='X_IMAGE', format='1E', unit='pixel',
                     disp='F11.4', array=cat_data['X_IMAGE'])
    # Object position along y
    c3 = fits.Column(name='Y_IMAGE', format='1E', unit='pixel',
                     disp='F11.4', array=cat_data['Y_IMAGE'])
    # Variance along x
    c4 = fits.Column(name='X2_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_data['X2_IMAGE'])
    # Variance along y
    c5 = fits.Column(name='Y2_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_data['Y2_IMAGE'])
    # Covariance between x and y
    c6 = fits.Column(name='XY_IMAGE', format='1D', unit='pixel**2',
                     disp='E18.10', array=cat_data['XY_IMAGE'])
    # Isophotal area above Analysis threshold
    c7 = fits.Column(name='ISOAREA_IMAGE', format='1J', unit='pixel**2',
                     disp='I9', array=cat_data['ISOAREA_IMAGE'])
    # Background at centroid position
    c8 = fits.Column(name='BACKGROUND', format='1E', unit='count',
                     disp='G12.7', array=cat_data['BACKGROUND'])
    # Detection threshold above background
    c9 = fits.Column(name='THRESHOLD', format='1E', unit='count',
                     disp='G12.7', array=cat_data['THRESHOLD'])
    # Peak flux above background
    c10 = fits.Column(name='FLUX_MAX', format='1E', unit='count',
                      disp='G12.7', array=cat_data['FLUX_MAX'])
    # Profile RMS along major axis
    c11 = fits.Column(name='A_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_data['A_IMAGE'])
    # Profile RMS along minor axis
    c12 = fits.Column(name='B_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_data['B_IMAGE'])
    # Position angle(CCW / x)
    c13 = fits.Column(name='THETA_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_data['THETA_IMAGE'])
    # RMS position error along major axis
    c14 = fits.Column(name='ERRA_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_data['ERRA_IMAGE'])
    # RMS position error along minor axis
    c15 = fits.Column(name='ERRB_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_data['ERRB_IMAGE'])
    # Isophotal flux
    c16 = fits.Column(name='FLUX_ISO', format='1E', unit='count',
                      disp='G12.7', array=cat_data['FLUX_ISO'])
    # RMS error for isophotal flux
    c17 = fits.Column(name='FLUXERR_ISO', format='1E', unit='count',
                      disp='G12.7', array=cat_data['FLUXERR_ISO'])
    # Isophotal magnitude
    c18 = fits.Column(name='MAG_ISO', format='1E', unit='mag',
                      disp='F8.4', array=cat_data['MAG_ISO'])
    # RMS error for isophotal magnitude
    c19 = fits.Column(name='MAGERR_ISO', format='1E', unit='mag',
                      disp='F8.4', array=cat_data['MAGERR_ISO'])
    # Flux vector within fixed circular aperture(s)
    c20 = fits.Column(name='FLUX_APER', format='1E', unit='count',
                      disp='G12.7', array=cat_data['FLUX_APER'])
    # RMS error vector for aperture flux(es)
    c21 = fits.Column(name='FLUXERR_APER', format='1E', unit='count',
                      disp='G12.7', array=cat_data['FLUXERR_APER'])
    # Fixed aperture magnitude vector
    c22 = fits.Column(name='MAG_APER', format='1E', unit='mag',
                      disp='F8.4', array=cat_data['MAG_APER'])
    # RMS error vector for fixed aperture mag
    c23 = fits.Column(name='MAGERR_APER', format='1E', unit='mag',
                      disp='F8.4', array=cat_data['MAGERR_APER'])
    # Right ascension of barycenter (native)
    c24 = fits.Column(name='ALPHA_SKY', format='1D', unit='deg',
                      disp='F11.7', array=cat_data['ALPHA_SKY'])
    # Declination of barycenter (native)
    c25 = fits.Column(name='DELTA_SKY', format='1D', unit='deg',
                      disp='F11.7', array=cat_data['DELTA_SKY'])
    # Error ellipse position angle (CCW/x)
    c26 = fits.Column(name='ERRTHETA_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_data['ERRTHETA_IMAGE'])
    # Peak surface brightness above background
    c27 = fits.Column(name='MU_MAX', format='1E',
                      unit='mag * arcsec**(-2)', disp='F8.4',
                      array=cat_data['MU_MAX'])
    # FWHM assuming a gaussian core
    c28 = fits.Column(name='FWHM_IMAGE', format='1E', unit='pixel',
                      disp='F8.2', array=cat_data['FWHM_IMAGE'])
    # S / G classifier output
    c29 = fits.Column(name='CLASS_STAR', format='1E', disp='F6.3',
                      array=cat_data['CLASS_STAR'])
    # Fraction-of-light radii
    c30 = fits.Column(name='FLUX_RADIUS', format='1E', unit='pixel',
                      disp='F10.3', array=cat_data['FLUX_RADIUS'])
    # A_IMAGE/B_IMAGE
    c31 = fits.Column(name='ELONGATION', format='1E', disp='F8.3',
                      array=cat_data['ELONGATION'])
    # 1 - B_IMAGE/A_IMAGE
    c32 = fits.Column(name='ELLIPTICITY', format='1E', disp='F8.3',
                      array=cat_data['ELLIPTICITY'])
    # Cxx object ellipse parameter
    c33 = fits.Column(name='CXX_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_data['CXX_IMAGE'])
    # Cxy object ellipse parameter
    c34 = fits.Column(name='CXY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_data['CXY_IMAGE'])
    # Cyy object ellipse parameter
    c35 = fits.Column(name='CYY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='E15.7', array=cat_data['CYY_IMAGE'])
    # Cxx error ellipse parameter
    c36 = fits.Column(name='ERRCXX_IMAGE', format='1E',
                      unit='pixel**(-2)', disp='G12.7',
                      array=cat_data['ERRCXX_IMAGE'])
    # Cxy error ellipse parameter
    c37 = fits.Column(name='ERRCXY_IMAGE', format='1E',
                      unit='pixel**(-2)', disp='G12.7',
                      array=cat_data['ERRCXY_IMAGE'])
    # Cyy error ellipse parameter
    c38 = fits.Column(name='ERRCYY_IMAGE', format='1E', unit='pixel**(-2)',
                      disp='G12.7', array=cat_data['ERRCYY_IMAGE'])
    # Kron-like elliptical aperture magnitude
    c39 = fits.Column(name='MAG_AUTO', format='1E', unit='mag',
                      disp='F8.4', array=cat_data['MAG_AUTO'])
    # Windowed position estimate along x
    c40 = fits.Column(name='XWIN_IMAGE', format='1D', unit='pixel',
                      disp='F11.4', array=cat_data['XWIN_IMAGE'])
    # Windowed position estimate along y
    c41 = fits.Column(name='YWIN_IMAGE', format='1D', unit='pixel',
                      disp='F11.4', array=cat_data['YWIN_IMAGE'])
    # Flux within a Kron-like elliptical aperture
    c42 = fits.Column(name='FLUX_AUTO', format='1E', unit='count',
                      disp='G12.7', array=cat_data['FLUX_AUTO'])
    # RMS error for AUTO flux
    c43 = fits.Column(name='FLUXERR_AUTO', format='1E', unit='count',
                      disp='G12.7', array=cat_data['FLUXERR_AUTO'])
    # RMS error for AUTO magnitude
    c44 = fits.Column(name='MAGERR_AUTO', format='1E', unit='mag',
                      disp='F8.4', array=cat_data['MAGERR_AUTO'])
    # Gaussian-weighted SNR
    c45 = fits.Column(name='SNR_WIN', format='1E', disp='G10.4',
                      array=cat_data['SNR_WIN'])
    # Right ascension of barycenter (J2000)
    c46 = fits.Column(name='ALPHA_J2000', format='1D', unit='deg',
                      disp='F11.7', array=cat_data['ALPHA_J2000'])
    # Declination of barycenter (J2000)
    c47 = fits.Column(name='DELTA_J2000', format='1D', unit='deg',
                      disp='F11.7', array=cat_data['DELTA_J2000'])
    # Barycenter position along world x axis
    c48 = fits.Column(name='X_WORLD', format='1D', unit='deg', disp='E18.10',
                      array=cat_data['X_WORLD'])
    # Barycenter position along world y axis
    c49 = fits.Column(name='Y_WORLD', format='1D', unit='deg', disp='E18.10',
                      array=cat_data['Y_WORLD'])
    # Variance of position along X-WORLD (alpha)
    c50 = fits.Column(name='ERRX2_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_data['ERRX2_WORLD'])
    # Variance of position along Y-WORLD (delta)
    c51 = fits.Column(name='ERRY2_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_data['ERRY2_WORLD'])
    # Covariance of position X-WORLD/Y-WORLD
    c52 = fits.Column(name='ERRXY_WORLD', format='1D', unit='deg**2',
                      disp='E18.10', array=cat_data['ERRXY_WORLD'])
    # Windowed profile RMS along major axis
    c53 = fits.Column(name='AWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_data['AWIN_IMAGE'])
    # Windowed profile RMS along minor axis
    c54 = fits.Column(name='BWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.3', array=cat_data['BWIN_IMAGE'])
    # Windowed position angle (CCW/x)
    c55 = fits.Column(name='THETAWIN_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_data['THETAWIN_IMAGE'])
    # RMS windowed position error along major axis
    c56 = fits.Column(name='ERRAWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_data['ERRAWIN_IMAGE'])
    # RMS windowed position error along minor axis
    c57 = fits.Column(name='ERRBWIN_IMAGE', format='1E', unit='pixel',
                      disp='F9.5', array=cat_data['ERRBWIN_IMAGE'])
    # Windowed error ellipse position angle (CCW/x)
    c58 = fits.Column(name='ERRTHETAWIN_IMAGE', format='1E', unit='deg',
                      disp='F6.2', array=cat_data['ERRTHETAWIN_IMAGE'])
    # Extraction flags
    c59 = fits.Column(name='FLAGS', format='1I', disp='I3',
                      array=cat_data['FLAGS'])
    # FWHM assuming a gaussian core
    c60 = fits.Column(name='FWHM_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_data['FWHM_WORLD'])
    # World RMS position error along major axis
    c61 = fits.Column(name='ERRA_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_data['ERRA_WORLD'])
    # World RMS position error along minor axis
    c62 = fits.Column(name='ERRB_WORLD', format='1E', unit='deg',
                      disp='G12.7', array=cat_data['ERRB_WORLD'])

    col_defs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12,
                             c13, c14, c15, c16, c17, c18, c19, c20, c21, c22,
                             c23, c24, c25, c26, c27, c28, c29, c30, c31, c32,
                             c33, c34, c35, c36, c37, c38, c39, c40, c41, c42,
                             c43, c44, c45, c46, c47, c48, c49, c50, c51, c52,
                             c53, c54, c55, c56, c57, c58, c59, c60, c61, c62])

    tb_hdu = fits.BinTableHDU.from_columns(col_defs)

    coadd_cat[2] = tb_hdu
    coadd_cat[2].header['EXTNAME'] = 'LDAC_OBJECTS'

    newcat_name = '{}/coadd2.cat'.format(prfs_d['fits_dir'])
    coadd_cat.writeto(newcat_name, overwrite=True)

    return coadd_cat


if __name__ == "__main__":
    prfs_dict = extract_settings_elvis()

    # Create a FPA catalogue
    single_cat = single_cat(prfs_dict)
    # Checks FPA catalogue against CCD catalogues
    check_cats(prfs_dict, single_cat)
