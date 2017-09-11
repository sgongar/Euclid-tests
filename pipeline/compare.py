#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Versions:
*

In order to improve legibilidad algunas abreviaturas han sido creadas
c = catalog
n = name
d = dictionary
loc = location

Todo:
    * Improve log messages
    * Improve docstring
    * Create versions history
"""

from multiprocessing import Process
from os import remove

from astropy.io import fits
from astropy.table import Table
from pandas import concat, DataFrame, read_csv

from cats_management import cut_catalog, get_output_catalog
from misc import extract_settings, get_fits, check_distance

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
"""
__license__ = "GPL"
"""
__version__ = "0.2"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class Compare:

    def __init__(self):
        """ __init__ creates the basic variables for our study
            it takes nothing from outside

        """
        prfs_d = extract_settings()

        # Hardcoded
        folder_sex = '2_1.35_1.35_0.1_4'
        folder_scmp = '150_1.2_5_0.033'

        # Creates a catalog
        self.merge_catalog(prfs_d, folder_sex, folder_scmp)
        # Perfoms an analysis
        # self.perform_analysis(prfs_d, folder_sex, folder_scmp)

    def perform_analysis(self, prfs_d, folder_sex, folder_scmp):
        """

        @param prfs_d:
        @param folder_sex:
        @param folder_scmp:

        @return True: if everything is alright
        """
        cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'], folder_sex,
                                     folder_scmp)

        single = False

        fits_files = get_fits(unique=False)

        for idx in range(0, len(fits_files), 5):
            compare_j = []
            for proc in range(0, 5, 1):
                idx_p = proc + idx
                if idx_proc < len(fits_files):
                    if single:
                        full_c = '{}/f_{}_1.cat'.format(cats_dir,
                                                        fits_files[idx_p][2:-5])
                    else:
                        pass

                    fits_n = '{}/{}/{}.cat'.format(prfs_d['fits_dir'],
                                                   folder_sex,
                                                   fits_files[idx_proc][:-5]) 

                    compare_p = Process(target=self.perform_analysis_thread,
                                        args=(full_c, fits_n, idx_proc,
                                              prfs_d,))
                    compare_j.append(compare_p)
                    compare_p.start()

            active_compare = list([job.is_alive() for job in compare_j])
            while True in active_compare:
                active_compare = list([job.is_alive() for job in compare_j])
                pass

        self.merge_stats(fits_files, prfs_d)

    def perform_analysis_thread(self, full_c, fits_n, idx, prfs_d):
        """

        @param full_c:
        @param fits_n:
        @param idx:
        """
        # Set-up indexes for this particular case
        idx_detected = 0
        idx_repeated = 0
        idx_lost = 0
        # Set-up dictionaries for this particular case
        stats_d = {}
        sources_d = {}
        stats_d, sources_d = self.populate_dict(stats_d, sources_d)

        stats_d['CCD'].append(fits_n[-12:-4])
        stats_d['dither'].append(fits_n[-5:-4])

        # Opens catalog file
        cat_file = fits.open(full_c)
        cat_data = Table(cat_file[2].data)
        cat_table = cat_data.to_pandas()
        # Removing zero catalog references
        cat_table = cat_table.loc[~cat_table['CATALOG_NUMBER'].isin([0])]

        stats_d['total'].append(len(cat_table['SOURCE_NUMBER'].tolist()))

        fits_file = fits.open(fits_n)
        fits_data = Table(fits_file[2].data)
        fits_table = fits_data.to_pandas()

        # All variables starting by 'cat' are referred to scamp output catalogs
        # All variables starting by 'fits' are referred to sextractor catalogs
        for cat_source in cat_table['SOURCE_NUMBER'].tolist():
            # Reset variables associated to each input source
            distances_cache = []
            duplicated_cache = []
            i_alpha_cache = []
            i_delta_cache = []
            o_alpha_cache = []
            o_delta_cache = []

            # Gets input data associated to selectes source
            cat_t = cat_table[cat_table['SOURCE_NUMBER'].isin([cat_source])]
            cat_ra = float(cat_t['ALPHA_J2000'])
            cat_dec = float(cat_t['DELTA_J2000'])

            # In order to improve comparation speed output catalog is reduced
            margins = [[cat_ra, 'ALPHA_J2000'], [cat_dec, 'DELTA_J2000']]
            margin = 2 * prfs_d['tolerance']
            fits_table_cut = cut_catalog(fits_table, margins, margin)

            # Takes a look over output sources
            # Each source has only one reference
            for fits_source in fits_table_cut['NUMBER'].tolist():
                mask = fits_table_cut['NUMBER'].isin([fits_source])
                fits_t = fits_table_cut[mask]
                fits_ra = float(fits_t['ALPHA_J2000'])
                fits_dec = float(fits_t['DELTA_J2000'])
                # Compare cat_ra, cat_dec against fits_ra, fits_dec
                close, distance = check_distance(cat_ra, fits_ra,
                                                 cat_dec, fits_dec,
                                                 prfs_d['tolerance'])
                if close:
                    close_flag = True
                    distance = distance * 3600  # From degrees to seconds
                    distances_cache.append(distance)
                    i_alpha_cache.append(cat_ra)
                    i_delta_cache.append(cat_dec)
                    o_alpha_cache.append(fits_ra)
                    o_delta_cache.append(fits_dec)

            if len(distances_cache) > 1 and close_flag == True:
                for idx_cache, distance_ in enumerate(distances_cache):
                    sources_d['CCD'].append(fits_n[-12:-4])
                    sources_d['dither'].append(fits_n[-5:-4])
                    sources_d['distance'].append(distance_)
                    sources_d['duplicated'].append(True)
                    sources_d['i_ALPHA_J2000'].append(i_alpha_cache[idx_cache])
                    sources_d['i_DELTA_J2000'].append(i_delta_cache[idx_cache])
                    sources_d['o_ALPHA_J2000'].append(o_alpha_cache[idx_cache])
                    sources_d['o_DELTA_J2000'].append(o_delta_cache[idx_cache])
                idx_repeated += 1
            elif len(distances_cache) == 1 and close_flag == True:
                for idx_cache, distance_ in enumerate(distances_cache):
                    sources_d['CCD'].append(fits_n[-12:-4])
                    sources_d['dither'].append(fits_n[-5:-4])
                    sources_d['distance'].append(distance_)
                    sources_d['duplicated'].append(False)
                    sources_d['i_ALPHA_J2000'].append(i_alpha_cache[idx_cache])
                    sources_d['i_DELTA_J2000'].append(i_delta_cache[idx_cache])
                    sources_d['o_ALPHA_J2000'].append(o_alpha_cache[idx_cache])
                    sources_d['o_DELTA_J2000'].append(o_delta_cache[idx_cache])
                idx_detected += 1
            elif len(distances_cache) == 0:
                idx_lost += 1
            else:
                raise Exception
        
        stats_d['detected'].append(idx_detected)
        stats_d['repeated'].append(idx_repeated)
        stats_d['lost'].append(idx_lost)

        # Creates a DataFrame from stats dictionary
        stats_df = DataFrame(stats_d)
        stats_df.to_csv('{}/{}.csv'.format(prfs_d['tmp_out'], fits_n[-12:-4]),
                        columns=['CCD', 'dither', 'total', 'detected',
                                 'repeated', 'lost'])

        # Creates a DataFrame from sources dictionary
        sources_df = DataFrame(sources_d)
        sources_df.to_csv('{}/sources_{}.csv'.format(prfs_d['tmp_out'],
                                                     fits_n[-12:-4]),
                          columns=['CCD', 'dither', 'distance', 'duplicated',
                                   'i_ALPHA_J2000', 'i_DELTA_J2000',
                                   'o_ALPHA_J2000', 'o_DELTA_J2000'])

    def populate_dict(self, stats_d, sources_d):
        """ populates dictionaries with selected keys

        @param stats_d:
        @param sources_d:

        @retun stats_d, sources_d
        """
        stats_d = {'CCD': [], 'dither': [], 'total': [], 'detected': [],
                   'repeated': [], 'lost': []}

        sources_d = {'CCD': [], 'dither': [], 'distance': [], 'duplicated': [],
                     'i_ALPHA_J2000': [], 'i_DELTA_J2000': [],
                     'o_ALPHA_J2000': [], 'o_DELTA_J2000': []}

        return stats_d, sources_d

    def merge_catalog(self, prfs_d, folder_sex, folder_scmp):
        """

        @param prfs_d:
        @param folder_sex:
        @param folder_scmp:

        @return cat_n: a reference in string format.
        """
        cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'], folder_sex,
                                     folder_scmp)
        fits_files = get_fits(unique=False)
        
        # Using a mutable dictionary instead a fixed lists will allow us
        # to change catalog columns without a problem
        cat_dict = {}
        for idx in range(0, len(fits_files), 1):
            cat_loc = '{}/{}/{}.cat'.format(prfs_d['fits_dir'], folder_sex,
                                            fits_files[idx][:-5])
            # Gets data from CCD catalog
            cat_data = get_output_catalog(prfs_d, cat_loc)
            cat_n = fits_files[idx][:-5]
            cat_dict[cat_n] = {}  # Creates a dict for selected CCD
            for column_ in cat_data.columns:
                cat_dict[cat_n][column_] = cat_data[column_].tolist()

        cat_final = {}
        for ccd_key in cat_dict:
            for param_key in cat_dict[ccd_key].keys():
                for value_ in cat_dict[ccd_key][param_key]:
                    try:
                        cat_final[param_key].append(value_)
                    except KeyError:
                        cat_final[param_key] = []
                        cat_final[param_key].append(value_)
        
        cat_df = DataFrame(cat_final)
        self.rewrite_catalog(cat_df)
        # cat_df.to_csv('test.csv')

    def rewrite_catalog(self, cat_df):
        """ writes a new catalogue only fulfilled with stars

        @param logger: a logger object
        @param prfs_d:
        @param catalogue_table:

        @return True: if everything goes alright
        """
        c1 = fits.Column(name='NUMBER', format='1J', disp='I10',
                         array=cat_df['NUMBER'])
        c2 = fits.Column(name='X_IMAGE', format='1E', unit='pixel',
                         disp='F11.4', array=cat_df['X_IMAGE'])
        c3 = fits.Column(name='Y_IMAGE', format='1E', unit='pixel',
                         disp='F11.4', array=cat_df['Y_IMAGE'])
        c4 = fits.Column(name='X2_IMAGE', format='1D', unit='pixel**2',
                         disp='E18.10', array=cat_df['X2_IMAGE'])
        c5 = fits.Column(name='Y2_IMAGE', format='1D', unit='pixel**2',
                         disp='E18.10', array=cat_df['Y2_IMAGE'])
        c6 = fits.Column(name='XY_IMAGE', format='1D', unit='pixel**2',
                         disp='E18.10', array=cat_df['XY_IMAGE'])
        c7 = fits.Column(name='ISOAREA_IMAGE', format='1J', unit='pixel**2',
                         disp='I9', array=cat_df['ISOAREA_IMAGE'])
        c8 = fits.Column(name='BACKGROUND', format='1E', unit='count',
                         disp='G12.7', array=cat_df['BACKGROUND'])
        c9 = fits.Column(name='THRESHOLD', format='1E', unit='count',
                         disp='G12.7', array=cat_df['THRESHOLD'])
        c10 = fits.Column(name='FLUX_MAX', format='1E', unit='count',
                          disp='G12.7', array=cat_df['FLUX_MAX'])
        c11 = fits.Column(name='A_IMAGE', format='1E', unit='pixel',
                          disp='F9.3', array=cat_df['A_IMAGE'])
        c12 = fits.Column(name='B_IMAGE', format='1E', unit='pixel',
                          disp='F9.3', array=cat_df['B_IMAGE'])
        c13 = fits.Column(name='THETA_IMAGE', format='1E', unit='deg',
                          disp='F6.2', array=cat_df['THETA_IMAGE'])
        c14 = fits.Column(name='ERRA_IMAGE', format='1E', unit='pixel',
                          disp='F9.5', array=cat_df['ERRA_IMAGE'])
        c15 = fits.Column(name='ERRB_IMAGE', format='1E', unit='pixel',
                          disp='F9.5', array=cat_df['ERRB_IMAGE'])
        c16 = fits.Column(name='FLUX_ISO', format='1E', unit='count',
                          disp='G12.7', array=cat_df['FLUX_ISO'])
        c17 = fits.Column(name='FLUXERR_ISO', format='1E', unit='count',
                          disp='G12.7', array=cat_df['FLUXERR_ISO'])
        c18 = fits.Column(name='MAG_ISO', format='1E', unit='mag',
                          disp='F8.4', array=cat_df['MAG_ISO'])
        c19 = fits.Column(name='MAGERR_ISO', format='1E', unit='mag',
                          disp='F8.4', array=cat_df['MAGERR_ISO'])
        c20 = fits.Column(name='FLUX_APER', format='1E', unit='count',
                          disp='G12.7', array=cat_df['FLUX_APER'])
        c21 = fits.Column(name='FLUXERR_APER', format='1E', unit='count',
                          disp='G12.7', array=cat_df['FLUXERR_APER'])
        c22 = fits.Column(name='MAG_APER', format='1E', unit='mag',
                          disp='F8.4', array=cat_df['MAG_APER'])
        c23 = fits.Column(name='MAGERR_APER', format='1E', unit='mag',
                          disp='F8.4', array=cat_df['MAGERR_APER'])
        c24 = fits.Column(name='ALPHA_SKY', format='1D', unit='deg',
                          disp='F11.7', array=cat_df['ALPHA_SKY'])
        c25 = fits.Column(name='DELTA_SKY', format='1D', unit='deg',
                          disp='F11.7', array=cat_df['DELTA_SKY'])
        c26 = fits.Column(name='ERRTHETA_IMAGE', format='1E', unit='deg',
                          disp='F6.2', array=cat_df['ERRTHETA_IMAGE'])
        c27 = fits.Column(name='MU_MAX', format='1E',
                          unit='mag * arcsec**(-2)', disp='F8.4',
                          array=cat_df['MU_MAX'])
        c28 = fits.Column(name='FWHM_IMAGE', format='1E', unit='pixel',
                          disp='F8.2', array=cat_df['FWHM_IMAGE'])
        c29 = fits.Column(name='FLUX_RADIUS', format='1E', unit='pixel',
                          disp='F10.3', array=cat_df['FLUX_RADIUS'])
        c30 = fits.Column(name='ELONGATION', format='1E', disp='F8.3',
                          array=cat_dfe['ELONGATION'])
        c31 = fits.Column(name='ELLIPTICITY', format='1E', disp='F8.3',
                          array=cat_df['ELLIPTICITY'])
        c32 = fits.Column(name='CXX_IMAGE', format='1E', unit='pixel**(-2)',
                          disp='E15.7', array=cat_df['CXX_IMAGE'])
        c33 = fits.Column(name='CXY_IMAGE', format='1E', unit='pixel**(-2)',
                          disp='E15.7', array=cat_df['CXY_IMAGE'])
        c34 = fits.Column(name='CYY_IMAGE', format='1E', unit='pixel**(-2)',
                          disp='E15.7', array=cat_df['CYY_IMAGE'])
        c35 = fits.Column(name='ERRCXX_IMAGE', format='1E',
                          unit='pixel**(-2)', disp='G12.7',
                          array=cat_df['ERRCXX_IMAGE'])
        c36 = fits.Column(name='ERRCXY_IMAGE', format='1E',
                          unit='pixel**(-2)', disp='G12.7',
                          array=cat_df['ERRCXY_IMAGE'])
        c37 = fits.Column(name='ERRCYY_IMAGE', format='1E', unit='pixel**(-2)',
                          disp='G12.7', array=cat_df['ERRCYY_IMAGE'])
        c38 = fits.Column(name='MAG_AUTO', format='1E', unit='mag',
                          disp='F8.4', array=cat_df['MAG_AUTO'])
        c39 = fits.Column(name='XWIN_IMAGE', format='1D', unit='pixel',
                          disp='F11.4', array=cat_df['XWIN_IMAGE'])
        c40 = fits.Column(name='YWIN_IMAGE', format='1D', unit='pixel',
                          disp='F11.4', array=cat_df['YWIN_IMAGE'])
        c41 = fits.Column(name='FLUX_AUTO', format='1E', unit='count',
                          disp='G12.7', array=cat_df['FLUX_AUTO'])
        c42 = fits.Column(name='FLUXERR_AUTO', format='1E', unit='count',
                          disp='G12.7', array=cat_df['FLUXERR_AUTO'])
        c43 = fits.Column(name='MAGERR_AUTO', format='1E', unit='mag',
                          disp='F8.4', array=cat_df['MAGERR_AUTO'])
        c44 = fits.Column(name='SNR_WIN', format='1E', disp='G10.4',
                          array=cat_df['SNR_WIN'])
        c45 = fits.Column(name='ALPHA_J2000', format='1D', unit='deg',
                          disp='F11.7', array=cat_df['ALPHA_J2000'])
        c46 = fits.Column(name='DELTA_J2000', format='1D', unit='deg',
                          disp='F11.7', array=cat_df['DELTA_J2000'])
        c47 = fits.Column(name='X_WORLD', format='1D', unit='deg',
                          disp='E18.10', array=cat_df['X_WORLD'])
        c48 = fits.Column(name='Y_WORLD', format='1D', unit='deg',
                          disp='E18.10', array=cat_df['Y_WORLD'])
        c49 = fits.Column(name='ERRX2_WORLD', format='1D', unit='deg**2',
                          disp='E18.10', array=cat_df['ERRX2_WORLD'])
        c50 = fits.Column(name='ERRY2_WORLD', format='1D', unit='deg**2',
                          disp='E18.10', array=cat_df['ERRY2_WORLD'])
        c51 = fits.Column(name='ERRXY_WORLD', format='1D', unit='deg**2',
                          disp='E18.10', array=cat_df['ERRXY_WORLD'])
        c52 = fits.Column(name='AWIN_IMAGE', format='1E', unit='pixel',
                          disp='F9.3', array=cat_df['AWIN_IMAGE'])
        c53 = fits.Column(name='BWIN_IMAGE', format='1E', unit='pixel',
                          disp='F9.3', array=cat_df['BWIN_IMAGE'])
        c54 = fits.Column(name='THETAWIN_IMAGE', format='1E', unit='deg',
                          disp='F6.2', array=cat_df['THETAWIN_IMAGE'])
        c55 = fits.Column(name='ERRAWIN_IMAGE', format='1E', unit='pixel',
                          disp='F9.5', array=cat_df['ERRAWIN_IMAGE'])
        c56 = fits.Column(name='ERRBWIN_IMAGE', format='1E', unit='pixel',
                          disp='F9.5', array=cat_df['ERRBWIN_IMAGE'])
        c57 = fits.Column(name='ERRTHETAWIN_IMAGE', format='1E', unit='deg',
                          disp='F6.2', array=cat_df['ERRTHETAWIN_IMAGE'])
        c58 = fits.Column(name='FLAGS', format='1I', disp='I3',
                          array=cat_df['FLAGS'])
        c59 = fits.Column(name='FWHM_WORLD', format='1E', unit='deg',
                          disp='G12.7', array=cat_df['FWHM_WORLD'])
        c60 = fits.Column(name='ERRA_WORLD', format='1E', unit='deg',
                          disp='G12.7', array=cat_df['ERRA_WORLD'])
        c61 = fits.Column(name='ERRB_WORLD', format='1E', unit='deg',
                          disp='G12.7', array=cat_df['ERRB_WORLD'])

        coldefs = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,
                                c12, c13, c14, c15, c16, c17, c18, c19, c20,
                                c21, c22, c23, c24, c25, c26, c27, c28, c29,
                                c30, c31, c32, c33, c34, c35, c36, c37, c38,
                                c39, c40, c41, c42, c43, c44, c45, c46, c47,
                                c48, c49, c50, c51, c52, c53, c54, c55, c56,
                                c57, c58, c59, c60, c61])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        cat_loc = '{}/{}/{}/catalog.cat'.format(prfs_d['catalogs_dir'],
                                                folder_sex, folder_scmp)

        """
        catalogue = fits.open(cat_location)
        catalogue[2] = tbhdu
        catalogue[2].header['EXTNAME'] = 'LDAC_OBJECTS'

        cat_loc = '{}/{}/{}/catalog.cat'.format(prfs_d['catalogs_dir'],
                                                folder_sex, folder_scmp)
        catalogue.writeto(cat_new_loc, clobber=True)
        """
        return True

    def merge_stats(self, fits_files, prfs_d):
        """

        @param fits_files:
        @param prfs_d:
        """
        cat_list = []

        for fits_ in fits_files:
            fits_file = '{}/{}.csv'.format(prfs_d['tmp_out'], fits_[-13:-5])
            cat_ = read_csv(fits_file, index_col=0)
            cat_list.append(cat_)

        for fits_ in fits_files:
            fits_file = '{}/{}.csv'.format(prfs_d['tmp_out'], fits_[-13:-5])
            remove(fits_file)

        stats_cat = concat(cat_list, axis=0)
        stats_cat.to_csv('stats.csv')


if __name__ == '__main__':
    test = Compare()
