#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for performance
Mide las detecciones de cada tipo

Versions:
- 0.1

Todo:
    * Improve log messages
    * Use same naming convention for all variables

"""
from astropy.io import fits
from astropy.table import Table
from pandas import concat, DataFrame, read_csv, Series

from cats_management import check_source, get_input_dicts
from misc import extract_settings, get_dither, get_cats


def redo_dict_keys():
    """ Creates a dictionary

    :return: dict_
    """
    dict_ = {'X_IMAGE': [], 'Y_IMAGE': [], 'X2_IMAGE': [], 'Y2_IMAGE': [],
             'XY_IMAGE': [], 'ISOAREA_IMAGE': [], 'BACKGROUND': [],
             'THRESHOLD': [], 'FLUX_MAX': [], 'A_IMAGE': [], 'B_IMAGE': [],
             'THETA_IMAGE': [], 'ERRA_IMAGE': [], 'ERRB_IMAGE': [],
             'FLUX_ISO': [], 'FLUXERR_ISO': [], 'MAG_ISO': [],
             'MAGERR_ISO': [], 'FLUX_APER': [], 'FLUXERR_APER': [],
             'MAG_APER': [], 'MAGERR_APER': [], 'ALPHA_SKY': [],
             'DELTA_SKY': [], 'ERRTHETA_IMAGE': [], 'MU_MAX': [],
             'FWHM_IMAGE': [], 'CLASS_STAR': [], 'FLUX_RADIUS': [],
             'ELONGATION': [], 'ELLIPTICITY': [], 'CXX_IMAGE': [],
             'CXY_IMAGE': [], 'CYY_IMAGE': [], 'ERRCXX_IMAGE': [],
             'ERRCXY_IMAGE': [], 'ERRCYY_IMAGE': [], 'MAG_AUTO': [],
             'XWIN_IMAGE': [], 'YWIN_IMAGE': [], 'FLUX_AUTO': [],
             'FLUXERR_AUTO': [], 'MAGERR_AUTO': [], 'SNR_WIN': [],
             'ALPHA_J2000': [], 'DELTA_J2000': [], 'X_WORLD': [],
             'Y_WORLD': [], 'ERRX2_WORLD': [], 'ERRY2_WORLD': [],
             'ERRXY_WORLD': [], 'AWIN_IMAGE': [], 'BWIN_IMAGE': [],
             'FLAGS': [], 'FWHM_WORLD': [], 'A_WORLD': [], 'B_WORLD': [],
             'THETA_WORLD': [], 'ERRA_WORLD': [], 'ERRB_WORLD': [],
             'THETAWIN_IMAGE': [], 'ERRAWIN_IMAGE': [], 'ERRBWIN_IMAGE': [],
             'ERRTHETAWIN_IMAGE': []}

    return dict_


def redo_dict_keys_scamp():
    """ Creates a dictionary

    :return: dict_
    """
    dict_ = {'SOURCE_NUMBER': [], 'CATALOG_NUMBER': [], 'X_IMAGE': [],
             'Y_IMAGE': [], 'ERRA_IMAGE': [], 'ERRB_IMAGE': [],
             'ERRTHETA_IMAGE': [], 'ALPHA_J2000': [], 'DELTA_J2000': [],
             'ERRA_WORLD': [], 'ERRB_WORLD': [], 'ERRTHETA_WORLD': [],
             'EPOCH': [], 'MAG': [], 'MAGERR': [], 'FLAGS_EXTRACTION': [],
             'FLAGS_SCAMP': [], 'FLAGS_IMA': [], 'PM': [], 'PMERR': [],
             'PMALPHA': [], 'PMDELTA': [], 'PMALPHAERR': [], 'PMDELTAERR': [],
             'ISOAREA_IMAGE': [], 'A_IMAGE': [], 'B_IMAGE': [],
             'THETA_IMAGE': [], 'FWHM_IMAGE': [], 'FLUX_RADIUS': [],
             'ELONGATION': [], 'ELLIPTICITY': [], 'MAG_ISO': [],
             'CLASS_STAR': []}

    return dict_


class FilterByPM:

    def __init__(self, pm_):
        """

        """
        # Variables definition
        self.dir_ = '/media/sf_CarpetaCompartida'  # change if necessary
        self.pm = pm_

        # Process
        self.sso_d = self.open_csv()
        self.filter_d = self.filter_cat()
        self.save_df()

    def open_csv(self):
        """

        :return: sso_d
        """
        csv_name = '{}/sso.csv'.format(self.dir_)
        sso_d = read_csv(csv_name, index_col=0)

        return sso_d

    def filter_cat(self):
        """

        :return: filter_d
        """
        filter_d = self.sso_d
        filter_d = filter_d[filter_d['input_pm'] < 0.3]

        return filter_d

    def save_df(self):
        """

        :return:
        """
        self.filter_d.to_csv('{}/sso_{}.csv'.format(self.dir_, self.pm))

        return True


class SlipFullSextractorCatalog:

    def __init__(self):
        """

        """
        self.dir_ = '/media/sf_CarpetaCompartida'
        self.prfs_d = extract_settings()

        self.data_d = redo_dict_keys()  # Creates a dictionary for statistics
        self.data_d['object'] = []
        self.data_d['input_pm'] = []

        self.full_catalog = self.merge_ccd_catalogs()
        self.get()  # Gets data from catalogs
        self.save_df()

    def merge_ccd_catalogs(self):
        """

        :return:
        """
        d = '{}/luca_data/v18/20-21/CCDs/30_1.5_1.5_0.01_4/'.format(self.dir_)
        data_list = []
        cats = get_cats()
        for key_ in cats.keys():
            for cat_ in cats[key_]:
                # Creates CCD and dither column
                cat_location = '{}{}'.format(d, cat_)
                hdu_list = fits.open(cat_location)
                data = Table(hdu_list[2].data).to_pandas()

                # Number of rows
                l_length = int(data['NUMBER'].size)
                ccd = cat_[-12:-7]
                ccd_l = [ccd] * l_length
                dither = cat_[-5:-4]
                dither_l = [dither] * l_length

                # Creates Series
                ccd_s = Series(ccd_l, name='CCD')
                data['CCD'] = ccd_s
                dither_s = Series(dither_l, name='DITHER')
                data['DITHER'] = dither_s

                # Appends to list
                data_list.append(data)

        data_df = concat(data_list, axis=0)
        data_df = data_df.reset_index(drop=True)

        return data_df

    def check_source(self, catalog_n, i_df, o_alpha, o_delta):
        """

        :param catalog_n:
        :param i_df:0
        :param o_alpha:
        :param o_delta:
        :return:
        """

        if catalog_n != 0:
            i_df = i_df[i_df['catalog'].isin([catalog_n])]

        i_df = i_df[i_df['alpha_j2000'] + self.prfs_d['tolerance'] > o_alpha]
        i_df = i_df[o_alpha > i_df['alpha_j2000'] - self.prfs_d['tolerance']]
        i_df = i_df[i_df['delta_j2000'] + self.prfs_d['tolerance'] > o_delta]
        i_df = i_df[o_delta > i_df['delta_j2000'] - self.prfs_d['tolerance']]

        return i_df

    def save_df(self):
        """

        :return: True
        """
        for object_ in ['star', 'galaxy', 'sso', 'empty']:
            object_df = DataFrame(self.data_d)
            object_df = object_df[object_df['object'].isin([object_])]
            object_df.to_csv('{}/full_{}.csv'.format(self.dir_, object_))

        return True

    def get(self):
        """

        :return:
        """
        mag = '20-21'
        # Gets dataframes for each type of source
        i_df = get_input_dicts(mag)

        keys_l = self.full_catalog.columns.tolist()

        idx_test = 0
        # Loops over unique sources of filtered file
        for idx_, row in enumerate(self.full_catalog.itertuples(), 1):
            tmp_d = {'object': 'empty'}
            catalog_n = 0

            print('idx value {}'.format(idx_test))
            idx_test += 1

            for idx_k, key_ in enumerate(keys_l):
                tmp_d[key_] = row[idx_k + 1]

            # look for object type!
            dither = row.DITHER

            i_stars_df = i_df['stars']
            i_stars_d_df = i_stars_df[
                i_stars_df['dither_values'].isin([dither])]

            i_galaxies_df = i_df['galaxies']
            i_galaxies_d_df = i_galaxies_df[
                i_galaxies_df['dither_values'].isin([dither])]

            i_ssos_df = i_df['ssos']
            i_ssos_d_df = i_ssos_df[
                i_ssos_df['dither_values'].isin([dither])]

            out_df = self.check_source(catalog_n, i_stars_d_df,
                                       tmp_d['ALPHA_J2000'],
                                       tmp_d['DELTA_J2000'])
            if out_df.empty is not True:
                # Gets proper motion
                tmp_d['object'] = 'star'

            out_df = self.check_source(catalog_n, i_galaxies_d_df,
                                       tmp_d['ALPHA_J2000'],
                                       tmp_d['DELTA_J2000'])
            if out_df.empty is not True:
                # Gets proper motion
                tmp_d['object'] = 'galaxy'

            out_df = self.check_source(catalog_n, i_ssos_d_df,
                                       tmp_d['ALPHA_J2000'],
                                       tmp_d['DELTA_J2000'])
            if out_df.empty is not True:
                # Gets proper motion
                tmp_d['object'] = 'sso'
                tmp_d['input_pm'] = out_df['pm_values'].iloc[0]
                print(out_df['pm_values'].iloc[0])

            # Appends data to temporal dictionary
            for key_ in self.data_d.keys():
                if key_ != 'input_pm':
                    self.data_d[key_].append(tmp_d[key_])
                else:
                    if tmp_d['object'] == 'sso':
                        self.data_d['input_pm'].append(tmp_d['input_pm'])
                    else:
                        self.data_d['input_pm'].append(0)


class SlipFilteredScampCatalog:

    def __init__(self):
        """

        """
        self.dir_ = '/media/sf_CarpetaCompartida/'
        self.prfs_d = extract_settings()

        self.data_d = redo_dict_keys_scamp()  # Creates a dict for statistics
        self.data_d['object'] = []
        self.data_d['input_pm'] = []

        self.get()  # Gets data from catalogs
        self.save_df()

    def save_df(self):
        """

        :return: True
        """
        for object_ in ['star', 'galaxy', 'sso', 'empty']:
            object_df = DataFrame(self.data_d)
            object_df = object_df[object_df['object'].isin([object_])]
            object_df.to_csv('{}{}_filt_6.csv'.format(self.dir_, object_))

        return True

    def get(self):
        """

        :return:
        """
        mag = '20-21'
        # Gets dataframes for each type of source
        i_df = get_input_dicts(mag)

        # Gets the name of filtered file
        filter_o_n = 'filt_10_1.1_0.5_0.04_20-21_6.csv'
        # Opens filtered file
        o_cat = read_csv('{}{}'.format(self.dir_, filter_o_n), index_col=0)

        keys_l = o_cat.columns.tolist()

        print('SOURCE {}'.format(len(o_cat['SOURCE_NUMBER'].tolist())))
        idx_test = 0
        # Loops over unique sources of filtered file
        for idx_, row in enumerate(o_cat.itertuples(), 1):
            tmp_d = {'object': 'empty'}
            catalog_n = 0

            print('idx value {}'.format(idx_test))
            idx_test += 1

            for idx_k, key_ in enumerate(keys_l):
                tmp_d[key_] = row[idx_k + 1]

            ccd, dither = get_dither(row.CATALOG_NUMBER)

            i_stars_df = i_df['stars']
            i_stars_d_df = i_stars_df[
                i_stars_df['dither_values'].isin([dither])]

            i_galaxies_df = i_df['galaxies']
            i_galaxies_d_df = i_galaxies_df[
                i_galaxies_df['dither_values'].isin([dither])]

            i_ssos_df = i_df['ssos']
            i_ssos_d_df = i_ssos_df[
                i_ssos_df['dither_values'].isin([dither])]

            out_df = check_source(catalog_n, i_stars_d_df,
                                  tmp_d['ALPHA_J2000'], tmp_d['DELTA_J2000'])
            if out_df.empty is not True:
                # Gets proper motion
                tmp_d['object'] = 'star'

            out_df = check_source(catalog_n, i_galaxies_d_df,
                                  tmp_d['ALPHA_J2000'], tmp_d['DELTA_J2000'])
            if out_df.empty is not True:
                # Gets proper motion
                tmp_d['object'] = 'galaxy'

            out_df = check_source(catalog_n, i_ssos_d_df,
                                  tmp_d['ALPHA_J2000'], tmp_d['DELTA_J2000'])
            if out_df.empty is not True:
                # Gets proper motion
                tmp_d['object'] = 'sso'
                tmp_d['input_pm'] = out_df['pm_values'].iloc[0]
                print(out_df['pm_values'].iloc[0])

            # Appends data to temporal dictionary
            for key_ in self.data_d.keys():
                if key_ != 'input_pm':
                    self.data_d[key_].append(tmp_d[key_])
                else:
                    if tmp_d['object'] == 'sso':
                        self.data_d['input_pm'].append(tmp_d['input_pm'])
                    else:
                        self.data_d['input_pm'].append(0)


if __name__ == "__main__":
    # Comment as necessary
    SlipFilteredScampCatalog()
    # SlipFullSextractorCatalog()
    # pm = argv[1]
    # FilterByPM(pm)
