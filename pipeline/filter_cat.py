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
from numpy import array
from pandas import DataFrame, read_csv

from cats_management import check_source, get_input_dicts
from misc import extract_settings


def create_catalog():
    pass
    # To get out!
    # idx_cat = o_cat[o_cat['NUMBER'].isin([1])]
    # idx_list = idx_cat['IDX'].tolist()
    #
    # dither_n = 1
    # for i, idx_ in enumerate(idx_list):
    #     if i != (len(idx_list) - 1):
    #         lower_value = idx_
    #         upper_value = idx_list[i + 1]
    #         for dither_tmp in range(0, upper_value - lower_value, 1):
    #             dither_list.append(dither_n)
    #     else:
    #         lower_value = idx_
    #         upper_value = len(o_cat['NUMBER'].tolist())
    #         for dither_tmp in range(0, upper_value - lower_value, 1):
    #             dither_list.append(dither_n)
    #     if dither_n == 4:
    #         dither_n = 0
    #     dither_n += 1
    #
    # print('len_1 test {}'.format(len(o_cat['NUMBER'].tolist())))
    # print('len_2 test {}'.format(len(dither_list)))
    #
    # print(o_cat.columns)
    # print(o_cat['NUMBER'].size)
    #
    # o_cat['DITHER'] = dither_list
    # o_cat.to_csv('/media/sf_CarpetaCompartida/catalog_n.csv')


def function_returner(p, input_data):
    input_data = array(input_data)
    x = input_data[0]
    z = input_data[1]

    return p[0] * x + p[1] * z + p[2]


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


class SlipFullSextractorCatalog:

    def __init__(self):
        """

        """
        self.prfs_d = extract_settings()

        self.data_d = redo_dict_keys()  # Creates a dictionary for statistics
        self.data_d['object'] = []
        self.data_d['input_pm'] = []

        self.get()  # Gets data from catalogs
        self.save_df()

    def save_df(self):
        """

        :return: True
        """
        dir_ = '/media/sf_CarpetaCompartida/'
        for object_ in ['star', 'galaxy', 'sso', 'empty']:
            object_df = DataFrame(self.data_d)
            object_df = object_df[object_df['object'].isin([object_])]
            object_df.to_csv('{}{}.csv'.format(dir_, object_))

        return True

    def get(self):
        """

        :return:
        """
        mag = '20-21'
        # Gets dataframes for each type of source
        i_df = get_input_dicts(mag)

        # Gets the name of filtered file
        filter_o_n = '/media/sf_CarpetaCompartida/catalog_n.csv'
        # Opens filtered file
        o_cat = read_csv('{}'.format(filter_o_n), index_col=0)

        keys_l = o_cat.columns.tolist()

        idx_test = 0
        # Loops over unique sources of filtered file
        for idx_, row in enumerate(o_cat.itertuples(), 1):
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
    SlipFullSextractorCatalog()
