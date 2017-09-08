#!/usr/bin/python
# -*- coding: utf-8 -*-

from multiprocessing import Manager, Process
from os import listdir

from astropy.io import fits
from astropy.wcs import WCS
from numpy import float64, genfromtxt, isclose, logical_and, random
from pandas import concat, read_csv, Series

from misc import extract_settings


class Check_PM:

    def __init__(self, logger):

        # Configuration parameters
        self.save = True

        prfs_d = extract_settings()

        out_dict, cat_files = self.check_position(logger, prfs_d)
        self.check_pm(logger, prfs_d, out_dict, cat_files)

    def check_pm(self, logger, prfs_d, out_dict, cat_files):
        """

        @param logger:
        @param prfs_d:
        """
        filters = []
        for cat_ in cat_files:
            cat_ = cat_.replace(prfs_d['results_dir'] + '/', '')
            filters.append(cat_)

        # append series to dataframe, dither
        # create three dataframes
        df_dict = {}

        for filter_ in filters:
            df_dict[filter_] = []

        for key_ in out_dict.keys():
            for filter_ in filters:
                if filter_ in key_:
                    df_dict[filter_].append(out_dict[key_])

        for key_ in df_dict.keys():
            df_dict[key_] = concat(df_dict[key_])
            df_dict[key_].to_csv('{}.csv'.format(key_))

        for key_ in df_dict.keys():
            m = df_dict[key_]['PM'] / df_dict[key_]['PMIN']
            df_dict[key_]['DIFF'] = m
            df_dict[key_].to_csv('test_{}.csv'.format(key_))

    def get_input_cat(self, logger, prfs_d, cat_file):
        """

        @param logger:
        @param prfs_d:
        @param cat_file:

        @return True: if everything goes alright.
        """

        catalogue = genfromtxt(cat_file)

        list_x = catalogue[:, 0]
        list_y = catalogue[:, 1]
        list_mag = catalogue[:, 2]
        list_pm = catalogue[:, 3]

        speed_0_001 = range(prfs_d['first_sso'], 137447, 75)
        speed_0_003 = range(prfs_d['first_sso'] + 10, 137457, 75)
        speed_0_01 = range(prfs_d['first_sso'] + 20, 137467, 75)
        speed_0_03 = range(prfs_d['first_sso'] + 30, 137477, 75)
        speed_0_1 = range(prfs_d['first_sso'] + 40, 137487, 75)
        speed_0_3 = range(prfs_d['first_sso'] + 50, 137497, 75)
        speed_1 = range(prfs_d['first_sso'] + 60, 137507, 75)

        speed_3 = range(prfs_d['first_sso'] + 66, 137512, 75)
        speed_10 = range(prfs_d['first_sso'] + 67, 137513, 75)
        speed_30 = range(prfs_d['first_sso'] + 68, 137514, 75)
        speed_100 = range(prfs_d['first_sso'] + 69, 137515, 75)
        speed_300 = range(prfs_d['first_sso'] + 70, 137516, 75)

        for index in speed_0_001:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 0.001
        for index in speed_0_003:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 0.003
        for index in speed_0_01:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 0.01
        for index in speed_0_03:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 0.03
        for index in speed_0_1:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 0.1
        for index in speed_0_3:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 0.3
        for index in speed_1:
            list_mag[index] = list_mag[index] - 2.5
            list_pm[index] = 1
        for index in speed_3:
            list_pm[index] = list_pm[index] - 1000
        for index in speed_10:
            list_pm[index] = list_pm[index] - 1000
        for index in speed_30:
            list_pm[index] = list_pm[index] - 1000
        for index in speed_100:
            list_pm[index] = list_pm[index] - 1000
        for index in speed_300:
            list_pm[index] = list_pm[index] - 1000

        indexes = (speed_0_001 + speed_0_003 + speed_0_01 + speed_0_03 +
                   speed_0_1 + speed_0_3 + speed_1 + speed_3 + speed_10 +
                   speed_30 + speed_100 + speed_300)
        indexes = sorted(indexes)

        s1 = Series(list_x, name='X_IMAGE', dtype=float64)
        s2 = Series(list_y, name='Y_IMAGE', dtype=float64)
        s3 = Series(list_mag, name='mag_values', dtype=float64)
        s4 = Series(list_pm, name='PM_INPUT', dtype=float64)

        sources_df = concat([s1, s2, s3, s4], axis=1)
        sources_df = sources_df.iloc[indexes, :]

        return sources_df

    def check_position(self, logger, prfs_d):
        """

        @param logger:
        @param prfs_d:
        """

        cat_files = []
        cat_files_d = listdir(prfs_d['results_dir'])
        for cat_file in cat_files_d:
            if cat_file[-5:] == '6.csv':
                cat_files.append(prfs_d['results_dir'] + '/' + cat_file)

        i_cats = {}
        i_cat_files = listdir(prfs_d['input_cats'])
        # Read  input catalog files
        for cat in i_cat_files:
            # All dithers, catalog should start with 'mag' string
            if cat[-4:] == '.dat':
                for mag in prfs_d['mags']:
                    if cat[-12:-7] == mag:
                        cat_name = prfs_d['input_cats'] + '/' + cat
                        catalog = self.get_input_cat(logger, prfs_d, cat_name)
                        i_cats['{}'.format(int(cat[-5:-4]))] = catalog

        fits_loc = '{}/m_{}_x0_y0_d1_copy.fits'.format(prfs_d['fits_dir'], mag)
        hdulist = fits.open(fits_loc)
        w = WCS(hdulist[0].header)

        i_cats_t = {}

        for dither in range(1, 5, 1):
            input_table = i_cats[str(dither)]
            regions_list = []
            for source_num in range(input_table['X_IMAGE'].as_matrix().size):
                x_value = input_table['X_IMAGE'].as_matrix()[source_num]
                y_value = input_table['Y_IMAGE'].as_matrix()[source_num]
                regions_list.append([x_value, y_value])

            logger.debug('transforming x/y values to ra/dec ones')
            input_regions = w.wcs_pix2world(regions_list, 1)

            logger.debug('creating new list with ra/dec values')
            alpha_list = []
            delta_list = []
            for regions in input_regions:
                alpha_list.append(regions[0])
                delta_list.append(regions[1])

            mag_values_list = input_table['mag_values'].tolist()
            pm_values_list = input_table['PM_INPUT'].tolist()

            mag_values = Series(mag_values_list, name='mag_values')
            pm_values = Series(pm_values_list, name='PM_INPUT')

            logger.debug('creating new Pandas series for ra/dec values')
            alpha_j2000 = Series(alpha_list, name='alpha_j2000')
            delta_j2000 = Series(delta_list, name='delta_j2000')

            logger.debug('concatenating series of ra and dec')
            coords_table = concat([alpha_j2000, delta_j2000,
                                   mag_values, pm_values], axis=1)

            """
            coords_table.to_csv('dither_{}.csv'.format(dither), index=False,
                                header=False, sep=" ")
            """
            i_cats_t[dither] = coords_table

        cat_n = {1: [1, 2, 3, 4, 5, 6, 7, 8, 9],
                 2: [10, 11, 12, 13, 14, 15, 16, 17, 18],
                 3: [19, 20, 21, 22, 23, 24, 25, 26, 27],
                 4: [28, 29, 30, 31, 32, 33, 34, 35, 36]}

        o_cats = {}
        for cat_file in cat_files:
            # cat_file[70:] id!
            filter_n = cat_file[70:]  # harcoded!
            cat = read_csv(cat_file, index_col=0)
            for dither in range(1, 5, 1):
                cat_d = cat[cat['CATALOG_NUMBER'].isin(cat_n[dither])]

                o_cats['{}_{}'.format(filter_n, dither)] = cat_d

        manager = Manager()
        out_dict = manager.dict()

        check_j = []
        for proc in range(1, 5, 1):
            # Each dither in a different thread

            temp_keys = []

            for key in o_cats.keys():
                if key[-1:] == str(proc):
                    temp_keys.append(key)

            o_cat_t = {new_key: o_cats[new_key] for new_key in temp_keys}

            check_p = Process(target=self.check_thread,
                              args=(logger, prfs_d, i_cats_t[proc],
                                    o_cat_t, proc, out_dict))
            check_j.append(check_p)
            check_p.start()

        active_check = list([job.is_alive() for job in check_j])
        while True in active_check:
            active_check = list([job.is_alive() for job in check_j])
            pass

        return out_dict, cat_files

    def check_thread(self, logger, prfs_d, i_cats_t, o_cats, dither, out_dict):
        """

        @param logger:
        @param prfs_d:
        @param i_cats_t:
        @param o_cats:

        @return True:
        """

        # counts
        stats_dict = {}
        for d in range(1, 5, 1):
            for o_cat in o_cats.keys():
                stats_dict['{}_{}'.format(dither, o_cat)] = {}

        columns = ['SOURCE_NUMBER', 'CATALOG_NUMBER', 'X_IMAGE', 'Y_IMAGE',
                   'ERRA_IMAGE', 'ERRB_IMAGE', 'ERRTHETA_IMAGE', 'ALPHA_J2000',
                   'DELTA_J2000', 'ERRA_WORLD', 'ERRB_WORLD', 'ERRTHETA_WORLD',
                   'EPOCH', 'PM', 'PMERR', 'PMALPHA', 'PMDELTA',
                   'PMALPHAERR', 'PMDELTAERR']

        for key_ in stats_dict.keys():
            for column in columns:
                stats_dict[key_][column] = []
            stats_dict[key_]['PMIN'] = []
            stats_dict[key_]['DITHER'] = []

        sources_number = i_cats_t['alpha_j2000'].size
        for source_num in range(sources_number):
            # FIXME no uso valores absolutos!
            source_alpha = i_cats_t['alpha_j2000'].iloc[source_num]
            source_delta = i_cats_t['delta_j2000'].iloc[source_num]

            for o_cat in o_cats.keys():
                alpha_m = isclose(o_cats[o_cat]['ALPHA_J2000'].values[:, None],
                                  [source_alpha],
                                  atol=prfs_d['alpha_t']).any(axis=1)
                delta_m = isclose(o_cats[o_cat]['DELTA_J2000'].values[:, None],
                                  [source_delta],
                                  atol=prfs_d['delta_t']).any(axis=1)

                m = logical_and(alpha_m, delta_m)
                df = o_cats[o_cat][m]
                if df.size == 27:
                    dict_name = '{}_{}'.format(dither, o_cat)
                    for column in columns:
                        stats_dict[dict_name][column].append(df[column].iloc[0])
                    pm_input = i_cats_t['PM_INPUT'].iloc[source_num]
                    stats_dict[dict_name]['PMIN'].append(pm_input)
                    stats_dict[dict_name]['DITHER'].append(dither)

        # Creates a collection of Pandas Series from dicts
        series_dict = {}
        for key_ in stats_dict.keys():
            series_dict[key_] = {}
            for column in columns:
                series_dict[key_][column] = Series(stats_dict[key_][column],
                                                   name=column)
            series_dict[key_]['PMIN'] = Series(stats_dict[key_]['PMIN'],
                                               name='PMIN')
            series_dict[key_]['DITHER'] = Series(stats_dict[key_]['DITHER'],
                                                 name='DITHER')

        for key_ in series_dict.keys():
            tmp_list = []
            for column in columns:
                tmp_list.append(series_dict[key_][column])
            tmp_list.append(series_dict[key_]['PMIN'])
            tmp_list.append(series_dict[key_]['DITHER'])
            out_dict[key_] = concat(tmp_list, axis=1)
