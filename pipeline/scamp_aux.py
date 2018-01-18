#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.
- 0.2: Functions reorganised to two different classes.
- 0.3: New organization for filter Class
- 0.4: Filter now gets areas

"""

from subprocess import Popen

from astropy.io import fits
from astropy.table import Table
from pandas import concat, read_csv, Series

from misc import pm_compute, pm_filter, extract_settings
from misc import motion_filter, create_folder, sn_filter

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.4"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class Scamp:

    def __init__(self, logger, mag, scmp_d, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_d:
        :param scmp_cf:
        :param sex_d:
        """
        self.prfs_d = extract_settings()

        self.scamp_process(logger, mag, scmp_d, scmp_cf, sex_d)

    def scamp_process(self, logger, mag, scmp_d, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_d:
        :param scmp_cf:
        :param sex_d:
        :return:
        """
        sex_cf = '{}_{}_{}_{}_{}'.format(sex_d['deblend_nthresh'],
                                         sex_d['analysis_thresh'],
                                         sex_d['detect_thresh'],
                                         sex_d['deblend_mincount'],
                                         sex_d['detect_minarea'])

        logger.info('scamp process for magnitude {}'.format(mag))
        logger.info('sextractor configuration {}'.format(sex_cf))
        logger.info('scamp configuration {}'.format(scmp_cf))

        scmp_p_1 = 'scamp -c {}'.format(self.prfs_d['conf_scamp'])
        sex_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], mag,
                                         sex_cf)
        sex_output = 'mag_{}_CCD_x?_y?_d?.cat'.format(mag)
        scmp_p_2 = ' {}/{}'.format(sex_loc, sex_output)
        scmp_p_3 = ' -ASTREFCAT_NAME'
        cat_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], mag,
                                         sex_cf)
        cat_input = 'catalog_{}.cat'.format(mag)
        scmp_p_4 = ' {}/{}'.format(cat_loc, cat_input)
        scmp_p_5 = ' -PIXSCALE_MAXERR {}'.format(scmp_d['pixscale_maxerr'])
        scmp_p_6 = ' -POSANGLE_MAXERR {}'.format(scmp_d['posangle_maxerr'])
        scmp_p_7 = ' -POSITION_MAXERR {}'.format(scmp_d['position_maxerr'])
        scmp_p_8 = ' -CROSSID_RADIUS {}'.format(scmp_d['crossid_radius'])
        # Output catalogs location
        cats_dir = '{}/{}/{}/{}'.format(self.prfs_d['catalogs_dir'], mag,
                                        sex_cf, scmp_cf)
        merged_cat = '{}/merged_{}_{}.cat'.format(cats_dir, scmp_cf, mag)
        scmp_p_9 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat)
        full_cat = '{}/full_{}_{}.cat'.format(cats_dir, scmp_cf, mag)
        scmp_p_10 = ' -FULLOUTCAT_NAME {}'.format(full_cat)
        scmp_p = scmp_p_1 + scmp_p_2 + scmp_p_3 + scmp_p_4 + scmp_p_5
        scmp_p = scmp_p + scmp_p_6 + scmp_p_7 + scmp_p_8 + scmp_p_9
        scmp_p = scmp_p + scmp_p_10

        # print('scmp_p {}'.format(scmp_p))

        create_folder(logger, cats_dir)

        process_scamp = Popen(scmp_p, shell=True)
        process_scamp.wait()

        logger.info('Scamp process finished. Data is ready to be analysed.')

        return True


class ScampFilter:  # TODO Split scamp_filter method into single methodss

    def __init__(self, logger, mag, scmp_cf, sex_d):
        """

        :param logger:
        :param mag:
        :param scmp_cf:
        :param sex_d:
        """
        self.prfs_d = extract_settings()
        self.logger = logger
        self.mag = mag
        self.scmp_cf = scmp_cf
        self.sex_cf = '{}_{}_{}_{}_{}'.format(sex_d['deblend_nthresh'],
                                              sex_d['analysis_thresh'],
                                              sex_d['detect_thresh'],
                                              sex_d['deblend_mincount'],
                                              sex_d['detect_minarea'])

        self.save = True
        (merged_db, full_db, filter_o_n) = self.scamp_filter()
        full_db = self.compute_pm(merged_db, full_db, filter_o_n)
        full_db = self.get_areas(full_db, filter_o_n)
        # self.filter_pm(full_db, filter_o_n)

    def check_source(self, o_df, o_alpha, o_delta):
        """

        :param o_df:
        :param o_alpha:
        :param o_delta:
        :return:
        """
        tolerance = 0.0001389  # 0.5 arcsecond

        o_df = o_df[o_df['ALPHA_J2000'] + tolerance > o_alpha]
        o_df = o_df[o_alpha > o_df['ALPHA_J2000'] - tolerance]
        o_df = o_df[o_df['DELTA_J2000'] + tolerance > o_delta]
        o_df = o_df[o_delta > o_df['DELTA_J2000'] - tolerance]

        return o_df

    def get_cat(self, catalog_n):
        """

        :param catalog_n:
        :return: cat_file
        """
        cats = [['x0_y0', 1, 1], ['x0_y0', 2, 2], ['x0_y0', 3, 3],
                ['x0_y0', 4, 4], ['x0_y1', 1, 5], ['x0_y1', 2, 6],
                ['x0_y1', 3, 7], ['x0_y1', 4, 8], ['x0_y2', 1, 9],
                ['x0_y2', 2, 10], ['x0_y2', 3, 11], ['x0_y2', 4, 12],
                ['x1_y0', 1, 13], ['x1_y0', 2, 14], ['x1_y0', 3, 15],
                ['x1_y0', 4, 16], ['x1_y1', 1, 17], ['x1_y1', 2, 18],
                ['x1_y1', 3, 19], ['x1_y1', 4, 20], ['x1_y2', 1, 21],
                ['x1_y2', 2, 22], ['x1_y2', 3, 23], ['x1_y2', 4, 24],
                ['x2_y0', 1, 25], ['x2_y0', 2, 26], ['x2_y0', 3, 27],
                ['x2_y0', 4, 28], ['x2_y1', 1, 29], ['x2_y1', 2, 30],
                ['x2_y1', 3, 31], ['x2_y1', 4, 32], ['x2_y2', 1, 33],
                ['x2_y2', 2, 34], ['x2_y2', 3, 35], ['x2_y2', 4, 36]]

        ccd = ''
        dither = ''
        for cat_ in cats:
            if cat_[2] == catalog_n:
                ccd = cat_[0]
                dither = cat_[1]

        cat_loc = '{}/{}/CCDs/{}/'.format(self.prfs_d['fits_dir'], self.mag,
                                          self.sex_cf)
        cat_name = 'mag_{}_CCD_{}_d{}.cat'.format(self.mag, ccd, dither)
        cat_file = '{}{}'.format(cat_loc, cat_name)

        return cat_file

    def scamp_filter(self):
        """

        :return:
        """
        self.logger.info("Filtering scamp's output")

        filter_dir = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'], self.mag,
                                          self.sex_cf, self.scmp_cf)
        create_folder(self.logger, filter_dir)

        self.logger.debug('Scamp configuration {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration {}'.format(self.sex_cf))

        # Full catalog name
        full_n_dir = '{}/{}/{}/{}/'.format(self.prfs_d['catalogs_dir'],
                                           self.mag, self.sex_cf, self.scmp_cf)
        full_n_cat = 'full_{}_{}_1.cat'.format(self.scmp_cf, self.mag)
        full_n = '{}{}'.format(full_n_dir, full_n_cat)
        # Filtered catalog name
        filt_n = 'filt_{}_{}'.format(self.scmp_cf, self.mag)

        self.logger.debug('opening full catalog {}'.format(full_n))
        full_cat = fits.open(full_n)
        full_db = Table(full_cat[2].data)
        self.logger.debug('converting full catalog to Pandas format')
        full_db = full_db.to_pandas()

        # Getting merge catalog
        mrgd_n_dir = '{}/{}/{}/{}/'.format(self.prfs_d['catalogs_dir'],
                                           self.mag, self.sex_cf, self.scmp_cf)
        mrgd_n_cat = 'merged_{}_{}_1.cat'.format(self.scmp_cf, self.mag)
        mrgd_n = '{}{}'.format(mrgd_n_dir, mrgd_n_cat)

        self.logger.debug('opening merged catalog {}'.format(mrgd_n))
        merged_cat = fits.open(mrgd_n)
        self.logger.debug('converting merged catalog to Pandas format')
        merged_db = Table(merged_cat[2].data)

        filter_o_n = '{}/{}'.format(filter_dir, filt_n)

        # Removing 0 catalog detections
        self.logger.debug('removing 0 catalog detections')
        full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

        full_db = concat(g for _, g in full_db.groupby("SOURCE_NUMBER")
                         if len(g) >= int(self.prfs_d['detections']))

        if self.save:
            full_db.to_csv('{}_1.csv'.format(filter_o_n))

        return merged_db, full_db, filter_o_n

    def compute_pm(self, merged_db, full_db, filter_o_n):
        """

        :param merged_db:
        :param full_db:
        :param filter_o_n:
        :return:
        """
        # Computing pm
        self.logger.debug('computing proper motion')
        full_db = pm_compute(self.logger, merged_db, full_db)
        if self.save:
            self.logger.debug('saving output to {}_2.csv'.format(filter_o_n))
            full_db.to_csv('{}_2.csv'.format(filter_o_n))

        return full_db

    def get_areas(self, full_db, filter_o_n):
        """

        :param full_db:
        :param filter_o_n:
        :return:
        """
        self.logger.debug('getting areas')

        # Opens filtered file
        filter_cat = read_csv('{}_2.csv'.format(filter_o_n), index_col=0)
        # Gets unique sources from filtered file
        unique_sources = list(set(filter_cat['SOURCE_NUMBER'].tolist()))
        print(len(unique_sources))

        isoarea_l = []
        a_image_l = []
        b_image_l = []
        thetha_image_l = []
        erra_image_l = []
        errb_image_l = []
        fwhm_image_l = []
        flux_radius_l = []
        elongation_l = []
        ellipticity_l = []
        # Loops over unique sources of filtered file
        for idx, source_ in enumerate(unique_sources):
            print(idx)
            o_df = filter_cat[filter_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                catalog_n = row.CATALOG_NUMBER
                cat_file = self.get_cat(catalog_n)
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000

                # self.logger.debug('opening CCD catalog {}'.format(cat_file))
                ccd_cat = fits.open(cat_file)
                ccd_df = Table(ccd_cat[2].data)
                # self.logger.debug('converting CCD catalog to Pandas format')
                ccd_df = ccd_df.to_pandas()

                cat_df = self.check_source(ccd_df, o_alpha, o_delta)
                if cat_df.empty:
                    isoarea_l.append('nan')
                else:
                    isoarea = cat_df['ISOAREA_IMAGE'].iloc[0]
                    isoarea_l.append(isoarea)
                    a_image = cat_df['A_IMAGE'].iloc[0]
                    a_image_l.append(a_image)
                    b_image = cat_df['B_IMAGE'].iloc[0]
                    b_image_l.append(b_image)
                    thetha_image = cat_df['THETA_IMAGE'].iloc[0]
                    thetha_image_l.append(thetha_image)
                    erra_image = cat_df['ERRA_IMAGE'].iloc[0]
                    erra_image_l.append(erra_image)
                    errb_image = cat_df['ERRB_IMAGE'].iloc[0]
                    errb_image_l.append(errb_image)
                    fwhm_image = cat_df['FWHM_IMAGE'].iloc[0]
                    fwhm_image_l.append(fwhm_image)
                    flux_radius = cat_df['FLUX_RADIUS'].iloc[0]
                    flux_radius_l.append(flux_radius)
                    elongation = cat_df['ELONGATION'].iloc[0]
                    elongation_l.append(elongation)
                    ellipticity = cat_df['ELLIPTICITY'].iloc[0]
                    ellipticity_l.append(ellipticity)

        # print(isoarea_l)
        isoarea_s = Series(isoarea_l)
        a_image_s = Series(a_image_l)
        b_image_s = Series(b_image_l)
        theta_image_s = Series(thetha_image_l)
        erra_image_s = Series(erra_image_l)
        errb_image_s = Series(errb_image_l)
        fwhm_image_s = Series(fwhm_image_l)
        flux_radius_s = Series(flux_radius_l)
        elongation_s = Series(elongation_l)
        ellipticity_s = Series(ellipticity_l)

        for i in set(full_db['SOURCE_NUMBER']):
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'ISOAREA_IMAGE'] = isoarea_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'A_IMAGE'] = a_image_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'B_IMAGE'] = b_image_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'THETA_IMAGE'] = theta_image_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'ERRA_IMAGE'] = erra_image_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'ERRB_IMAGE'] = errb_image_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'FWHM_IMAGE'] = fwhm_image_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'FLUX_RADIUS'] = flux_radius_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'ELONGATION'] = elongation_s.loc[i - 1]
            full_db.loc[full_db['SOURCE_NUMBER'] == i,
                        'ELLIPTICITY'] = ellipticity_s.loc[i - 1]

        if self.save:
            self.logger.debug('saving output to {}_3.csv'.format(filter_o_n))
            full_db.to_csv('{}_3.csv'.format(filter_o_n))

        return full_db

    def filter_pm(self, full_db, filter_o_n):
        """

        :param full_db:
        :param filter_o_n:
        :return:
        """
        self.logger.debug('after filtering detections')
        full_db = pm_filter(full_db, self.prfs_d['pm_low'],
                            self.prfs_d['pm_up'])

        if self.save:
            self.logger.debug('saving output to {}_4.csv'.format(filter_o_n))
            full_db.to_csv('{}_4.csv'.format(filter_o_n))

        full_db = sn_filter(full_db, self.prfs_d['pm_sn'])
        if self.save:
            self.logger.debug('saving output to {}_5.csv'.format(filter_o_n))
            full_db.to_csv('{}_5.csv'.format(filter_o_n))

        self.logger.debug('after proper motion')
        full_db = motion_filter(full_db, self.prfs_d['r_fit'])
        if self.save:
            self.logger.debug('saving output to {}_6.csv'.format(filter_o_n))
            full_db.to_csv('{}_6.csv'.format(filter_o_n))
