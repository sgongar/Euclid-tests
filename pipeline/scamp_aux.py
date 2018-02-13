#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Versions:
- 0.1: Initial release.
- 0.2: Functions reorganised to two different classes.
- 0.3: New organization for filter Class
- 0.4: Filter now gets areas
- 0.5: check_source moved to another file

Todo:
    * Improve documentation
    * Speed-up areas process
    * sex_cf and scmp_cf should come from main file
"""

from subprocess import Popen

from astropy.io import fits
from astropy.table import Table
from pandas import concat, read_csv, Series

from misc import pm_compute, pm_filter, extract_settings, check_source
from misc import confidence_filter, create_folder, get_dither

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.5"
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

        self.logger = logger
        self.mag = mag
        self.scmp_d = scmp_d
        self.scmp_cf = scmp_cf
        self.sex_d = sex_d  # TODO

        self.scamp_process()

    def scamp_process(self):
        """

        :return:
        """
        sex_cf = '{}_{}_{}_{}_{}'.format(self.sex_d['deblend_nthresh'],
                                         self.sex_d['analysis_thresh'],
                                         self.sex_d['detect_thresh'],
                                         self.sex_d['deblend_mincount'],
                                         self.sex_d['detect_minarea'])

        self.logger.info('scamp process for magnitude {}'.format(self.mag))
        self.logger.info('Sextractor configuration: {}'.format(sex_cf))
        self.logger.info('Scamp configuration: {}'.format(self.scmp_cf))

        scmp_1 = 'scamp -c {}'.format(self.prfs_d['conf_scamp'])
        sex_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], self.mag,
                                         sex_cf)
        sex_output = 'mag_{}_CCD_x?_y?_d?.cat'.format(self.mag)
        scmp_2 = ' {}/{}'.format(sex_loc, sex_output)
        scmp_3 = ' -ASTREFCAT_NAME'
        cat_loc = '{}/{}/CCDs/{}'.format(self.prfs_d['fits_dir'], self.mag,
                                         sex_cf)
        cat_input = 'catalog_{}.cat'.format(self.mag)
        scmp_4 = ' {}/{}'.format(cat_loc, cat_input)
        scmp_5 = ' -PIXSCALE_MAXERR {}'.format(self.scmp_d['pixscale_maxerr'])
        scmp_6 = ' -POSANGLE_MAXERR {}'.format(self.scmp_d['posangle_maxerr'])
        scmp_7 = ' -POSITION_MAXERR {}'.format(self.scmp_d['position_maxerr'])
        scmp_8 = ' -CROSSID_RADIUS {}'.format(self.scmp_d['crossid_radius'])
        # Output catalogs location
        cats_dir = '{}/{}/{}/{}'.format(self.prfs_d['catalogs_dir'], self.mag,
                                        sex_cf, self.scmp_cf)
        merged_cat = '{}/merged_{}_{}.cat'.format(cats_dir, self.scmp_cf,
                                                  self.mag)
        scmp_9 = ' -MERGEDOUTCAT_NAME {}'.format(merged_cat)
        full_cat = '{}/full_{}_{}.cat'.format(cats_dir, self.scmp_cf, self.mag)
        scmp_10 = ' -FULLOUTCAT_NAME {}'.format(full_cat)
        scmp_p = scmp_1 + scmp_2 + scmp_3 + scmp_4 + scmp_5
        scmp_p = scmp_p + scmp_6 + scmp_7 + scmp_8 + scmp_9
        scmp_p = scmp_p + scmp_10

        create_folder(self.logger, cats_dir)

        process_scamp = Popen(scmp_p, shell=True)
        process_scamp.wait()

        self.logger.info('Scamp process finished.')

        return True


class ScampFilter:  # TODO Split scamp_filter method into single methods

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

        # Saves _1.csv
        # (merged_db, full_db, filter_o_n) = self.scamp_filter()
        # Saves _2.csv
        # full_db = self.compute_pm(merged_db, full_db, filter_o_n)

        # Filtered catalog dir
        filter_dir = '{}/{}/{}/{}'.format(self.prfs_d['filter_dir'], self.mag,
                                          self.sex_cf, self.scmp_cf)
        # Filtered catalog name
        filt_n = 'filt_{}_{}'.format(self.scmp_cf, self.mag)

        filter_o_n = '{}/{}'.format(filter_dir, filt_n)
        full_db = read_csv(('{}_2.csv'.format(filter_o_n)), index_col=0)

        full_db = self.get_areas(full_db, filter_o_n)
        # full_db = self.filter_pm(full_db, filter_o_n)
        # full_db = self.filter_class(full_db, filter_o_n)
        # self.filter_coherence(full_db, filter_o_n)

    def get_cat(self, catalog_n):
        """

        :param catalog_n:
        :return: cat_file
        """
        ccd, dither = get_dither(catalog_n)

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

        self.logger.debug('Scamp configuration: {}'.format(self.scmp_cf))
        self.logger.debug('Sextractor configuration: {}'.format(self.sex_cf))

        # Full catalog name
        full_n_dir = '{}/{}/{}/{}/'.format(self.prfs_d['catalogs_dir'],
                                           self.mag, self.sex_cf, self.scmp_cf)
        full_n_cat = 'full_{}_{}_1.cat'.format(self.scmp_cf, self.mag)
        full_n = '{}{}'.format(full_n_dir, full_n_cat)
        # Filtered catalog name
        filt_n = 'filt_{}_{}'.format(self.scmp_cf, self.mag)

        self.logger.debug('Opens full catalog')
        self.logger.debug('Dir: {}'.format(full_n_dir))
        self.logger.debug('Name: {}'.format(full_n_cat))
        full_cat = fits.open(full_n)
        full_db = Table(full_cat[2].data)
        self.logger.debug('Converts full Astropy catalog to Pandas format')
        full_db = full_db.to_pandas()

        # Getting merge catalog
        mrgd_n_dir = '{}/{}/{}/{}/'.format(self.prfs_d['catalogs_dir'],
                                           self.mag, self.sex_cf, self.scmp_cf)
        mrgd_n_cat = 'merged_{}_{}_1.cat'.format(self.scmp_cf, self.mag)
        mrgd_n = '{}{}'.format(mrgd_n_dir, mrgd_n_cat)

        self.logger.debug('Opens merged catalog')
        self.logger.debug('Dir: {}'.format(mrgd_n_dir))
        self.logger.debug('Name: {}'.format(mrgd_n_cat))
        merged_cat = fits.open(mrgd_n)
        self.logger.debug('Converts merged Astropy catalog to Pandas format')
        merged_db = Table(merged_cat[2].data)

        filter_o_n = '{}/{}'.format(filter_dir, filt_n)

        # Removing 0 catalog detections
        self.logger.debug('removes 0 catalog detections')
        full_db = full_db.loc[~full_db['CATALOG_NUMBER'].isin([0])]

        full_db = concat(g for _, g in full_db.groupby("SOURCE_NUMBER")
                         if len(g) >= int(self.prfs_d['detections']))

        if self.save:
            self.logger.debug('saves output to {}_1.csv'.format(filter_o_n))
            full_db.to_csv('{}_1.csv'.format(filter_o_n))

        return merged_db, full_db, filter_o_n

    def compute_pm(self, merged_db, full_db, filter_o_n):
        """

        :param merged_db:
        :param full_db:
        :param filter_o_n:
        :return: full_db
        """
        # Computing pm
        self.logger.debug('computing proper motion')
        full_db = pm_compute(self.logger, merged_db, full_db)
        if self.save:
            self.logger.debug('saves output to {}_2.csv'.format(filter_o_n))
            full_db.to_csv('{}_2.csv'.format(filter_o_n))

        return full_db

    def get_areas(self, full_db, filter_o_n):
        """

        :param full_db:
        :param filter_o_n:
        :return: full_db
        """
        self.logger.debug('runs areas filter')

        # Opens filtered file
        filter_cat = read_csv('{}_2.csv'.format(filter_o_n), index_col=0)
        # Gets unique sources from filtered file
        unique_sources = list(set(filter_cat['SOURCE_NUMBER'].tolist()))

        o_source_number_l = []
        o_catalog_number_l = []
        o_x_image_l = []
        o_y_image_l = []
        o_a_image_l = []
        o_erra_image_l = []
        o_b_image_l = []
        o_errb_image_l = []
        o_theta_image_l = []
        o_errtheta_image_l = []
        o_isoarea_l = []
        o_fwhm_image_l = []
        o_flux_radius_l = []
        o_mag_iso_l = []
        o_elongation_l = []
        o_ellipticity_l = []
        o_class_star_l = []
        o_alpha_j2000_l = []
        o_delta_j2000_l = []
        o_erra_world_l = []
        o_errb_world_l = []
        o_errtheta_world_l = []
        o_epoch_l = []
        o_mag_l = []
        o_magerr_l = []
        o_flags_extraction_l = []
        o_flags_scamp_l = []
        o_flags_ima_l = []
        o_pm_l = []
        o_pmerr_l = []
        o_pmalpha_l = []
        o_pmdelta_l = []
        o_pmalphaerr_l = []
        o_pmdeltaerr_l = []

        # Loops over unique sources of filtered file
        print(len(unique_sources))
        for idx, source_ in enumerate(unique_sources):
            print(idx)
            o_df = filter_cat[filter_cat['SOURCE_NUMBER'].isin([source_])]
            for i, row in enumerate(o_df.itertuples(), 1):
                source_number = row.SOURCE_NUMBER
                catalog_n = row.CATALOG_NUMBER
                cat_file = self.get_cat(catalog_n)
                x_image = row.X_IMAGE
                y_image = row.Y_IMAGE
                erra_image = row.ERRA_IMAGE
                errb_image = row.ERRB_IMAGE
                errtheta_image = row.ERRTHETA_IMAGE
                o_alpha = row.ALPHA_J2000
                o_delta = row.DELTA_J2000
                erra_world = row.ERRA_WORLD
                errb_world = row.ERRB_WORLD
                errtheta_world = row.ERRTHETA_WORLD
                epoch = row.EPOCH
                mag = row.MAG
                magerr = row.MAGERR
                flags_extraction = row.FLAGS_EXTRACTION
                flags_scamp = row.FLAGS_SCAMP
                flags_ima = row.FLAGS_IMA
                pm = row.PM
                pmerr = row.PMERR
                pmalpha = row.PMALPHA
                pmdelta = row.PMDELTA
                pmalphaerr = row.PMALPHAERR
                pmdeltaerr = row.PMDELTAERR

                # self.logger.debug('opening CCD catalog {}'.format(cat_file))
                ccd_cat = fits.open(cat_file)
                ccd_df = Table(ccd_cat[2].data)
                # self.logger.debug('converting CCD catalog to Pandas format')
                ccd_df = ccd_df.to_pandas()

                cat_df = check_source(ccd_df, o_alpha, o_delta)
                if cat_df.empty:
                    # FIXME problemas entre sextractor y scamp
                    # Como el objeto no esta en las imagenes originales de
                    # sextractor lo borro del catalogo de scamp
                    o_source_number_l.append('nan')
                    o_catalog_number_l.append('nan')
                    o_x_image_l.append('nan')
                    o_y_image_l.append('nan')
                    o_a_image_l.append('nan')
                    o_erra_image_l.append('nan')
                    o_b_image_l.append('nan')
                    o_errb_image_l.append('nan')
                    o_theta_image_l.append('nan')
                    o_errtheta_image_l.append('nan')
                    o_alpha_j2000_l.append('nan')
                    o_delta_j2000_l.append('nan')
                    o_isoarea_l.append('nan')
                    o_fwhm_image_l.append('nan')
                    o_flux_radius_l.append('nan')
                    o_mag_iso_l.append('nan')
                    o_elongation_l.append('nan')
                    o_ellipticity_l.append('nan')
                    o_class_star_l.append('nan')
                    o_erra_world_l.append('nan')
                    o_errb_world_l.append('nan')
                    o_errtheta_world_l.append('nan')
                    o_epoch_l.append('nan')
                    o_mag_l.append('nan')
                    o_magerr_l.append('nan')
                    o_flags_extraction_l.append('nan')
                    o_flags_scamp_l.append('nan')
                    o_flags_ima_l.append('nan')
                    o_pm_l.append('nan')
                    o_pmerr_l.append('nan')
                    o_pmalpha_l.append('nan')
                    o_pmdelta_l.append('nan')
                    o_pmalphaerr_l.append('nan')
                    o_pmdeltaerr_l.append('nan')
                else:
                    o_source_number_l.append(source_number)
                    o_catalog_number_l.append(catalog_n)
                    o_x_image_l.append(x_image)
                    o_y_image_l.append(y_image)
                    a_image = cat_df['A_IMAGE'].iloc[0]
                    o_a_image_l.append(a_image)
                    o_erra_image_l.append(erra_image)
                    b_image = cat_df['B_IMAGE'].iloc[0]
                    o_b_image_l.append(b_image)
                    o_errb_image_l.append(errb_image)
                    theta_image = cat_df['THETA_IMAGE'].iloc[0]
                    o_theta_image_l.append(theta_image)
                    o_errtheta_image_l.append(errtheta_image)
                    o_alpha_j2000_l.append(o_alpha)
                    o_delta_j2000_l.append(o_delta)
                    isoarea_image = cat_df['ISOAREA_IMAGE'].iloc[0]
                    o_isoarea_l.append(isoarea_image)
                    fwhm_image = cat_df['FWHM_IMAGE'].iloc[0]
                    o_fwhm_image_l.append(fwhm_image)
                    flux_radius = cat_df['FLUX_RADIUS'].iloc[0]
                    o_flux_radius_l.append(flux_radius)
                    mag_iso = cat_df['MAG_ISO'].iloc[0]
                    o_mag_iso_l.append(mag_iso)
                    elongation = cat_df['ELONGATION'].iloc[0]
                    o_elongation_l.append(elongation)
                    ellipticity = cat_df['ELLIPTICITY'].iloc[0]
                    o_ellipticity_l.append(ellipticity)
                    class_star = cat_df['CLASS_STAR'].iloc[0]
                    o_class_star_l.append(class_star)
                    o_erra_world_l.append(erra_world)
                    o_errb_world_l.append(errb_world)
                    o_errtheta_world_l.append(errtheta_world)
                    o_epoch_l.append(epoch)
                    o_mag_l.append(mag)
                    o_magerr_l.append(magerr)
                    o_flags_extraction_l.append(flags_extraction)
                    o_flags_scamp_l.append(flags_scamp)
                    o_flags_ima_l.append(flags_ima)
                    o_pm_l.append(pm)
                    o_pmerr_l.append(pmerr)
                    o_pmalpha_l.append(pmalpha)
                    o_pmdelta_l.append(pmdelta)
                    o_pmalphaerr_l.append(pmalphaerr)
                    o_pmdeltaerr_l.append(pmdeltaerr)

        o_source_number_s = Series(o_source_number_l, name='SOURCE_NUMBER')
        o_catalog_number_s = Series(o_catalog_number_l, name='CATALOG_NUMBER')
        o_x_image_s = Series(o_x_image_l, name='X_IMAGE')
        o_y_image_s = Series(o_y_image_l, name='Y_IMAGE')
        o_a_image_s = Series(o_a_image_l, name='A_IMAGE')
        o_erra_image_s = Series(o_erra_image_l, name='ERRA_IMAGE')
        o_b_image_s = Series(o_b_image_l, name='B_IMAGE')
        o_errb_image_s = Series(o_errb_image_l, name='ERRB_IMAGE')
        o_theta_image_s = Series(o_theta_image_l, name='THETA_IMAGE')
        o_errtheta_image_s = Series(o_errtheta_image_l, name='ERRTHETA_IMAGE')
        o_alpha_j2000_s = Series(o_alpha_j2000_l, name='ALPHA_J2000')
        o_delta_j2000_s = Series(o_delta_j2000_l, name='DELTA_J2000')
        o_isoarea_s = Series(o_isoarea_l, name='ISOAREA_IMAGE')
        o_fwhm_image_s = Series(o_fwhm_image_l, name='FWHM_IMAGE')
        o_flux_radius_s = Series(o_flux_radius_l, name='FLUX_RADIUS')
        o_mag_iso_s = Series(o_mag_iso_l, name='MAG_ISO')
        o_elongation_s = Series(o_elongation_l, name='ELONGATION')
        o_ellipticity_s = Series(o_ellipticity_l, name='ELLIPTICITY')
        o_class_star_s = Series(o_class_star_l, name='CLASS_STAR')
        o_erra_world_s = Series(o_erra_world_l, name='ERRA_WORLD')
        o_errb_world_s = Series(o_errb_world_l, name='ERRB_WORLD')
        o_errtheta_world_s = Series(o_errtheta_world_l, name='ERRTHETA_WORLD')
        o_epoch_s = Series(o_epoch_l, name='EPOCH')
        o_mag_s = Series(o_mag_l, name='MAG')
        o_magerr_s = Series(o_magerr_l, name='MAGERR')
        o_flags_extraction_s = Series(o_flags_extraction_l,
                                      name='FLAGS_EXTRACTION')
        o_flags_scamp_s = Series(o_flags_scamp_l, name='FLAGS_SCAMP')
        o_flags_ima_s = Series(o_flags_ima_l, name='FLAGS_IMA')
        o_pm_s = Series(o_pm_l, name='PM')
        o_pmerr_s = Series(o_pmerr_l, name='PMERR')
        o_pmalpha_s = Series(o_pmalpha_l, name='PMALPHA')
        o_pmdelta_s = Series(o_pmdelta_l, name='PMDELTA')
        o_pmalphaerr_s = Series(o_pmalphaerr_l, name='PMALPHAERR')
        o_pmdeltaerr_s = Series(o_pmdeltaerr_l, name='PMDELTAERR')

        full_db = concat([o_source_number_s, o_catalog_number_s, o_x_image_s,
                          o_y_image_s, o_a_image_s, o_erra_image_s,
                          o_b_image_s, o_errb_image_s, o_theta_image_s,
                          o_errtheta_image_s, o_alpha_j2000_s, o_delta_j2000_s,
                          o_isoarea_s, o_fwhm_image_s, o_flux_radius_s,
                          o_mag_iso_s, o_elongation_s, o_ellipticity_s,
                          o_class_star_s, o_erra_world_s, o_errb_world_s,
                          o_errtheta_world_s, o_epoch_s, o_mag_s, o_magerr_s,
                          o_flags_extraction_s, o_flags_scamp_s, o_flags_ima_s,
                          o_pm_s, o_pmerr_s, o_pmalpha_s, o_pmdelta_s,
                          o_pmalphaerr_s, o_pmdeltaerr_s], axis=1)

        if self.save:
            self.logger.debug('saves output to {}_3.csv'.format(filter_o_n))
            full_db.to_csv('{}_3.csv'.format(filter_o_n))

        return full_db

    def filter_pm(self, full_db, filter_o_n):
        """

        :param full_db:
        :param filter_o_n:
        :return: full_db
        """
        self.logger.debug('runs proper motion filter')
        full_db = pm_filter(full_db, self.prfs_d['pm_low'],
                            self.prfs_d['pm_up'])

        if self.save:
            self.logger.debug('saves output to {}_4.csv'.format(filter_o_n))
            full_db.to_csv('{}_4.csv'.format(filter_o_n))

        return full_db
    
    def filter_class(self, full_db, filter_o_n):
        """

        :param full_db:
        :param filter_o_n:
        :return:
        """
        self.logger.debug('runs class star filter')
        full_db = full_db[full_db['CLASS_STAR'] > 0.95]

        if self.save:
            self.logger.debug('saves output to {}_5.csv'.format(filter_o_n))
            full_db.to_csv('{}_5.csv'.format(filter_o_n))

        return full_db

    def filter_coherence(self, full_db, filter_o_n):
        """

        :param full_db:
        :param filter_o_n:
        :return: full_db
        """
        self.logger.debug('runs coherence motion filter')
        full_db = confidence_filter(full_db, self.prfs_d['r_fit'])
        if self.save:
            self.logger.debug('saves output to {}_6.csv'.format(filter_o_n))
            full_db.to_csv('{}_6.csv'.format(filter_o_n))
