#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Python script to compare catalogs created with scamp against catalogs
    created with sextractor.

Versions:
* 0.1 - First version. Compare the catalog of each CCD with the same CCD in
        the single-star versions.
* 0.2 - Second version. Compare the catalog of each CCD with a catalog made
        with all the CCDs.
        This version is not correct as we don't have the correct header data.
* 0.3 - Third version.

In order to improve legibilidad algunas abreviaturas han sido creadas
* c = catalog
* n = name
* d = dictionary
* loc = location
* cmp = comparation

Todo:
    * Improve log messages
    * Improve docstring
    * Create versions history
    * License??
"""

from multiprocessing import Process
from os import remove

from astropy.io import fits
from astropy.table import Table
from pandas import concat, DataFrame, read_csv

from images_management import get_fits_limits
from cats_management import cut_catalog
from misc import extract_settings, get_fits, check_distance

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.3"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class Compare:

    def __init__(self):
        """ __init__ creates the basic variables for our study
            it takes nothing from outside

        """
        prfs_d = extract_settings()

        # Perfoms an analysis
        self.perform_analysis(prfs_d)

    def perform_analysis(self, prfs_d):
        """

        @param prfs_d:
        @param folder_sex:
        @param folder_scmp:

        @return True: if everything is alright
        """
        # Hardcoded
        folder_sex = '2_1.35_1.35_0.1_4'
        folder_scmp = '150_1.2_5_0.033'

        # Creates a dict
        cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'], folder_sex,
                                     folder_scmp)
        cats = self.load_catalogs(prfs_d, cats_dir)

        fits_files = get_fits(unique=False)

        print fits_files

        for idx in range(0, len(fits_files), 17):
            compare_j = []
            for proc in range(0, 17, 1):
                idx_p = proc + idx
                if idx_p < len(fits_files):
                    # fits_n is a sextracted catalog
                    fits_n = '{}/{}/{}.cat'.format(prfs_d['fits_dir'],
                                                   folder_sex,
                                                   fits_files[idx_p][:-5])

                    compare_p = Process(target=self.perform_analysis_thread,
                                        args=(fits_n, idx_p, prfs_d, cats,))
                    compare_j.append(compare_p)
                    compare_p.start()

            active_compare = list([job.is_alive() for job in compare_j])
            while True in active_compare:
                active_compare = list([job.is_alive() for job in compare_j])
                pass

        self.merge_stats(fits_files, prfs_d)

    def perform_analysis_thread(self, fits_n, idx, prfs_d, cats):
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
        ccd_borders_d = self.define_ccd_borders(prfs_d)

        stats_d['CCD'].append(fits_n[-12:-4])
        stats_d['dither'].append(fits_n[-5:-4])

        fits_file = fits.open(fits_n)
        fits_data = Table(fits_file[2].data)
        fits_table = fits_data.to_pandas()
        # Removing zero catalog references
        # fits_table = fits_table.loc[~fits_table['CATALOG_NUMBER'].isin([0])]
        stats_d['total'].append(len(fits_table['NUMBER'].tolist()))

        # All variables starting by 'cat' are referred to scamp output catalogs
        # All variables starting by 'fits' are referred to sextractor catalogs
        fits_sources = fits_table['NUMBER'].tolist()
        for fits_source in fits_sources:
            # Reset variables associated to each input source
            distances_cache = []
            i_alpha_cache = []
            i_delta_cache = []
            o_alpha_cache = []
            o_delta_cache = []

            # Gets input data associated to selectes source
            fits_t = fits_table[fits_table['NUMBER'].isin([fits_source])]
            fits_ra = float(fits_t['ALPHA_J2000'])
            fits_dec = float(fits_t['DELTA_J2000'])

            ccd_ = self.look_for_ccd(prfs_d, fits_ra, fits_dec, ccd_borders_d)

            if ccd_ is not None:
                cat_table = cats[ccd_]

                # To improve comparation speed output catalog is reduced
                margins = [[fits_ra, 'ALPHA_J2000'],
                           [fits_dec, 'DELTA_J2000']]
                margin = 2 * prfs_d['tolerance']
                cat_table_cut = cut_catalog(cat_table, margins, margin)

                # Takes a look over output sources
                # Each source has only one reference
                for cat_source in cat_table_cut['SOURCE_NUMBER'].tolist():
                    mask = cat_table_cut['SOURCE_NUMBER'].isin([cat_source])
                    cat_t = cat_table_cut[mask]
                    cat_ra = float(cat_t['ALPHA_J2000'])
                    cat_dec = float(cat_t['DELTA_J2000'])
                    # Compare cat_ra, cat_dec against fits_ra, fits_dec
                    close, distance = check_distance(fits_ra, cat_ra,
                                                     fits_dec, cat_dec,
                                                     prfs_d['tolerance'])
                    if close:
                        close_flag = True
                        distance = distance * 3600  # From degrees to seconds
                        distances_cache.append(distance)
                        i_alpha_cache.append(cat_ra)
                        i_delta_cache.append(cat_dec)
                        o_alpha_cache.append(fits_ra)
                        o_delta_cache.append(fits_dec)

                if len(distances_cache) > 1 and close_flag is True:
                    idx_cache = distances_cache.index(min(distances_cache))
                    sources_d['CCD'].append(fits_n[-12:-4])
                    sources_d['cat'].append(ccd_)
                    sources_d['dither'].append(fits_n[-5:-4])
                    distance_ = distances_cache[idx_cache]
                    sources_d['distance'].append(distance_)
                    sources_d['duplicated'].append(True)
                    i_alpha = i_alpha_cache[idx_cache]
                    sources_d['i_ALPHA_J2000'].append(i_alpha)
                    i_delta = i_delta_cache[idx_cache]
                    sources_d['i_DELTA_J2000'].append(i_delta)
                    o_alpha = o_alpha_cache[idx_cache]
                    sources_d['o_ALPHA_J2000'].append(o_alpha)
                    o_delta = o_delta_cache[idx_cache]
                    sources_d['o_DELTA_J2000'].append(o_delta)
                    idx_repeated += 1
                elif len(distances_cache) == 1 and close_flag is True:
                    idx_cache = distances_cache.index(min(distances_cache))
                    sources_d['CCD'].append(fits_n[-12:-4])
                    sources_d['cat'].append(ccd_)
                    sources_d['dither'].append(fits_n[-5:-4])
                    distance_ = distances_cache[idx_cache]
                    sources_d['distance'].append(distance_)
                    sources_d['duplicated'].append(False)
                    i_alpha = i_alpha_cache[idx_cache]
                    sources_d['i_ALPHA_J2000'].append(i_alpha)
                    i_delta = i_delta_cache[idx_cache]
                    sources_d['i_DELTA_J2000'].append(i_delta)
                    o_alpha = o_alpha_cache[idx_cache]
                    sources_d['o_ALPHA_J2000'].append(o_alpha)
                    o_delta = o_delta_cache[idx_cache]
                    sources_d['o_DELTA_J2000'].append(o_delta)
                    idx_detected += 1
                elif len(distances_cache) == 0:
                    sources_d['CCD'].append(fits_n[-12:-4])
                    sources_d['cat'].append(ccd_)
                    sources_d['dither'].append(fits_n[-5:-4])
                    distance_ = ''
                    sources_d['distance'].append(distance_)
                    sources_d['duplicated'].append('')
                    i_alpha = i_alpha_cache[idx_cache]
                    sources_d['i_ALPHA_J2000'].append(i_alpha)
                    i_delta = i_delta_cache[idx_cache]
                    sources_d['i_DELTA_J2000'].append(i_delta)
                    o_alpha = ''
                    sources_d['o_ALPHA_J2000'].append(o_alpha)
                    o_delta = ''
                    sources_d['o_DELTA_J2000'].append(o_delta)
                    idx_lost += 1
                else:
                    raise Exception
            else:
                idx_lost += 1

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
                          columns=['CCD', 'cat', 'dither', 'distance',
                                   'duplicated',
                                   'i_ALPHA_J2000', 'i_DELTA_J2000',
                                   'o_ALPHA_J2000', 'o_DELTA_J2000'])

    def define_ccd_borders(self, prfs_d):
        """

        @return ccd_borders_d:
        """

        fits_files = get_fits(unique=True)

        ccd_borders_d = {}
        for fits_ in fits_files:
            fits_loc = '{}/{}'.format(prfs_d['fits_dir'], fits_)
            fits_n = fits_[-13:-5]
            ccd_borders_d[fits_n] = get_fits_limits(fits_loc)

        # The arrangement of the elements in the dictionary makes it
        # difficult to use a DataFrame object.
        # ccd_borders = DataFrame.from_dict(ccd_borders_d)

        return ccd_borders_d

    def look_for_ccd(self, prfs_d, ra, dec, ccd_borders_d):
        """

        """
        for key_ in ccd_borders_d.keys():
            ccd_borders = ccd_borders_d[key_]
            ra_cmp = ccd_borders['below_ra'] < ra < ccd_borders['above_ra']
            dec_cmp = ccd_borders['below_dec'] < dec < ccd_borders['above_dec']

            if ra_cmp and dec_cmp:
                return key_

    def load_catalogs(self, prfs_d, cats_dir):
        """

        @param prfs_d:
        @param cats_dir:

        @retun cat_d;
        """
        cat_d = {}

        ccd_borders_d = self.define_ccd_borders(prfs_d)
        for ccd_ in ccd_borders_d.keys():
            full_c = '{}/f_20-21_{}_1.cat'.format(cats_dir, ccd_)

            # Opens catalog file
            cat_file = fits.open(full_c)
            cat_data = Table(cat_file[2].data)
            cat_table = cat_data.to_pandas()
            # Removing zero catalog references
            mask_zero = ~cat_table['CATALOG_NUMBER'].isin([0])
            cat_table = cat_table.loc[mask_zero]
            cat_d[ccd_] = cat_table

        return cat_d

    def populate_dict(self, stats_d, sources_d):
        """ populates dictionaries with selected keys

        @param stats_d:
        @param sources_d:

        @retun stats_d, sources_d
        """
        stats_d = {'CCD': [], 'dither': [],
                   'total': [], 'detected': [], 'repeated': [], 'lost': []}

        sources_d = {'CCD': [], 'dither': [], 'cat': [],
                     'distance': [], 'duplicated': [],
                     'i_ALPHA_J2000': [], 'i_DELTA_J2000': [],
                     'o_ALPHA_J2000': [], 'o_DELTA_J2000': []}

        return stats_d, sources_d

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
