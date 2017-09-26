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
* 0.4 - 
* 0.5 - Creates regions for input and output catalogs.

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
from cats_management import cut_catalog, shrink_catalog
from misc import extract_settings, get_fits, check_distance

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.5"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class Compare:

    def __init__(self):
        """ __init__ creates the basic variables for our study
            it takes nothing from outside

        """
        prfs_d = extract_settings()

        # Perfoms an analysis TODO Create a custom Exception!
        self.perform_analysis(prfs_d)

    def perform_analysis(self, prfs_d):
        """ This method obtains the catalogues under study and launches
        different threads for each pair of catalogues.

        TODO For now location is hard-coded. In posterior versions of this
        script this should be fixed in order to allow a more complete analysis
        along different sextractor and scamp configurations.

        @param prfs_d: a regular preferences dictionary.

        @return True: if everything is alright.
        """
        # Hardcoded references.
        folder_sex = '20_100_100_0.1_4'
        # folder_scmp = '150_1.2_5_0.033'

        # Load scamp's catalogs.
        # scamp_cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'],
        #                                    folder_sex, folder_scmp)
        # scamp_cats = self.load_scamp_catalogs(prfs_d, scamp_cats_dir)

        # Load sextractor's catalogs. Option to filter
        # sex_cats_dir = '{}/{}'.format(prfs_d['fits_dir'], folder_sex)
        # sex_cats = self.load_sex_catalogs(prfs_d, sex_cats_dir, unique=True)

        # Load sextractor's reference catalogs.
        sex_cats_dir = '{}/{}'.format(prfs_d['fits_dir'], folder_sex)
        sex_cats_i = self.load_sex_catalogs(prfs_d, sex_cats_dir,
                                            unique=True, filter=False)

        # Load sextractor's full populated catalogs.
        fits_dir = '/media/sf_CarpetaCompartida/luca_data/v14/CCDs'
        sex_cats_dir = '{}/{}'.format(fits_dir, folder_sex)
        sex_cats_o = self.load_sex_catalogs(prfs_d, sex_cats_dir,
                                            unique=True, filter=False)

        # Unique flag allow us to get a list of references populated only by
        # images of first dither.
        fits_files = get_fits(unique=True)

        # In order to improve script reuse, references to catalogues
        # will be dynamic.
        input_cats = [sex_cats_i, ['NUMBER', 'ALPHA_J2000', 'DELTA_J2000']]
        output_cats = [sex_cats_o, ['NUMBER', 'ALPHA_J2000', 'DELTA_J2000']]

        # TODO '17' value is harcoded to this (24 cores/48 threads) CPU.
        # Should be change to a relative value related to microprocessor used.
        for idx in range(0, len(fits_files), 17):
            compare_j = []  # A list for active jobs.
            for proc in range(0, 17, 1):
                idx_p = proc + idx
                if idx_p < len(fits_files):
                    # fits_n is a sextracted catalog
                    fits_n = '{}/{}/{}.cat'.format(prfs_d['fits_dir'],
                                                   folder_sex,
                                                   fits_files[idx_p][:-5])

                    compare_p = Process(target=self.perform_analysis_thread,
                                        args=(fits_n, prfs_d,
                                              input_cats, output_cats,))
                    compare_j.append(compare_p)
                    compare_p.start()

            active_compare = list([job.is_alive() for job in compare_j])
            while True in active_compare:
                active_compare = list([job.is_alive() for job in compare_j])
                pass

        # TODO Create a directory if this doesn't exist  
        # Creates a file with all the statistics obtained.
        if not self.merge_stats(fits_files, prfs_d):
            raise Exception  # TODO Creates a custom Exception.

        return True

    def perform_analysis_thread(self, fits_n, prfs_d, input_cats, output_cats):
        """

        @param fits_n:
        @param prfs_d: A dictionary with all script preferences.
        @param input_cats:
        @param output_cats:
        """
        # Set-up indexes for this particular case
        idx_detected = 0
        idx_repeated = 0
        idx_lost = 0
        # Set-up dictionaries for this particular case
        stats_d = {}
        sources_d = {}
        i_comp_d = {}
        o_comp_d = {}
        (stats_d, sources_d,
         i_comp_d, o_comp_d) = self.populate_dict(stats_d, sources_d,
                                                  i_comp_d, o_comp_d)

        ccd_borders_d = self.define_ccd_borders(prfs_d)

        distorsion = 0.05  # Set-up variable for distorsion allowed.

        stats_d['CCD'].append(fits_n[-12:-4])
        stats_d['dither'].append(fits_n[-5:-4])

        input_table = input_cats[0][fits_n[-12:-4]]
        # Removing zero catalog references
        # fits_table = fits_table.loc[~fits_table['CATALOG_NUMBER'].isin([0])]
        stats_d['total'].append(len(input_table['NUMBER'].tolist()))

        # All variables starting by 'cat' are referred to scamp output catalogs
        # All variables starting by 'fits' are referred to sextractor catalogs
        input_sources = input_table['NUMBER'].tolist()
        for input_source in input_sources:
            # Reset variables associated to each input source
            distances_cache = []
            i_alpha_cache = []
            i_delta_cache = []
            o_alpha_cache = []
            o_delta_cache = []

            i_flux_iso_cache = []
            i_flux_auto_cache = []
            o_flux_iso_cache = []
            o_flux_auto_cache = []

            # Gets input data associated to selectes source
            input_t = input_table[input_table['NUMBER'].isin([input_source])]
            input_ra = float(input_t['ALPHA_J2000'])
            input_dec = float(input_t['DELTA_J2000'])
            input_flux_iso = float(input_t['FLUX_ISO'])
            input_flux_auto = float(input_t['FLUX_AUTO'])

            ccd_ = self.look_for_ccd(prfs_d, input_ra, input_dec,
                                     ccd_borders_d)

            if ccd_ is not None:
                output_table = output_cats[0][ccd_]

                # To improve comparation speed output catalog is reduced
                margins = [[input_ra, 'ALPHA_J2000'],
                           [input_dec, 'DELTA_J2000']]
                margin = 2 * prfs_d['tolerance']
                output_table_cut = cut_catalog(output_table, margins, margin)

                # Takes a look over output sources
                # Each source has only one reference
                for output_source in output_table_cut['NUMBER'].tolist():
                    mask = output_table_cut['NUMBER'].isin([output_source])
                    output_t = output_table_cut[mask]
                    output_ra = float(output_t['ALPHA_J2000'])
                    output_dec = float(output_t['DELTA_J2000'])
                    output_flux_iso = float(output_t['FLUX_ISO'])
                    output_flux_auto = float(output_t['FLUX_AUTO'])
                    # Compare output_ra, output_dec against input_ra, input_dec
                    close, distance = check_distance(input_ra, output_ra,
                                                     input_dec, output_dec,
                                                     prfs_d['tolerance'])
                    if close:
                        close_flag = True
                        distance = distance * 3600  # From degrees to seconds
                        distances_cache.append(distance)
                        i_alpha_cache.append(input_ra)
                        i_delta_cache.append(input_dec)
                        o_alpha_cache.append(output_ra)
                        o_delta_cache.append(output_dec)
                        i_flux_iso_cache.append(input_flux_iso)
                        i_flux_auto_cache.append(input_flux_auto)
                        o_flux_iso_cache.append(output_flux_iso)
                        o_flux_auto_cache.append(output_flux_auto)

                if len(distances_cache) > 1 and close_flag is True:
                    # These values will be the same whatever the output is.
                    sources_d['CCD'].append(fits_n[-12:-4])
                    sources_d['cat'].append(ccd_)
                    sources_d['dither'].append(fits_n[-5:-4])

                    idx_cache = distances_cache.index(min(distances_cache))
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

                    i_flux_iso = i_flux_iso_cache[idx_cache]
                    sources_d['i_FLUX_ISO'].append(i_flux_iso)
                    i_flux_auto = i_flux_auto_cache[idx_cache]
                    sources_d['i_FLUX_AUTO'].append(i_flux_auto)
                    o_flux_iso = o_flux_iso_cache[idx_cache]
                    sources_d['o_FLUX_ISO'].append(o_flux_iso)
                    o_flux_auto = o_flux_auto_cache[idx_cache]
                    sources_d['o_FLUX_AUTO'].append(o_flux_auto)

                    idx_repeated += 1
                elif len(distances_cache) == 1 and close_flag is True:
                    # These values will be the same whatever the output is.
                    sources_d['CCD'].append(fits_n[-12:-4])
                    sources_d['cat'].append(ccd_)
                    sources_d['dither'].append(fits_n[-5:-4])

                    idx_cache = distances_cache.index(min(distances_cache))
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

                    i_flux_iso = i_flux_iso_cache[idx_cache]
                    sources_d['i_FLUX_ISO'].append(i_flux_iso)
                    i_flux_auto = i_flux_auto_cache[idx_cache]
                    sources_d['i_FLUX_AUTO'].append(i_flux_auto)
                    o_flux_iso = o_flux_iso_cache[idx_cache]
                    sources_d['o_FLUX_ISO'].append(o_flux_iso)
                    o_flux_auto = o_flux_auto_cache[idx_cache]
                    sources_d['o_FLUX_AUTO'].append(o_flux_auto)

                    idx_detected += 1
                elif len(distances_cache) == 0:
                    idx_lost += 1

                if distance_ > distorsion and len(distances_cache) >= 1:
                    idx_cache = distances_cache.index(min(distances_cache))
                    i_alpha = i_alpha_cache[idx_cache]
                    i_comp_d['i_ALPHA_J2000'].append(i_alpha)
                    i_delta = i_delta_cache[idx_cache]
                    i_comp_d['i_DELTA_J2000'].append(i_delta)
                    o_alpha = o_alpha_cache[idx_cache]
                    o_comp_d['o_ALPHA_J2000'].append(o_alpha)
                    o_delta = o_delta_cache[idx_cache]
                    o_comp_d['o_DELTA_J2000'].append(o_delta)

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
                                   'o_ALPHA_J2000', 'o_DELTA_J2000',
                                   'i_FLUX_ISO', 'i_FLUX_AUTO',
                                   'o_FLUX_ISO', 'o_FLUX_AUTO'])

        i_comp_df = DataFrame(i_comp_d)
        i_comp_df.to_csv('{}/i_regions_{}.reg'.format(prfs_d['tmp_out'],
                                                      fits_n[-12:-4]),
                         columns=['i_ALPHA_J2000', 'i_DELTA_J2000'],
                         index=False, header=False, sep=" ")

        o_comp_df = DataFrame(o_comp_d)
        o_comp_df.to_csv('{}/o_regions_{}.reg'.format(prfs_d['tmp_out'],
                                                      fits_n[-12:-4]),
                         columns=['o_ALPHA_J2000', 'o_DELTA_J2000'],
                         index=False, header=False, sep=" ")

    def define_ccd_borders(self, prfs_d):
        """

        @param prfs_d: A dictionary with all script preferences.

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

        @param prfs_d: A dictionary with all script preferences.
        @param ra:
        @param dec:
        @param ccd_borders_d:

        @return key_:
        """
        for key_ in ccd_borders_d.keys():
            ccd_borders = ccd_borders_d[key_]
            ra_cmp = ccd_borders['below_ra'] < ra < ccd_borders['above_ra']
            dec_cmp = ccd_borders['below_dec'] < dec < ccd_borders['above_dec']

            if ra_cmp and dec_cmp:
                return key_

    def load_scamp_catalogs(self, prfs_d, cats_dir):
        """

        @param prfs_d: A dictionary with all script preferences.
        @param cats_dir:

        @retun cat_d;
        """
        cat_d = {}  # A dictionary for catalogs.

        # Get keys for catalogs.
        ccd_borders_d = self.define_ccd_borders(prfs_d)
        # A for loop for open and filter catalogs.
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

    def load_sex_catalogs(self, prfs_d, sex_cats_dir, unique, filter):
        """

        @param prfs_d: A dictionary with all script preferences.
        @param sex_cats_dir: A string with the catalog location.
        @param unique: If active, only images from the first dither
                       will be used.
        @param filter: If active catalogs will be shrink to dither one limits.

        @return cat_d: A dictionary populated by sextractor's catalogs.
        """
        cat_d = {}  # A dictionary for catalogs.

        if filter:
            margins_d = self.get_catalog_limits(prfs_d, unique)

        fits_files = get_fits(unique=True)
        for fits_ in fits_files:
            print "fits_", fits_[:-5]
            sex_c = '{}/{}.cat'.format(sex_cats_dir, fits_[:-5])
            print "sex_c", sex_c

            cat_file = fits.open(sex_c)
            cat_data = Table(cat_file[2].data)
            cat_table = cat_data.to_pandas()

            if filter:
                cat_d[fits_[-13:-5]] = shrink_catalog(cat_table, margins_d)
            else:
                cat_d[fits_[-13:-5]] = cat_table

        return cat_d

    def get_catalog_limits(self, prfs_d, unique):
        """
        @param prfs_d:
        @param unique:

        @return margins_d:
        """
        limits_d = {}  # A dictionay for limits.

        # Populates dictionary with lists
        limits_d['above_ra'] = []
        limits_d['below_ra'] = []
        limits_d['above_dec'] = []
        limits_d['below_dec'] = []

        # Gets fits files references
        fits_files = get_fits(unique=True)

        # Get limits for dither one.
        for fits_ in fits_files:
            fits_loc = '{}/{}'.format(prfs_d['fits_dir'], fits_)
            limits = get_fits_limits(fits_loc)
            limits_d['above_ra'].append(limits['above_ra'])
            limits_d['below_ra'].append(limits['below_ra'])
            limits_d['above_dec'].append(limits['above_dec'])
            limits_d['below_dec'].append(limits['below_dec'])

        above_ra = max(limits_d['above_ra'])
        below_ra = min(limits_d['below_ra'])
        above_dec = max(limits_d['above_dec'])
        below_dec = min(limits_d['below_dec'])

        margins_d = {'above_ra': above_ra, 'below_ra': below_ra,
                     'above_dec': above_dec, 'below_dec': below_dec}

        return margins_d

    def populate_dict(self, stats_d, sources_d, i_comp_d, o_comp_d):
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
                     'o_ALPHA_J2000': [], 'o_DELTA_J2000': [],
                     'i_FLUX_ISO': [], 'i_FLUX_AUTO': [],
                     'o_FLUX_ISO': [], 'o_FLUX_AUTO': []}

        i_comp_d = {'i_ALPHA_J2000': [], 'i_DELTA_J2000': []}
        o_comp_d = {'o_ALPHA_J2000': [], 'o_DELTA_J2000': []}

        return stats_d, sources_d, i_comp_d, o_comp_d

    def merge_stats(self, fits_files, prfs_d):
        """

        @param fits_files:
        @param prfs_d:

        @return True: if everything goes alright.
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

        return True


if __name__ == '__main__':
    test = Compare()
