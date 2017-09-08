#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

In order to improve legibilidad algunas abreviaturas han sido creadas
c = catalog
n = name
d = dictionary

Todo:
    * Improve log messages

"""

from multiprocessing import Process
from os import remove

from astropy.io import fits
from astropy.table import Table
from pandas import concat, DataFrame, read_csv

from cats_management import cut_catalog
from misc import extract_settings, get_fits, check_distance

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
"""
__license__ = "GPL"
"""
__version__ = "0.1"
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

        self.perform_analysis(prfs_d, folder_sex, folder_scmp)

    def perform_analysis(self, prfs_d, folder_sex, folder_scmp):
        """

        @param prfs_d:
        @param folder_sex:
        @param folder_scmp:

        @return True: if everything is alright
        """
        cats_dir = '{}/{}/{}'.format(prfs_d['catalogs_dir'], folder_sex,
                                     folder_scmp)

        fits_files = get_fits(unique=False)
        
        for idx in range(0, len(fits_files), 5):
            compare_j = []
            for proc in range(0, 5, 1):
                idx_proc = proc + idx
                if idx_proc < len(fits_files):
                    full_c = '{}/f_{}_1.cat'.format(cats_dir,
                                                      fits_files[idx_proc][2:-5])
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

            # Gets input data associated to selectes source
            cat_t = cat_table[cat_table['SOURCE_NUMBER'].isin([cat_source])]
            cat_ra = float(cat_t['ALPHA_J2000'])
            cat_dec = float(cat_t['DELTA_J2000'])

            # In order to improve comparation speed output catalog is reduced
            margins = [[cat_ra, 'ALPHA_J2000'], [cat_dec, 'DELTA_J2000']]
            margin = 2 * prfs_d['tolerance']
            fits_table_cut = cut_catalog(fits_table, margins, margin)

            # Takes a look over output sources
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
                    distances_cache.append(distance)

            if len(distances_cache) > 1 and close_flag == True:
                for distance_ in distances_cache:
                    sources_d['CCD'].append(fits_n[-12:-4])
                    sources_d['dither'].append(fits_n[-5:-4])
                    sources_d['distance'].append(distance_)
                idx_repeated += 1
            elif len(distances_cache) == 1 and close_flag == True:
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
        stats_df.to_csv('{}/{}.csv'.format(prfs_d['tmp_out'], fits_n[-12:-4]))

        for key_ in sources_d.keys():
            print key_, len(sources_d[key_])

        # Creates a DataFrame from sources dictionary
        sources_df = DataFrame(sources_d)
        sources_df.to_csv('{}/sources_{}.csv'.format(prfs_d['tmp_out'],
                                                     fits_n[-12:-4]))

    def populate_dict(self, stats_d, sources_d):
        """ populates dictionaries with selected keys

        @param stats_d:
        @param sources_d:

        @retun stats_d, sources_d
        """
        stats_d = {'CCD': [], 'dither': [], 'total': [], 'detected': [],
                   'repeated': [], 'lost': []}

        sources_d = {'CCD': [], 'dither': [], 'distance': [],
                     'i_ALPHA_J2000': [], "i_DELTA_J2000": [],
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
