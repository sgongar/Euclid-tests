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

Todo:
    * Improve log messages

"""

from multiprocessing import Process

from astropy.io import fits
from astropy.table import Table

from misc import extract_settings, get_fits

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

        for idx in range(0, len(fits_files), 9):
            full_cat = '{}/f_{}_1.cat'.format(cats_dir, fits_files[idx][2:-5])
            fits_name = '{}/{}/{}.cat'.format(prfs_d['fits_dir'],
                                              folder_sex, fits_files[idx][:-5])

            compare_j = []
            for proc in range(0, 9, 1):
                idx_proc = proc + idx

                compare_p = Process(target=self.perform_analysis_thread,
                                    args=(full_cat, fits_name, idx_proc, ))
                compare_j.append(compare_p)
                compare_p.start()

            active_compare = list([job.is_alive() for job in compare_j])
            while True in active_compare:
                active_compare = list([job.is_alive() for job in compare_j])
                pass

    def perform_analysis_thread(self, full_cat, fits_name, idx):
        """

        """
        # open catalog
        # print "idx", idx
        cat_file = fits.open(full_cat)
        cat_data = Table(cat_file[2].data)
        cat_table = cat_data.to_pandas()

        fits_file = fits.open(fits_name)
        fits_data = Table(fits_file[2].data)
        fits_table = fits_data.to_pandas()

        print cat_table.columns

        for source_ in cat_table['SOURCE_NUMBER'].tolist():
                t = sex_table[sex_table['NUMBER'].isin([sex_source])].iloc[0]
                ra = t['X_WORLD']
                dec = t['Y_WORLD']

        # open image

        # compare


if __name__ == '__main__':
    test = Compare()
