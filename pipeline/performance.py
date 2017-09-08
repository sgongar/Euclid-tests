#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

This module perfoms a test over sextractor output.

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
from os import listdir

from pandas import concat, Series

from images_management import get_fits_limits
from misc import create_configurations, get_cats
from regions import Create_regions


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


class SextractorPerformance:

    def __init__(self, logger, prfs_d):
        
        if not self.check(logger, prfs_d):
            raise Exception

    def check(self, logger, prfs_d):
        """

        @param logger:
        @param prfs_d:

        @return True: if everything goes alright.
        """
        cat_dict = get_cats(logger)
        mode = {'type': 'sextractor'}
        confs, total_confs = create_configurations(logger, prfs_d, mode)

        # Busca por todas las configuraciones posibles
        for conf_idx, conf_ in enumerate(confs):
            # Longitude of each folder
            conf_long = len(cat_dict[cat_dict.keys()[conf_idx]])

            for cat_idx in range(0, conf_long, prfs_d['cores_number']):
                check_j = []  # A list for all the works for each configuration
                for idx in range(0, prfs_d['cores_number'], 1):
                    catalog = cat_dict[cat_dict.keys()[conf_idx]][idx + cat_idx]
                    dither = catalog[-5:-4]

                    image = '{}.fits'.format(catalog[:-3])
                    ssos = self.get_ssos(logger, prfs_d, image, dither)

                    check_p = Process(target=self.check_thread,
                                      args=(logger, prfs_d, catalog, ssos,))
                    check_j.append(check_p)
                    check_p.start()

                active_check = list([job.is_alive() for job in check_j])
                while True in active_check:
                    active_check = list([job.is_alive() for job in check_j])
                    pass

        return True

    def check_thread(self, logger, prfs_d, catalog, ssos):
        """

        @param logger:
        @param prfs_d:
        @param catalog:

        @return True: if everything goes alright.
        """

        """
        # Creates a list populated with all single sources
        unique_sources = list(set(ssos['source'].tolist()))

        # Creates lists for final catalogue in out_dict
        # TODO out_dict should be a shared dict

        counter_source = 0
        # Loops over CCD sources - Custom catalog
        for source_number in unique_sources:

            logger.debug('look for source {}'.format(source_number))
            counter_source = counter_source + 1
            tmp_cat = ssos.loc[ssos['source'] == source_number]
            # Looks for duplicate sources in same dither
            if tmp_cat['source'].size is not tmp_cat['dither_values'].size:
                raise Exception

            logger.debug('checking {} of {}/{}'.format(source_number,
                                                       counter_source,
        """

    def get_ssos(self, logger, prfs_d, image, dither):
        """

        @param logger:
        @param prfs_d:
        @param image:
        @param dither:

        @return ssos:
        """
        save = False
        complete = True

        # Gets sky limits of each image
        d_limits = {}
        for d in range(1, 5, 1):
            i_image = '{}/{}{}{}'.format(prfs_d['fits_dir'], image[:15],
                                         d, image[-5:])
            d_limits['{:d}'.format(d)] = get_fits_limits(i_image)

        # Gets ra/dec coordinates for each dither
        # TODO Get only coordinates for selected dither
        d_ssos = {}
        for d in range(1, 5, 1):
            d_cat = prfs_d['input_cats'] + '/Cat_20-21_d{}.dat'.format(d)
            d_ssos['{}'.format(d)] = Create_regions(d_cat, prfs_d).luca(save,
                                                                        complete)

        cat_list = []
        # Substract not present SSOs in image selected
        limits = d_limits[dither]  # get fits limits
        out_cat = d_ssos['{}'.format(dither)]
        out_cat = out_cat[limits['above_ra'] > out_cat['alpha_j2000']]
        out_cat = out_cat[limits['below_ra'] < out_cat['alpha_j2000']]
        out_cat = out_cat[limits['above_dec'] > out_cat['delta_j2000']]
        out_cat = out_cat[limits['below_dec'] < out_cat['delta_j2000']]
        
        cat_list.append(out_cat)

        # Merge all catalog dither into a single one
        ssos = concat(cat_list, ignore_index=True)

        return ssos