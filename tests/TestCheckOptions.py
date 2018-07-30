#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script

Versions:

Problems:
    * setup has been ignored

Todo:
    * Improve log messages
    * 

"""

from os import getenv
from sys import argv, modules, path
from types import ModuleType

from unittest import TestCase, main
from mock import MagicMock, patch

home = getenv("HOME")
path.append('{}/build/sgongar/Euclid-tests/pipeline_elvis'.format(home))
path.append('{}/Dev/Euclid-tests/pipeline_elvis'.format(home))

from check_elvis import Check
from errors import BadSettings
from misc import extract_settings_elvis


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestCheckSuccessfulOptions(TestCase):
    """

    """

    def setup(self):
        pass

    @patch('check_elvis.Check.parameters')
    @patch('check_elvis.Check.full_pipeline')
    def test_full_option_chosen(self, parameters, full_pipeline):
        """

        :param parameters:
        :param full_pipeline:
        :return:
        """
        parameters.return_value = True
        full_pipeline.return_value = True
        argv[1] = '-full'
        test_check = Check()

        return self.assertTrue(test_check)

    @patch('check_elvis.Check.parameters')
    @patch('check_elvis.Check.clean')
    def test_clean_option_chosen(self, parameters, clean):
        """

        :param parameters:
        :param clean:
        :return:
        """
        parameters.return_value = True
        clean.return_value = True
        argv[1] = '-clean'
        test_check = Check()

        return self.assertTrue(test_check)

    @patch('check_elvis.Check.parameters')
    @patch('check_elvis.Check.split')
    def test_split_option_chosen(self, parameters, split):
        """

        :param parameters:
        :param split:
        :return:
        """
        parameters.return_value = True
        split.return_value = True
        argv[1] = '-split'
        test_check = Check()

        return self.assertTrue(test_check)

    @patch('check_elvis.Check.parameters')
    @patch('check_elvis.Check.sextractor')
    def test_sextractor_option_chosen(self, parameters, sextractor):
        """

        :param parameters:
        :param sextractor:
        :return:
        """
        parameters.return_value = True
        sextractor.return_value = True
        argv[1] = '-sextractor'
        test_check = Check()

        return self.assertTrue(test_check)

    @patch('check_elvis.Check.parameters')
    @patch('check_elvis.Check.scamp')
    def test_scamp_option_chosen(self, parameters, scamp):
        """

        :param parameters:
        :param scamp:
        :return:
        """
        parameters.return_value = True
        scamp.return_value = True
        argv[1] = '-scamp'
        test_check = Check()

        return self.assertTrue(test_check)

    @patch('check_elvis.Check.parameters')
    @patch('check_elvis.Check.filt')
    def test_filter_option_chosen(self, parameters, filt):
        """

        :param parameters:
        :param filt:
        :return:
        """
        parameters.return_value = True
        filt.return_value = True
        argv[1] = '-filter'
        test_check = Check()

        return self.assertTrue(test_check)

    @patch('check_elvis.Check.parameters')
    @patch('check_elvis.Check.restart')
    def test_scamp_option_chosen(self, parameters, restart):
        """

        :param parameters:
        :param restart:
        :return:
        """
        parameters.return_value = True
        restart.return_value = True
        argv[1] = '-restart'
        test_check = Check()

        return self.assertTrue(test_check)

    def teardrown(self):
        pass



# class TestCheckUnsuccessfulOptions(TestCase):
#     """
#
#     """
#
#     def setup(self):
#         pass
#
#     @patch('check_elvis.Check.parameters')
#     @patch('check_elvis.Check.full_pipeline')
#     def test_fulloptionchoosen(self, parameters, full_pipeline):
#         """
#
#         :param parameters:
#         :param full_pipeline:
#         :return:
#         """
#         parameters.return_value = True
#         full_pipeline.return_value = False
#         argv[1] = '-full'
#         test_check = Check()
#
#         return self.ra(test_check)
#
#     @patch('check_elvis.Check.parameters')
#     @patch('check_elvis.Check.split')
#     def test_splitoptionchoosen(self, parameters, split):
#         """
#
#         :param parameters:
#         :param split:
#         :return:
#         """
#         parameters.return_value = True
#         split.return_value = True
#         argv[1] = '-split'
#         test_check = Check()
#
#         return self.assertTrue(test_check)
#
#     def teardrown(self):
#         pass
