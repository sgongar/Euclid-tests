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
from errors import CleanFailed


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class MockedLogger:
    def __init__(self, text):
        """

        :param text:
        """
        pass

    def info(self, text):
        """

        :param text:
        :return:
        """
        pass


class TestCleanMethodFromCheckElvis(TestCase):
    """

    """

    def setup(self):
        """

        :return:
        """
        pass

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('cosmic_elvis.CosmicELViS')
    def test_clean_works(self, CosmicELViS, setting_logger,
                         extract_settings_elvis):
        """

        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        CosmicELViS.return_value = True
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = {'fits_dir': 'fits_dir_mock',
                                               'fpas_dir': 'fpas_dir_mock'}

        argv[1] = '-clean'

        self.assertTrue(Check)

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('cosmic_elvis.CosmicELViS')
    def test_clean_fails(self, CosmicELViS, setting_logger,
                         extract_settings_elvis):
        """

        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        CosmicELViS.return_value = False
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = {'fits_dir': 'fits_dir_mock',
                                               'fpas_dir': 'fpas_dir_mock'}

        argv[1] = '-clean'

        self.assertRaises(CleanFailed, Check)

    def teardrown(self):
        """

        :return:
        """
        pass
