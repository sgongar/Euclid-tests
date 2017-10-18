#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script

Versions:

Todo:
    * Improve log messages
    * 

"""

from os import getenv
from sys import argv, modules, path
from types import ModuleType

from unittest import TestCase, main
from mock import MagicMock, patch




__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def extract_settings_mock():
    """
    """
    return True

def setting_logger_mock():
    """
    """
    return True

class TestCheckOptions(TestCase):
    """

    """

    def setUp(self):
        print "helllloooo"
        home = getenv("HOME")
        path.append('{}/build/sgongar/Euclid-tests/pipeline'.format(home))
        # path.append('/mnt/g/dev/Euclid-tests/pipeline')

        from check import Check
        from errors import BadSettings
        import misc

    @patch('misc.extract_settings', side_effect=extract_settings_mock)
    @patch('misc.setting_logger', side_effect=setting_logger_mock)
    @patch.object(Check, 'full_pipeline', return_value=True)
    def test_fulloptionchoosen(self, extract_settings, setting_logger,
                               full_pipeline):
        """

        """
        argv[1] = '-full'

        return self.assertTrue(Check)

    def tearDrown(self):
        pass

"""
class TestCheckNoOptions(TestCase):

    @patch('misc.extract_settings', side_effect=extract_settings_mock)
    @patch('misc.setting_logger', side_effect=setting_logger_mock)
    def test_nooptionspassed(self, extract_settings, setting_logger):
        argv[1] = '-wrong'

        return self.assertRaises(BadSettings, Check)
"""
