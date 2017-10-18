#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script

Versions:

Todo:
    * Improve log messages
    * astropy.io should be imported as a mocked module not a real one

"""

from sys import modules, path
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
        """
        define setUp parameters
        """
        path.append("/home/travis/build/sgongar/Euclid-tests/pipeline")

    # TODO Improve side_effect
    @patch('misc.extract_settings', side_effect=extract_settings_mock)
    @patch('misc.setting_logger', side_effect=setting_logger_mock)
    # @object(SatNetUI, 'initLogo', return_value=True)
    def test_NoOptionsPassed(self, extract_settings, setting_logger):
        """

        """
        from check import Check
        import misc

        self.assertRaises(IndexError, Check)

    def tearDown(self):
        """
        define tearDown parameters
        """
        pass


if __name__ == '__main__':
    main()
