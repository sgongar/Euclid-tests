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

import os
import sys

from unittest import TestCase, main
from mock import MagicMock, patch

# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import check_elvis
from errors import BadSettings
import misc


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


class TestCheckSuccessfulOptions(TestCase):
    """

    """

    def setUp(self):
        """

        :return:
        """
        pass

    def test_full_option_chosen(self):
        """

        :param full_pipeline:
        :param extract_settings_elvis:
        :return:
        """
        misc.extract_settings_elvis = MagicMock(return_value=True)
        check_elvis.Check.full_pipeline = MagicMock(return_value=True)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)

        sys.argv.append('-full')

        return self.assertTrue(check_elvis.Check())

    def test_clean_option_chosen(self):
        """

        :return:
        """
        misc.extract_settings_elvis = MagicMock(return_value=True)
        check_elvis.Check.clean = MagicMock(return_value=True)
        misc.setting_logger = MagicMock(side_effect=MockedLogger)

        sys.argv.append('-clean')

        return self.assertTrue(check_elvis.Check())
    #
    # @patch('misc.extract_settings_elvis')
    # @patch('misc.setting_logger')
    # @patch('check_elvis.Check.split')
    # def test_split_option_chosen(self, split, setting_logger,
    #                              extract_settings_elvis):
    #     """
    #
    #     :param split:
    #     :param setting_logger:
    #     :param extract_settings_elvis:
    #     :return:
    #     """
    #     split.return_value = True
    #     setting_logger.side_effect = MockedLogger
    #     extract_settings_elvis.return_value = True
    #
    #     sys.argv[1] = '-split'
    #
    #     return self.assertTrue(Check())
    #
    # @patch('misc.extract_settings_elvis')
    # @patch('misc.setting_logger')
    # @patch('check_elvis.Check.sextractor')
    # def test_sextractor_option_chosen(self, sextractor, setting_logger,
    #                                   extract_settings_elvis):
    #     """
    #
    #     :param sextractor:
    #     :param setting_logger:
    #     :param extract_settings_elvis:
    #     :return:
    #     """
    #     sextractor.return_value = True
    #     setting_logger.side_effect = MockedLogger
    #     extract_settings_elvis.return_value = True
    #
    #     sys.argv[1] = '-sextractor'
    #
    #     return self.assertTrue(Check())
    #
    # @patch('misc.extract_settings_elvis')
    # @patch('misc.setting_logger')
    # @patch('check_elvis.Check.scamp')
    # def test_scamp_option_chosen(self, scamp, setting_logger,
    #                              extract_settings_elvis):
    #     """
    #
    #     :param scamp:
    #     :param setting_logger:
    #     :param extract_settings_elvis:
    #     :return:
    #     """
    #     scamp.return_value = True
    #     setting_logger.side_effect = MockedLogger
    #     extract_settings_elvis.return_value = True
    #
    #     sys.argv[1] = '-scamp'
    #
    #     return self.assertTrue(Check())
    #
    # @patch('misc.extract_settings_elvis')
    # @patch('misc.setting_logger')
    # @patch('check_elvis.Check.filt')
    # def test_filter_option_chosen(self, filt, setting_logger,
    #                               extract_settings_elvis):
    #     """
    #
    #     :param filt:
    #     :param setting_logger:
    #     :param extract_settings_elvis:
    #     :return:
    #     """
    #     filt.return_value = True
    #     setting_logger.side_effect = MockedLogger
    #     extract_settings_elvis.return_value = True
    #
    #     sys.argv[1] = '-filter'
    #
    #     return self.assertTrue(Check())
    #
    # @patch('misc.extract_settings_elvis')
    # @patch('misc.setting_logger')
    # @patch('check_elvis.Check.restart')
    # def test_scamp_option_chosen(self, restart, setting_logger,
    #                              extract_settings_elvis):
    #     """
    #
    #     :param restart:
    #     :param setting_logger:
    #     :param extract_settings_elvis:
    #     :return:
    #     """
    #     restart.return_value = True
    #     setting_logger.side_effect = MockedLogger
    #     extract_settings_elvis.return_value = True
    #
    #     sys.argv[1] = '-restart'
    #
    #     return self.assertTrue(Check())

    def tearDown(self):
        """

        :return:
        """
        pass


if __name__ == '__main__':
    main()
