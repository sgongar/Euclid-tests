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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from check_elvis import Check
from errors import FullPipelineFailed, CleanFailed, SplitFailed
from errors import SextractorFailed, ScampFailed, FiltFailed, RestartFailed
from errors import ChangeTimeFailed
import times_elvis

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


class TestFullPipelineUnsuccessfulSteps(TestCase):
    """

    """
    def setUp(self):
        """

        :return:
        """
        pass

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('check_elvis.Check.restart')
    def test_restart_fails(self, restart, setting_logger,
                           extract_settings_elvis):
        """

        :param restart:
        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        restart.return_value = False
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = True

        sys.argv[1] = '-full'

        return self.assertRaises(RestartFailed, Check)

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('check_elvis.Check.restart')
    @patch('check_elvis.Check.split')
    def test_split_fails(self, split, restart, setting_logger,
                         extract_settings_elvis):
        """

        :param split:
        :param restart:
        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        split.return_value = False
        restart.return_value = True
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = True

        sys.argv[1] = '-full'

        return self.assertRaises(SplitFailed, Check)

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('check_elvis.Check.restart')
    @patch('check_elvis.Check.split')
    @patch('times_elvis.change_times')
    def test_change_times_fails(self, change_times, split, restart,
                                setting_logger, extract_settings_elvis):
        """

        :param change_times:
        :param split:
        :param restart:
        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        change_times.return_value = False
        split.return_value = True
        restart.return_value = True
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = True

        sys.argv[1] = '-full'

        return self.assertRaises(ChangeTimeFailed, Check)

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('check_elvis.Check.restart')
    @patch('check_elvis.Check.split')
    @patch('times_elvis.change_times')
    @patch('check_elvis.Check.sextractor')
    def test_sextractor_fails(self, sextractor, change_times, split, restart,
                              setting_logger, extract_settings_elvis):
        """

        :param sextractor:
        :param change_times:
        :param split:
        :param restart:
        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        sextractor.return_value = False
        change_times.return_value = True
        split.return_value = True
        restart.return_value = True
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = True

        sys.argv[1] = '-sextractor'

        return self.assertRaises(SextractorFailed, Check)

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('check_elvis.Check.restart')
    @patch('check_elvis.Check.split')
    @patch('times_elvis.change_times')
    @patch('check_elvis.Check.sextractor')
    @patch('check_elvis.Check.scamp')
    def test_scamp_fails(self, scamp, sextractor, change_times, split,
                         restart, setting_logger, extract_settings_elvis):
        """

        :param scamp:
        :param sextractor:
        :param change_times:
        :param split:
        :param restart:
        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        scamp.return_value = False
        sextractor.return_value = True
        change_times.return_value = True
        split.return_value = True
        restart.return_value = True
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = True

        sys.argv[1] = '-scamp'

        return self.assertRaises(ScampFailed, Check)

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('check_elvis.Check.restart')
    @patch('check_elvis.Check.split')
    @patch('times_elvis.change_times')
    @patch('check_elvis.Check.sextractor')
    @patch('check_elvis.Check.scamp')
    @patch('check_elvis.Check.filt')
    def test_filter_fails(self, filt, scamp, sextractor, change_times, split,
                          restart, setting_logger, extract_settings_elvis):
        """

        :param filt:
        :param scamp:
        :param sextractor:
        :param change_times:
        :param split:
        :param restart:
        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        filt.return_value = False
        scamp.return_value = True
        sextractor.return_value = True
        change_times.return_value = True
        split.return_value = True
        restart.return_value = True
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = True

        sys.argv[1] = '-filter'

        return self.assertRaises(FiltFailed, Check)

    def tearDown(self):
        """

        :return:
        """
        pass

if __name__ == '__main__':
    main()
