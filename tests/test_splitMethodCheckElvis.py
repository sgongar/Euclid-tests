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


class MockedProcess:

    def __init__(self, target, args):
        """

        :param target:
        :param args:
        """
        pass

    def start(self):
        """

        :return:
        """
        return 'start'

    def is_alive(self):
        """

        :return:
        """
        return 'is_alive'


class TestSplitMethodFromCheckElvis(TestCase):
    """

    """

    def setup(self):
        """

        :return:
        """
        pass

    @patch('misc.extract_settings_elvis')
    @patch('misc.setting_logger')
    @patch('misc.get_fpa_elvis')
    @patch('multiprocessing.Process')
    def test_processes_created(self, Process, get_fpa_elvis, setting_logger,
                               extract_settings_elvis):
        """

        :param setting_logger:
        :param extract_settings_elvis:
        :return:
        """
        Process.side_effect = MockedProcess
        get_fpa_elvis.return_value = ['fpa_1', 'fpa_2', 'fpa_3', 'fpa_4']
        setting_logger.side_effect = MockedLogger
        extract_settings_elvis.return_value = {'fits_dir': 'fits_dir_mock',
                                               'fpas_dir': 'fpas_dir_mock'}

        argv[1] = '-split'

        self.assertTrue(Check)

    def teardrown(self):
        """

        :return:
        """
        pass
