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


def extract_settings_elvis_mock():
    """

    :return:
    """
    return True


def setting_logger_mock():
    """

    :return:
    """
    return True


def get_os_mock():
    """

    :return:
    """
    return True


def full_pipeline_mock():
    """

    :return:
    """
    return True


class TestCheckOptions(TestCase):
    """

    """

    def setup(self):
        pass

    @patch('check_elvis.Check.parameters')
    @patch('check_elvis.Check.full_pipeline')
    def test_fulloptionchoosen(self, parameters, full_pipeline):
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
    @patch('check_elvis.Check.split')
    def test_fulloptionchoosen(self, parameters, split):
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

    def teardrown(self):
        pass
