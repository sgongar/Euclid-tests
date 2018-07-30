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
from errors import FullPipelineFailed, CleanFailed, SplitFailed
from errors import SextractorFailed, ScampFailed, FiltFailed, RestartFailed
from misc import extract_settings_elvis


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def parameters_mock():
    """

    :return:
    """
    return True


def full_pipeline_mock():
    """

    :return:
    """
    return False


def clean_mock():
    """

    :return:
    """
    return False


def split_mock():
    """

    :return:
    """
    return False


def sextractor_mock():
    """

    :return:
    """
    return False


def scamp_mock():
    """

    :return:
    """
    return False


def filt_mock():
    """

    :return:
    """
    return False


def restart_mock():
    """

    :return:
    """
    return False


class TestCheckUnsuccessfulOptions(TestCase):
    """

    """

    def setup(self):
        pass

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.full_pipeline', side_effect=full_pipeline_mock)
    def test_full_option_chosen(self, full_pipeline, parameters):
        """

        :param parameters:
        :param full_pipeline:
        :return:
        """
        argv[1] = '-full'

        return self.assertRaises(FullPipelineFailed, Check)

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.clean', side_effect=clean_mock)
    def test_clean_option_chosen(self, clean, parameters):
        """

        :param parameters:
        :param clean:
        :return:
        """
        argv[1] = '-clean'

        return self.assertRaises(CleanFailed, Check)

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.split', side_effect=split_mock)
    def test_split_option_chosen(self, split, parameters):
        """

        :param parameters:
        :param split:
        :return:
        """
        argv[1] = '-split'

        return self.assertRaises(SplitFailed, Check)

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.sextractor', side_effect=sextractor_mock)
    def test_sextractor_option_chosen(self, sextractor, parameters):
        """

        :param parameters:
        :param sextractor:
        :return:
        """
        argv[1] = '-sextractor'

        return self.assertRaises(SextractorFailed, Check)

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.scamp', side_effect=scamp_mock)
    def test_scamp_option_chosen(self, scamp, parameters):
        """

        :param parameters:
        :param scamp:
        :return:
        """
        argv[1] = '-scamp'

        return self.assertRaises(ScampFailed, Check)

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.filt', side_effect=filt_mock)
    def test_filter_option_chosen(self, filt, parameters):
        """

        :param parameters:
        :param filt:
        :return:
        """
        argv[1] = '-filter'

        return self.assertRaises(FiltFailed, Check)

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.restart', side_effect=restart_mock)
    def test_scamp_option_chosen(self, restart, parameters):
        """

        :param parameters:
        :param restart:
        :return:
        """
        argv[1] = '-restart'

        return self.assertRaises(RestartFailed, Check)

    def teardrown(self):
        pass
