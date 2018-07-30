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
from errors import ChangeTimeFailed
from misc import extract_settings_elvis
import times_elvis

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def change_times_mock_false():
    """

    :return:
    """
    return False


def parameters_mock():
    """

    :return:
    """
    return True


def clean_mock():
    """

    :return:
    """
    return False


def split_mock_true():
    """

    :return:
    """
    return True


def split_mock_false():
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


def restart_mock_true():
    """

    :return:
    """
    return True


def restart_mock_false():
    """

    :return:
    """
    return False


class TestFullPipelineUnsuccessfulSteps(TestCase):
    """

    """

    def setup(self):
        pass

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.restart', side_effect=restart_mock_false)
    def test_restart_fails(self, restart, parameters):
        """

        :param parameters:
        :param restart:
        :return:
        """
        argv[1] = '-full'

        return self.assertRaises(RestartFailed, Check)

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.restart', side_effect=restart_mock_true)
    @patch('check_elvis.Check.split', side_effect=split_mock_false)
    def test_split_fails(self, split, restart, parameters):
        """

        :param parameters:
        :param split:
        :return:
        """
        argv[1] = '-full'

        return self.assertRaises(SplitFailed, Check)

    @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    @patch('check_elvis.Check.restart', side_effect=restart_mock_true)
    @patch('check_elvis.Check.split', side_effect=split_mock_true)
    @patch('times_elvis.change_times', side_effect=change_times_mock_false)
    def test_change_times_fails(self, change_times, split, restart, parameters):
        """

        :param parameters:
        :param split:
        :return:
        """
        argv[1] = '-full'

        return self.assertRaises(ChangeTimeFailed, Check)

    # @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    # @patch('check_elvis.Check.clean', side_effect=clean_mock)
    # def test_clean_option_chosen(self, parameters, clean):
    #     """
    #
    #     :param parameters:
    #     :param clean:
    #     :return:
    #     """
    #     argv[1] = '-full'
    #
    #     return self.assertRaises(CleanFailed, Check)
    #
    # @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    # @patch('check_elvis.Check.split', side_effect=split_mock)
    # def test_split_option_chosen(self, parameters, split):
    #     """
    #
    #     :param parameters:
    #     :param split:
    #     :return:
    #     """
    #     argv[1] = '-split'
    #
    #     return self.assertRaises(SplitFailed, Check)
    #
    # @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    # @patch('check_elvis.Check.sextractor', side_effect=sextractor_mock)
    # def test_sextractor_option_chosen(self, parameters, sextractor):
    #     """
    #
    #     :param parameters:
    #     :param sextractor:
    #     :return:
    #     """
    #     argv[1] = '-sextractor'
    #
    #     return self.assertRaises(SextractorFailed, Check)
    #
    # @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    # @patch('check_elvis.Check.scamp', side_effect=scamp_mock)
    # def test_scamp_option_chosen(self, parameters, scamp):
    #     """
    #
    #     :param parameters:
    #     :param scamp:
    #     :return:
    #     """
    #     argv[1] = '-scamp'
    #
    #     return self.assertRaises(ScampFailed, Check)
    #
    # @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    # @patch('check_elvis.Check.filt', side_effect=filt_mock)
    # def test_filter_option_chosen(self, parameters, filt):
    #     """
    #
    #     :param parameters:
    #     :param filt:
    #     :return:
    #     """
    #     argv[1] = '-filter'
    #
    #     return self.assertRaises(FiltFailed, Check)
    #
    # @patch('check_elvis.Check.parameters', side_effect=parameters_mock)
    # @patch('check_elvis.Check.restart', side_effect=restart_mock)
    # def test_scamp_option_chosen(self, parameters, restart):
    #     """
    #
    #     :param parameters:
    #     :param restart:
    #     :return:
    #     """
    #     argv[1] = '-full'
    #
    #     return self.assertRaises(RestartFailed, Check)
    #
    # def teardrown(self):
    #     pass
