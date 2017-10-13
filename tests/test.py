#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements


Todo:
    * Improve log messages

"""

from os import path

from unittest import TestCase, main
from mock import MagicMock, patch

path.append("/pcdisk/holly/sgongora/Documents/Euclid/Euclid-tests/")


__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
"""
__license__ = "GPL"
"""
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


class TestSimple(TestCase):

    """
    """

    def setUp(self):
        """
        define setUp parameters
        """
        pass

    def test(self):

        return True

    def tearDown(self):
        """
        define tearDown parameters
        """
        pass


if __name__ == '__main__':
    main()
