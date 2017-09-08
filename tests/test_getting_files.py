#!/usr/bin/python
# -*- coding: utf-8 -*-

from unittest import TestCase, main
from sys import path
import os
from os import listdir
from mock import MagicMock, patch


path.append("/pcdisk/holly/sgongora/Documents/Euclid/Euclid-tests/")

from errors import FolderNotPresent
from catalog_crossing import Catalog_crossing

"""
import os
import sys

# Dependencies for the tests
from mock import patch, Mock, MagicMock
from exceptions import KeyError
import exceptions
"""

def simple_urandom():
    return True

class TestFoldersAreNotPresent(TestCase):

    """
    Testing multiple client connections
    TDOD. Test multiple valid connections
    """

    def setUp(self):
        os.chdir('..')
        """
        os.listdir = MagicMock(return_value=['file1.txt', 'file2.txt',
                                          'file3.txt'])
        """
    @patch.object(Catalog_crossing.dirs, ['hello'])
    def test(self, test_function):
        # listdir('blabla')

        self.class_test = Catalog_crossing()

        return True

    def tearDown(self):
        print "out"


if __name__ == '__main__':
    main()
