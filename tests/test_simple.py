import unittest

from os import getenv
from sys import argv, modules, path
from types import ModuleType

from unittest import TestCase, main
from mock import MagicMock, patch

home = getenv("HOME")
path.append('{}/build/sgongar/Euclid-tests/pipeline_elvis'.format(home))
path.append('{}/Dev/Euclid-tests/pipeline_elvis'.format(home))

from misc import fun

class MyTest(unittest.TestCase):
    def test(self):
        self.assertEqual(fun(3), 4)