#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `NdexLoader` class."""

import os
import tempfile
import shutil

import unittest
from ndexstringloader.ndexloadstring import NDExSTRINGLoader


class FakeParam(object):
    pass


class TestNdexLoader(unittest.TestCase):
    """Tests for `NdexLoader` class."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_constructor(self):
        try:
            loader = NDExSTRINGLoader(None)
            self.fail('Expected AttributeError')
        except AttributeError as nt:
            pass
