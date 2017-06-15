# Module to run tests on I/O
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.table import Table

from ..io import load_summ

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_load_summ():
    """ Test loading the summary file
    """
    summ = load_summ()
    # Test
    assert len(summ) == 15
    assert isinstance(summ, Table)
