# Module to run tests on utils
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.table import Table, Row
from astropy.coordinates import SkyCoord

from ..utils import match_coord_to_summ
from ..utils import get_coord

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_coord():
    coord = get_coord('J0026+0015')
    #
    assert isinstance(coord, SkyCoord)

def test_match():
    """ Test loading the summary file
    """
    coord = get_coord('J0026+0015')
    row = match_coord_to_summ(coord)

    # Test
    assert 'RA_QSO' in row.colnames
    assert isinstance(row, Row)
