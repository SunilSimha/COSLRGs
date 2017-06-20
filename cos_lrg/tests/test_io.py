# Module to run tests on I/O
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.table import Table

from linetools.spectra.xspectrum1d import XSpectrum1D
from pyigm.abssys.igmsys import IGMSystem

#from ..io import load_abssys
#from ..io import load_summ
#from ..io import load_spectrum
#from ..utils import get_coord

from cos_lrg.utils import get_coord
from cos_lrg.io import load_abssys
from cos_lrg.io import load_summ
from cos_lrg.io import load_spectrum


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_load_spectrum():
    spec = load_spectrum('J0226+0015')
    assert isinstance(spec, XSpectrum1D)


def test_load_abssys():
    """ Load AbsSystem 
    """
    abssys, _ = load_abssys('J0226+0015')
    # Test
    assert isinstance(abssys, IGMSystem)
    #
    coord = get_coord('J0226+0015')
    abssys, filename = load_abssys(coord)
    # Test
    assert isinstance(abssys, IGMSystem)
    assert isinstance(filename, str)

def test_load_summ():
    """ Test loading the summary file
    """
    summ = load_summ()
    # Test
    assert len(summ) == 15
    assert isinstance(summ, Table)
