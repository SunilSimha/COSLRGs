
""" Make objects Galaxy, CGM, CGMAbsSurvey
"""

# import
import pdb
import os

import numpy as np
import glob as glob

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from linetools.isgm.abscomponent import AbsComponent as lt_abscomp
from linetools import utils as ltu
from linetools.spectralline import AbsLine, SpectralLine
from linetools.spectra.xspectrum1d import XSpectrum1D
from pyigm.igm.igmsightline import IGMSightline
from pyigm.field import galaxy as galaxy #.py:class Galaxy(object)
from pyigm.cgm import cgm as cgm #.py:class CGM(object)
from pyigm.cgm.cgmsurvey import CGMAbsSurvey as cgmsurvey
from pyigm.abssys.igmsys import IGMSystem

from pkg_resources import resource_filename

from cos_lrg.utils import match_coord_to_summ
from cos_lrg.utils import get_coord
from cos_lrg.io import load_abssys, load_summ


try:
    basestring
except NameError:  # For Python 3
    basestring = str



def make_gal(icoord):
    """
        Make Galaxy object

    Parameters
    ----------
    icoord : SkyCoord or str
      e.g. J0026+0015

    Returns
    -------
    gal : Galaxy object

    """

    irow = match_coord_to_summ(icoord)
    zlrg = irow['Z_GAL']
    ra = irow['RA_GAL']*u.deg
    dec = irow['DEC_GAL']*u.deg
    radec = (ra, dec)
    gal = galaxy.Galaxy(radec,z=zlrg)
    return gal


def make_cgm(icoord):
    """
        Make CGM object

    Parameters
    ----------
    icoord : SkyCoord or str
      e.g. J0026+0015

    Returns
    -------
    cgm1 : CGM object

    """
    irow = match_coord_to_summ(icoord)
    zlrg = irow['Z_GAL']
    R_lrg = irow['RP_MPC']*1000. #in kpc
    ra = irow['RA_GAL']*u.deg
    dec = irow['DEC_GAL']*u.deg
    radec = (ra, dec)
    gal = make_gal(icoord)
    cgm1 = cgm.CGM(gal)

    #radec_qso = (irow['RA_QSO']*u.deg, irow['DEC_QSO']*u.deg)   # e.g.
    #igmsys = IGMSystem('CGM', radec_qso, gal.z,
    #                   [-500, 500] * u.km / u.s)  # e.g. / or use existing one? #read from the file?
    ## should read from file make_system's file
    #igmsys, full_file = load_abssys(icoord,foldername='lrg_from_guesses')
    igmsys, full_file = load_abssys(icoord,foldername='lrg_xabssys')
    cgmabs_1 = cgm.CGMAbsSys(gal, igmsys)
    cgm1.rlim = R_lrg
    cgm1.cgm_abs = cgmabs_1  # = List of CGMAbsSys classes

    return cgm1



def make_survey(summ_file=None):
    """
        Make Survey object

    Parameters
    ----------
    summ_file : str
      summary file for the survey

    Returns
    -------
    lrgsurvey : Survey object

    """
    LRGsurvey = cgmsurvey()
    LRGsurvey.survey = 'LRG-QSO survey'
    lrgs_cgmabs = []

    summ = load_summ()
    #lrg_qso_coords = SkyCoord(ra=summ['RA_QSO'], dec=summ['DEC_QSO'], unit='deg')
    lrg_qso_coords = SkyCoord(ra=summ['RA_GAL'], dec=summ['DEC_GAL'], unit='deg')
    for icoords in lrg_qso_coords:
        #pdb.set_trace()
        name = 'J{:s}{:s}'.format(
            icoords.ra.to_string(unit=u.hour,
                sep='',pad=True, precision=2),icoords.dec.to_string(sep='',
                pad=True,alwayssign=True, precision=1))
        #pdb.set_trace()
        icoord = get_coord(name)
        icgm = make_cgm(icoord)
        lrgs_cgmabs.append(icgm.cgm_abs)
    LRGsurvey.cgm_abs = lrgs_cgmabs # list of CGMAbsSys objects - append for all LRGs

    # and write it in a file

    return LRGsurvey


#def measure_ew_from_survey():









#def ew_figure():

















