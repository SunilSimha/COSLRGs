
""" Make objects Galaxy, CGM, CGMAbsSurvey
"""

# import
import pdb
import os

import numpy as np
import glob as glob
import matplotlib.pyplot as plt

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
from cos_lrg.io import load_abssys, load_summ, write_file


try:
    basestring
except NameError:  # For Python 3
    basestring = str


chk_vel=False
chk_sep=False
chk_z=False


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


def make_cgmabs(icoord, towrite = True):
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

    igmsys, full_file = load_abssys(icoord,foldername='lrg_xabssys')
    cgmabs_1 = cgm.CGMAbsSys(gal, igmsys)

    # write in a file
    if towrite:
        write_file('cgmabs',cgmabs_1,iicoord = icoord, suff = 'cgmabs')

    return cgmabs_1



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
        name = 'J{:s}{:s}'.format(
            icoords.ra.to_string(unit=u.hour,
                sep='',pad=True, precision=2),icoords.dec.to_string(sep='',
                pad=True,alwayssign=True, precision=1))
        icoord = get_coord(name)
        icgm = make_cgmabs(icoord)
        lrgs_cgmabs.append(icgm)
    LRGsurvey.cgm_abs = lrgs_cgmabs # list of CGMAbsSys objects - append for all LRGs

    return LRGsurvey



def ew_figure(LRGsurvey,iline='HI 1215',summ_file=None):
    """
        Figure EW vs Rperp

    Parameters
    ----------
    LRGsurvey : survey object
    iline: str
      line
    summfile : str

    Returns
    -------
    measured EWs for a survey
    """
    #Rperp and MgII
    summ = load_summ(summ_file)
    Rperp = summ['RP_MPC']*1000.
    Mgiiind = summ['HAVE_MGII']
    Mgiicol = []
    for i in np.arange(len(Mgiiind)):
        if Mgiiind[i] == 1:
            Mgiicol.append('r')
        else:
            Mgiicol.append('b')

    # EWs
    abssystems = LRGsurvey.cgm_abs
    EW_iline = []
    EW_iline_sig = []
    for asys in abssystems:
        HIind = False
        for icomp in asys._components:
            if icomp.name[0:len(iline.split(' ')[0])] == iline.split(' ')[0]:
                for jcomp in icomp._abslines:
                    if jcomp.name[0:len(iline)] == iline:
                        EW_iline.append(jcomp.attrib['EW'].value)
                        EW_iline_sig.append(jcomp.attrib['sig_EW'].value)
                        HIind = True
        if HIind == False:
            EW_iline.append(0)
            EW_iline_sig.append(0)

    # figure
    for i in np.arange(len(Rperp)):
        plt.plot([Rperp[i],Rperp[i]],[EW_iline[i]-EW_iline_sig[i], EW_iline[i]+EW_iline_sig[i] ],color=Mgiicol[i])
        plt.scatter([Rperp[i]],[EW_iline[i]],c=Mgiicol[i])
    plt.xlabel('R [kpc]')
    plt.ylabel('EW('+iline+') ['+r'$\AA$'+']' )

    plt.show()

    return [Rperp, [EW_iline, EW_iline_sig]]





