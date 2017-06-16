""" Utilties for COS-LRG project
"""

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

from linetools.scripts.utils import coord_arg_to_coord
from linetools.utils import radec_to_coord


def get_coord(icoord, QSO=True):
    """ 
    Parameters
    ----------
    icoord : str or SkyCoord
      e.g. J0026+0015

    Returns
    -------
    coord : SkyCoord

    """
    # SkyCoord?
    if isinstance(icoord, SkyCoord):
        return icoord
    #
    if len(icoord) == 10:
        # Padding
        gd_coord = 'J{:s}00{:s}00'.format(icoord[1:5], icoord[-5:])
    else:
        gd_coord = icoord
    # Another manipulation
    ic = coord_arg_to_coord(gd_coord)
    # Finally a SkyCoord
    scoord = radec_to_coord(ic)
    # Now match!
    row = match_coord_to_summ(scoord)
    # Finally
    if QSO:
        coord = SkyCoord(ra=row['RA_QSO'], dec=row['DEC_QSO'], unit='deg')
    else:
        coord = SkyCoord(ra=row['RA_GAL'], dec=row['DEC_GAL'], unit='deg')
    # Return
    return coord


def match_coord_to_summ(coord, verbose=True):
    """ 
    Parameters
    ----------
    coord : SkyCoord

    Returns
    -------
    row : Row
      Row in the Summary file that is the closest match
      Could restrict to be within some tolerance (to avoid error)
    """
    from cos_lrg.io import load_summ # This needs to be imported here

    # Load Summary file
    summ = load_summ()
    lrg_qso_coords = SkyCoord(ra=summ['RA_QSO'], dec=summ['DEC_QSO'], unit='deg')

    # Match to closest
    idx = np.argmin(lrg_qso_coords.separation(coord))
    name = 'J{:s}{:s}'.format(
        lrg_qso_coords[idx].ra.to_string(unit=u.hour,
            sep='',pad=True, precision=2),lrg_qso_coords[idx].dec.to_string(sep='',
            pad=True,alwayssign=True, precision=1))
    if verbose:
        print("Associating your input LRG-QSO with {:s}".format(name))
    # Return
    return summ[idx]
