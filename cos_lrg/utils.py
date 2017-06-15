""" Utilties for COS-LRG project
"""

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

from cos_lrg.io import load_summ

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
