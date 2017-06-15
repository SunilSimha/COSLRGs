#!/usr/bin/env python

""" Measure EWs for one or more LRGs
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Measure EWs v0.1')
    parser.add_argument("lrg_name", type=str, help="LRG name (e.g. J022600+001500 or ALL")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False, **kwargs):
    """ Run
    """
    import numpy as np
    from cos_lrg.io import load_summ

    from astropy.coordinates import SkyCoord
    from astropy import units as u

    from linetools.scripts.utils import coord_arg_to_coord
    from linetools.utils import radec_to_coord

    # Load Summary file
    summ = load_summ()
    lrg_qso_coords = SkyCoord(ra=summ['RA_QSO'], dec=summ['DEC_QSO'], unit='deg')

    #
    if args.lrg_name == 'ALL':
        pdb.set_trace() # NOT READY FOR THIS YET
    else:
        ic = coord_arg_to_coord(args.lrg_name)
        icoords = [radec_to_coord(ic)]

    for icoord in icoords:
        # Find the closest
        idx = np.argmin(lrg_qso_coords.separation(icoord))
        name = 'J{:s}{:s}'.format(lrg_qso_coords[idx].ra.to_string(unit=u.hour,
            sep='',pad=True, precision=2),lrg_qso_coords[idx].dec.to_string(sep='',
            pad=True,alwayssign=True, precision=1))
        #
        print("Associating your input LRG-QSO with {:s}".format(name))
        # Continue





