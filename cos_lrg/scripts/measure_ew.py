#!/usr/bin/env python

""" Measure EWs for one or more LRGs
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Measure EWs v0.1')
    parser.add_argument("lrg_name", type=str, help="LRG name (e.g. J0226+0015 or ALL")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False, **kwargs):
    """ Run
    """
    import numpy as np

    from astropy.coordinates import SkyCoord


    from cos_lrg.utils import get_coord
    from cos_lrg.io import load_abssys


    #
    if args.lrg_name == 'ALL':
        pdb.set_trace() # NOT READY FOR THIS YET
    else:
        icoords = [get_coord(args.lrg_name)]

    for icoord in icoords:
        # Find the closest
        #
        row = match_coord_to_summ(icoord)
        coord = SkyCoord(ra=row['RA_QSO'], dec=row['DEC_QSO'], unit='deg')

        # Load the AbsSystem
        load_abssys(coord, zlrg=row['Z_GAL'])
        pdb.set_trace()





