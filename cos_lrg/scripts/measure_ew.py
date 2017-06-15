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

    from astropy.coordinates import SkyCoord

    from linetools.scripts.utils import coord_arg_to_coord
    from linetools.utils import radec_to_coord

    from cos_lrg.utils import match_coord_to_summ


    #
    if args.lrg_name == 'ALL':
        pdb.set_trace() # NOT READY FOR THIS YET
    else:
        ic = coord_arg_to_coord(args.lrg_name)
        icoords = [radec_to_coord(ic)]

    for icoord in icoords:
        # Find the closest
        #
        row = match_coord_to_summ(icoord)
        # Load the AbsSystem
        pdb.set_trace()





