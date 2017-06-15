#!/usr/bin/env python

""" Measure EWs for one or more LRGs
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Measure EWs v0.1')
    parser.add_argument("lrg_name", type=str, help="LRG name or ALL")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False, **kwargs):
    """ Run
    """
    import numpy as np
    from astropy import units as u
    from specdb.utils import load_db
    from linetools.scripts.utils import coord_arg_to_coord

    # Load Summary file

    #
    if args.lrg_name == 'ALL':
        pdb.set_trace() # NOT READY FOR THIS YET
    else:


