#!/usr/bin/env python

""" Measure EWs for one or more LRGs
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Measure EWs v0.1')
    parser.add_argument("lrg_name", type=str, help="LRG name (e.g. J0226+0015 or ALL")
    parser.add_argument("--skip_check", default=False, action='store_true', help="Skip set_trace that asks for permission")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False, **kwargs):
    """ Run
    """
    import numpy as np
    import warnings
    
    from cos_lrg.io import load_spectrum
    from cos_lrg.utils import get_coord
    from cos_lrg.io import load_abssys


    #
    if args.lrg_name == 'ALL':
        pdb.set_trace() # NOT READY FOR THIS YET
    else:
        icoords = [get_coord(args.lrg_name)]

    for icoord in icoords:
        # Load the Spectrum
        spec = load_spectrum(icoord)
        # Load the AbsSystem
        abssys, filename = load_abssys(icoord)#, zlrg=row['Z_GAL'])

        # Measure EW
        for iline in abssys.list_of_abslines():
            iline.analy['spec'] = spec
            iline.measure_restew(flg=1)  # Boxcar

        # Write to file
        if not args.skip_check:
            print("About to overwrite the file: {:s}".format(filename))
            warnings.warn("Continue only if you know what you are doing!!")
            pdb.set_trace()
        abssys.write_json(outfil=filename)






