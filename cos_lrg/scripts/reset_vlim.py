#!/usr/bin/env python

""" Measure EWs for one or more LRGs
"""

import pdb
from cos_lrg.io import load_abssys, load_summ, match_coord_to_summ
from astropy.coordinates import SkyCoord
import astropy.units as u


def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Reset vlim values of systems in JSON files v0.1')
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
    icoords = []
    names = []
    summ = load_summ()  # (summ_file=None)
    lrg_qso_coords = SkyCoord(ra=summ['RA_QSO'], dec=summ['DEC_QSO'], unit='deg')
    for iicoords in lrg_qso_coords:
        name = 'J{:s}{:s}'.format(
            iicoords.ra.to_string(unit=u.hour,
                                 sep='', pad=True, precision=2), iicoords.dec.to_string(sep='',
                                                                                       pad=True, alwayssign=True,
                                                                                       precision=1))
        names.append(name)
        icoords.append(get_coord(name))

    # Loop on the systems

    #for icoord,iicoord in icoords,lrg_qso_coords:
    for iname in names:
        # Load the AbsSystem
        abssys, filename = load_abssys(iname, chk_z=False)
        # Set coords to QSO coords (abssys.coord, components.coord, abslines.coord)
        ilrg_qso_coord = get_coord(iname)
        abssys.coord = ilrg_qso_coord
        for icomp in abssys._components:
            icomp.coord = ilrg_qso_coord
            for iline in icomp._abslines:
                iline.attrib['coord'] = ilrg_qso_coord
        print(abssys)
        # Set zem to QSO zem
        row = match_coord_to_summ(ilrg_qso_coord)
        zqso = row['Z_QSO']
        abssys.zem = zqso
        print("New zem = ", abssys.zem)
        # Reset vlim
        print("Original vlim = ", abssys.vlim)
        abssys.update_vlim()
        print("New vlim = ", abssys.vlim)
        # Write to file
        if not args.skip_check:
            print("About to overwrite the file: {:s}".format(filename))
            warnings.warn("Continue only if you know what you are doing!!")
            pdb.set_trace()
        abssys.write_json(outfil=filename)






