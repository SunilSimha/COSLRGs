#!/usr/bin/env python


import numpy as np
import os
import astropy.units as u
from linetools.spectralline import AbsLine, SpectralLine
from linetools.spectra.xspectrum1d import XSpectrum1D
from pyigm.abssys.igmsys import IGMSystem
from linetools.isgm import abscomponent as lt_abscomp
import pdb


def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Generate LRG system file from igmguesses v0.1')
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
    from cos_lrg.utils import get_coord, match_coord_to_summ
    from cos_lrg.io import load_abssys, load_guesses


    #
    if args.lrg_name == 'ALL':
        pdb.set_trace() # NOT READY FOR THIS YET
    else:
        icoords = [get_coord(args.lrg_name)]


    for icoord in icoords:
        row = match_coord_to_summ(icoord)
        zlrg = row['Z_GAL']
        # load guesses file and systems
        igm_systems, guesses_file = load_guesses(icoord)

        # find lrg system(s?)
        z_sys = [igm_sys.zabs for igm_sys in igm_systems]
        dzabssys = np.abs(np.array(z_sys) - zlrg)
        #dzabssys = []
        #for igm_sys in igm_systems:   ## igm_sys[0].zabs
        #    dzabssys.append(np.abs(igm_sys.zabs - zlrg))
        idx = np.argmin(dzabssys)
        lrgsys = igm_systems[idx]

        # write lrg system in a file

        # new folder
        lent = len(guesses_file.split('/')[-1]) + len(guesses_file.split('/')[-2]) + 1
        lrg_g_folder = guesses_file[0:-lent] + 'lrg_from_guesses/'
        pdb.set_trace()
        # Create new folder
        try:
            os.mkdir(lrg_g_folder)
        except OSError:  # likely already exists
            pass
        # filename
        # Build the filename
        coord = get_coord(icoord)
        ra = coord.ra.to_string(unit=u.hour, sep='', pad=True, precision=2)[0:4]
        dec = coord.dec.to_string(sep='', pad=True, alwayssign=True, precision=1)[0:5]
        # Full file
        lrg_g_file = 'LRG_guesses_J{:s}{:s}_z{:0.3f}.json'.format(ra, dec, zlrg)
        lrg_g_file_full = lrg_g_folder+lrg_g_file
        print("Writing/Overwriting file {:s} ?".format(lrg_g_file_full))
        pdb.set_trace()
        lrgsys.write_json(outfil=lrg_g_file_full)




