#!/usr/bin/env python


import numpy as np
import os
import astropy.units as u
from linetools.spectralline import AbsLine, SpectralLine
from linetools.spectra.xspectrum1d import XSpectrum1D
from pyigm.abssys.igmsys import IGMSystem
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


### could be ignored?
def abs_sys_ignore(allabslrgs):
    """ 
    Generate systems for all LRGs
    
    Parameters
    ----------
    allabslrgs - dict
        Contains info about all absorption lines in all LRGs
        
    Returns
    -------
    linesystems - dict
        Dict with systems
    
    """
    
    lrgs = allabslrgs.keys()
    linesystems = {}
    
    for ilrg in lrgs:
        ilrgallabs = allabslrgs[ilrg]
        linekeys = ilrgallabs.keys()
        sysnames = []
        linesys = {}  
        linesystems[ilrg]={}

        for iline in linekeys:
            sysname = iline.split(" ")[0]
            if sysname not in sysnames:
                sysnames.append(sysname)
                linesys[sysname]=[]
            linesys[sysname].append(ilrgallabs[iline])

        for sysname in sysnames: 
            abscomp = AbsComponent.from_abslines(linesys[sysname])
            lsys = IGMSystem.from_components([abscomp])
            linesystems[ilrg][sysname] = lsys

    return linesystems




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

    linesystems = {}

    for icoord in icoords:
        # Load the Spectrum
        spec = load_spectrum(icoord)
        # Load the AbsSystem
        abssys, filename = load_abssys(icoord)#, zlrg=row['Z_GAL'])

	#sysnames = []
        linesys = {}  
        linesystems[ilrg]={}

        for iline in abssys.list_of_abslines():
            iline.analy['spec'] = spec
            #iline.measure_restew(flg=1)  # Boxcar

            sysname = iline.analy['name'].split(" ")[0]
            #if sysname not in sysnames:
	    if sysname not in linesys.keys():
                #sysnames.append(sysname)
                linesys[sysname]=[]
            linesys[sysname].append(iline)


        #for sysname in sysnames: 
        for sysname in linesys.keys(): 
            abscomp = AbsComponent.from_abslines(linesys[sysname])
            lsys = IGMSystem.from_components([abscomp])
            linesystems[ilrg][sysname] = lsys

        #write to file
        #abssys.write_json(outfil=filename_sys) 


    return linesystems






        # Write to file
        #if not args.skip_check:
        #    print("About to overwrite the file: {:s}".format(filename))
        #    warnings.warn("Continue only if you know what you are doing!!")
        #    pdb.set_trace()
        #abssys.write_json(outfil=filename)



