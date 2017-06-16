
import numpy as np
import os
import astropy.units as u
from linetools.spectralline import AbsLine, SpectralLine
from linetools.spectra.xspectrum1d import XSpectrum1D
from pyigm.abssys.igmsys import IGMSystem
import pdb


def abs_sys(allabslrgs):
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
