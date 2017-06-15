
# import
import pdb
import numpy as np
import glob
import json

from linetools.spectra import io as lsio
from pyigm.abssys.lls import LLSSystem

import astropy.units as u
from linetools.spectralline import AbsLine, SpectralLine
from linetools import spectralline as ltsp
from linetools.spectra.xspectrum1d import XSpectrum1D

from pkg_resources import resource_filename

#xabsspth = resource_filename('cos_lrg', 'data/lrg_xabssys')
#contflspth = resource_filename('cos_lrg', 'data/spectra')

xabsspth = 'cos_lrg/data/lrg_xabssys/'
contflspth = 'cos_lrg/data/spectra/'


def measure_ew(xabsspath=xabsspth,contflspath=contflspth, do_aodm=False, verbose=True):
#def measure_ew(xabsspath,contflspath, do_aodm=False, verbose=True):
    """ Writing EWs to JSON file

    Parameters
    ----------
    xabsspath: str
        path to the abssys files
        The files should be in one folder which does not contain other files
    contflspath: str
        path to the continuum normalized spectra
        The files should be in one folder which does not contain other files
    do_aodm : bool, optional
    verbose : bool, optional
    
    Returns
    -------

    """

    xabssfiles=glob.glob(xabsspath+'*')
    contfls=glob.glob(contflspath+'*')
    #sort / match
    xabssfiles=np.sort(xabssfiles)
    contfls=np.sort(contfls) 
    
    
    
    for i in np.arange(len(xabssfiles)):
        if verbose:
            print('Calculating EWs for ',xabssfiles[i])
            
        xxfile = xabssfiles[i] 
        contfile = contfls[i]  
        with open(xxfile) as json_data:
            d = json.load(json_data)
        dd = d['components']
        if verbose:
            for ddkey in dd.keys():
                print(ddkey)

        for ddkey in dd.keys():
            ddlines=dd[ddkey]['lines']
            for line in ddlines.keys():
                if verbose:
                    print(ddkey,line)
                trans=ddlines[line]['analy']['name']
                iline = AbsLine(trans, z=ddlines[line]['limits']['z'])
                iline.attrib['z'] = ddlines[line]['limits']['z'] 
                iline.analy['vlim'] = ddlines[line]['limits']['vlim']['value']*u.km/u.s #format: [-250.,80.]*u.km/u.s
                # Set spectrum
                iline.analy['spec'] = XSpectrum1D.from_file(contfile)
                iline.limits.vlim[0]=ddlines[line]['limits']['vlim']['value'][0]*u.km/u.s
                iline.limits.wvlim[0]=ddlines[line]['limits']['wvlim']['value'][0]*u.AA
                iline.limits.zlim[0]=0. 
                #ddlines[line]['limits']['zlim']['value'][0]
                iline.limits.vlim[1]=ddlines[line]['limits']['vlim']['value'][1]*u.km/u.s
                iline.limits.wvlim[1]=ddlines[line]['limits']['wvlim']['value'][1]*u.AA
                iline.limits.zlim[1]=1. 
                #ddlines[line]['limits']['zlim']['value'][1]
                # Set analysis range
                iline.analy['wvlim'] = ddlines[line]['limits']['wvlim']['value']*u.AA   
                # EW
                ###iline.measure_ew() # Observer frame
                iline.measure_restew(flg=1)  # Boxcar
                # AODM
                if do_aodm:
                    iline.measure_aodm()

                if verbose:
                    print('EW = {:g} with error {:g}'.format(iline.attrib['EW'].value,iline.attrib['sig_EW'].value))
                    print(' ')
        if verbose:
            print('====================================')
            print('   ')
        # Write updated JSON file to disk

