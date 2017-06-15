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

from setup_package.py import get_package_data




#xabsspth = resource_filename('cos_lrg', 'data/lrg_xabssys')
#contflspth = resource_filename('cos_lrg', 'data/spectra')

xabsspth = 'cos_lrg/data/lrg_xabssys/'
contflspth = 'cos_lrg/data/spectra/'
outf = 'lines_ews.json'


def measure_ew(datafile,xabsspath=xabsspth,contflspath=contflspth,outfile=outf, do_aodm=False, verbose=True):
#def measure_ew(xabsspath,contflspath, do_aodm=False, verbose=True):
    """ Writing EWs to JSON file

    Parameters
    ----------
    datafile: str
        File with the list of objects (galaxies), and information about them
    xabsspath: str
        path to the abssys files
        The files should be in one folder which does not contain other files
    contflspath: str
        path to the continuum normalized spectra
        The files should be in one folder which does not contain other files
    outfile: str
        Output file with EWs
    do_aodm : bool, optional
    verbose : bool, optional
    
    Returns
    -------
    allews : dict
        EWs of all galaxies and selected transitions

    """
    
    xabssfiles=get_package_data(datafile,xabsspath)
    contfls=get_package_data(datafile,contflspath)
    
    allews = {}
    
    for i in np.arange(len(xabssfiles)):
        if verbose:
            print('Calculating EWs for ',xabssfiles[i])
 
        xabfile = xabssfiles[i] 
        contfile = contfls[i]  
        
        #name
        iname = xabfile[len(xabsspath):len(xabsspath)+5]
        allews[iname] = {}
        
        with open(xabfile) as json_data:
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
                    
                allews[iname][trans]={}
                allews[iname][trans]['ew']=iline.attrib['EW'] #.value
                allews[iname][trans]['sig_ew']=iline.attrib['sig_EW'] #.value
                
                # AODM
                if do_aodm:
                    iline.measure_aodm()
                 
                
                if verbose:
                    #print('EW = {:g} with error {:g}'.format(iline.attrib['EW'].value,iline.attrib['sig_EW'].value))
                    print('EW = {:g} with error {:g}'.format(iline.attrib['EW'],iline.attrib['sig_EW']))
                    print(' ')
        if verbose:
            print('====================================')
            print('   ')
            
    # Write updated JSON file to disk
    print('===================================')
    print('===================================')
        
    if verbose:
        print(allews)
        
    jdict = ltu.jsonify(allews)
    # Write
    ltu.savejson(outfile, jdict, easy_to_read=True, overwrite=True)   
        
    return allews
        