
# import
import pdb
from astropy import constants as const
import numpy as np
from linetools.isgm.abscomponent import AbsComponent
from linetools import utils as ltu
from pyigm.igm.igmsightline import IGMSightline
import astropy.units as u
from linetools.spectralline import AbsLine, SpectralLine
from linetools import spectralline as ltsp
from linetools.spectra.xspectrum1d import XSpectrum1D
import json
from pyigm.abssys.lls import LLSSystem
import glob
from linetools.spectra import io as lsio



def measure_ew(xabsspath,contflspath):
    """ Printing EWs

    Parameters
    ----------
    xabsspath: str
        path to the abssys files
        The files should be in one folder which does not contain other files
    contflspath: str
        path to the continuum normalized spectra
        The files should be in one folder which does not contain other files
    
    Returns
    -------

    """

    xabssfiles=glob.glob(xabsspath+'*')
    contfls=glob.glob(contflspath+'*')
    #sort / match
    xabssfiles=np.sort(xabssfiles)
    contfls=np.sort(contfls) 
    
    #for i in np.arange(len(xabssfiles)):
    #    print(xabssfiles[i],contfls[i])
    
    for i in np.arange(len(xabssfiles)):
        print(xabssfiles[i])
        lls = LLSSystem.from_json(xabssfiles[i], chk_vel=False)
        spec_file = contfls[i]
        abs_lines = lls.list_of_abslines()

        for absline in abs_lines:
            #Spectrum
            spec = lsio.readspec(spec_file)
            absline.analy['spec'] = spec
            # EW
            absline.measure_restew(flg=1)  # Boxcar
            # AODM
            absline.measure_aodm()

            print(absline)
            print('EW = {:g} with error {:g}'.format(absline.attrib['EW'],absline.attrib['sig_EW']))
            print(' ')
        print('====================================')
        print('   ')


