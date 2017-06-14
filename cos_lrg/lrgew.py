
# import
import pdb
import numpy as np
import glob
jj
from linetools.spectra import io as lsio
from pyigm.abssys.lls import LLSSystem


def measure_ew(xabsspath,contflspath, do_aodm=False, verbose=True):
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
    
    #for i in np.arange(len(xabssfiles)):
    #    print(xabssfiles[i],contfls[i])
    
    for i in np.arange(len(xabssfiles)):
        if verbose:
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
            if do_aodm:
                absline.measure_aodm()

            if verbose:
                print(absline)
                print('EW = {:g} with error {:g}'.format(absline.attrib['EW'],absline.attrib['sig_EW']))
                print(' ')
        if verbose:
            print('====================================')
            print('   ')
        # Write updated JSON file to disk


