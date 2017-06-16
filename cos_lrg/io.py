""" Module for I/O with COS-LRGs
"""

import pdb
import os

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from pkg_resources import resource_filename

def get_data():
    """
    Getting files with data

    Parameters
    ----------
    datafile : str
        File with the list of objects (galaxies), and information about them
    path_to_files : str
        Path to the data

    Returns
    -------
    fileslist : list
        List with files names

    """

    allfiles = np.sort(glob.glob(path_to_files+'*'))
    fileslist = []

    # read fits table
    datainfo = Table.read(datafile)
    coords = ltu.radec_to_coord((datainfo['RA_QSO'], datainfo['DEC_QSO']))
    names = []
    for coord in coords:
        name = 'J{:s}{:s}'.format(coord.ra.to_string(unit=u.hour,sep='',pad=True, precision=2),coord.dec.to_string(sep='',pad=True,alwayssign=True, precision=1))
        names.append(name)
    datainfo['NAME'] = names

    for i in np.arange(len(datainfo)):
        name_beginning=datainfo['NAME'][i][0:5]
        ifiles=[]
        for ifile in allfiles:
            if ifile[0+len(path_to_files):5+len(path_to_files)] == name_beginning:
                ifiles.append(ifile)
        if len(ifiles) == 0:
            print('Warning: Corresponding file is missing or naming is different than used here. Skipping file.')
        elif len(ifiles) >= 2:
            print('Warning: Multiple files with similar names, or naming is different than used here. Appending all files')
        for ifile in ifiles:
            fileslist.append(ifile)

    return fileslist

def load_abssys(coord, zlrg=None):
    """ 
    Parameters
    ----------
    coord : SkyCoord
    zlrg : float, optional

    Returns
    -------
    abssys : GenericIGMSystem

    """
    from pyigm.abssys.igmsys import IGMSystem
    # Match coord to table to grab z
    if zlrg is None:
        summ = load_summ()
        lrg_qso_coords = SkyCoord(ra=summ['RA_QSO'], dec=summ['DEC_QSO'], unit='deg')
        idx = np.argmin(lrg_qso_coords.separation(coord))
    # Build the filename
    ra = coord.ra.to_string(unit=u.hour,sep='',pad=True, precision=2)[0:4]
    dec = coord.dec.to_string(sep='',pad=True,alwayssign=True, precision=1)[0:5]
    # Full file
    filename = 'LRG_J{:s}{:s}_z{:0.3f}.json'.format(ra, dec, zlrg)
    full_file = resource_filename('cos_lrg', 'data/lrg_abs_syst_1/{:s}'.format(filename))
    # Check
    if not os.path.isfile(full_file):
        print("No file named {:s} found!!".format(full_file))
        print("You may need to rename your file")
        pdb.set_trace()
    # Load
    abssys = IGMSystem.from_json(full_file)
    # Return
    return abssys



def load_summ(summ_file=None):
    """ Load the Summary file 
    Returns
    -------
    summ : Table
    """
    if summ_file is None:
        summ_file = resource_filename('cos_lrg', 'data/hstselect_final.fits')
    # Load
    summ = Table.read(summ_file)
    #
    return summ

