
""" Make objects Galaxy, CGM, CGMAbsSurvey
"""

# import
import pdb
import os

import numpy as np
import glob as glob
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from linetools.isgm.abscomponent import AbsComponent as lt_abscomp
from linetools import utils as ltu
from linetools.analysis import plots as ltaplots
from linetools.spectralline import AbsLine, SpectralLine
from linetools.spectra.xspectrum1d import XSpectrum1D
from pyigm.igm.igmsightline import IGMSightline
from pyigm.field import galaxy as galaxy #.py:class Galaxy(object)
from pyigm.cgm import cgm as cgm #.py:class CGM(object)
from pyigm.cgm.cgmsurvey import CGMAbsSurvey as cgmsurvey
from pyigm.abssys.igmsys import IGMSystem

from pyigm.cgm import cos_halos as pch

from pkg_resources import resource_filename

from cos_lrg.utils import match_coord_to_summ
from cos_lrg.utils import get_coord
from cos_lrg.io import load_abssys, load_summ, write_cgmabs_file, load_guesses


try:
    basestring
except NameError:  # For Python 3
    basestring = str


chk_vel=False
chk_sep=False
chk_z=False


def make_gal(icoord):
    """
        Make Galaxy object

    Parameters
    ----------
    icoord : SkyCoord or str
      e.g. J0026+0015

    Returns
    -------
    gal : Galaxy object

    """

    irow = match_coord_to_summ(icoord)
    zlrg = irow['Z_GAL']
    ra = irow['RA_GAL']*u.deg
    dec = irow['DEC_GAL']*u.deg
    radec = (ra, dec)
    gal = galaxy.Galaxy(radec,z=zlrg)
    return gal


def make_cgmabs(icoord, towrite = False, fromfile = True, filename = None):
    """
        Make or read from file a CGMAbsSys object

    Parameters
    ----------
    icoord : SkyCoord or str
      e.g. J0026+0015
    towrite : bool
      write in files
    fromfile : bool
      read from file
    filename : str
      read CGM object from file

    Returns
    -------
    cgmabs_1 : CGMAbsSys object

    """
    icoord = get_coord(icoord)
    irow = match_coord_to_summ(icoord)
    zlrg = irow['Z_GAL']
    R_lrg = irow['RP_MPC']*1000. #in kpc
    ra = irow['RA_GAL']*u.deg
    dec = irow['DEC_GAL']*u.deg
    radec = (ra, dec)

    gal = make_gal(icoord)

    if fromfile == False:
        igmsys, full_file = load_abssys(icoord,foldername='lrg_xabssys', chk_z=False)
        cgmabs_1 = cgm.CGMAbsSys(gal, igmsys)
    else:
        #print('Not yet ready')
        #pass
        ### not correct:   cgmabs_1 = cgm.CGMAbsSys.from_dict(filepath)
        if filename == None:
            datafolder = 'cgmabs'
            suff = 'cgmabs'
            folderpath = resource_filename('cos_lrg', 'data/' + datafolder + '/')
            coord = get_coord(icoord)
            ra = coord.ra.to_string(unit=u.hour, sep='', pad=True, precision=2)[0:4]
            dec = coord.dec.to_string(sep='', pad=True, alwayssign=True, precision=1)[0:5]
            filename = folderpath + suff + '_J{:s}{:s}.json'.format(ra, dec)  # _z{:0.3f}
            # read dict from a file
        #print("Step1")
        cgmasdict = ltu.loadjson(filename)
        #print("Step2")
        # dict to CGMAbsSys
        #igmsys = cgmasdict
        #cgmasdict2 = {'galaxy': gal, 'igm_sys' : cgmasdict}
        #pdb.set_trace()
        #cgmabs_1 = {'galaxy': gal, 'igm_sys' : cgmasdict}
        cgmabs_1 = cgm.CGMAbsSys.from_dict(cgmasdict) #, chk_z=False)  ## ,chk_vel=False  ### here
        #print("Step3")

    # write in a file
    if towrite:
        write_cgmabs_file(cgmabs_1,datafolder='cgmabs',iicoord = icoord, suff = 'cgmabs')

    return cgmabs_1



def make_survey(summ_file=None, towrite = True, fromfile = False, filename = None):
    """
        Make Survey object
        option 1: Read from the CGMAbsSys JSON files
        option 2: First make CGMAbsSys JSON files

    Parameters
    ----------
    summ_file : str
      summary file for the survey
    towrite : bool
      write in files CGMAbsSys (see make_cgmabs)
    fromfile : bool
      read from file CGMAbsSys (see make_cgmabs)
    filename : str
      read CGM object from file CGMAbsSys (see make_cgmabs)

    Returns
    -------
    lrgsurvey : Survey object

    """
    LRGsurvey = cgmsurvey()
    LRGsurvey.survey = 'LRG-QSO survey'
    lrgs_cgmabs = []

    summ = load_summ()
    #lrg_qso_coords = SkyCoord(ra=summ['RA_QSO'], dec=summ['DEC_QSO'], unit='deg')
    lrg_qso_coords = SkyCoord(ra=summ['RA_GAL'], dec=summ['DEC_GAL'], unit='deg')
    for icoords in lrg_qso_coords:
        name = 'J{:s}{:s}'.format(
            icoords.ra.to_string(unit=u.hour,
                sep='',pad=True, precision=2),icoords.dec.to_string(sep='',
                pad=True,alwayssign=True, precision=1))
        #icoord = get_coord(name)
        icgm = make_cgmabs(name, towrite = towrite, fromfile = fromfile, filename = filename) # icoord
        lrgs_cgmabs.append(icgm)
    LRGsurvey.cgm_abs = lrgs_cgmabs # list of CGMAbsSys objects - append for all LRGs

    return LRGsurvey



def ew_figure(LRGsurvey,iline='HI 1215',summ_file=None):
    """
        Figure EW vs Rperp

    Parameters
    ----------
    LRGsurvey : survey object
    iline: str
      line
    summfile : str

    Returns
    -------
    measured EWs for a survey
    """
    #Rperp and MgII
    summ = load_summ(summ_file)
    Rperp = summ['RP_MPC']*1000.
    Mgiiind = summ['HAVE_MGII']
    Mgiicol = []
    for i in np.arange(len(Mgiiind)):
        if Mgiiind[i] == 1:
            Mgiicol.append('r')
        else:
            Mgiicol.append('b')

    tbl = LRGsurvey.trans_tbl(iline)
    EW_iline = tbl['EW']
    EW_iline_sig = tbl['sig_EW']

    # figure
    #for i in np.arange(len(Rperp)):
    #plt.plot([Rperp[i],Rperp[i]],[EW_iline[i]-EW_iline_sig[i], EW_iline[i]+EW_iline_sig[i] ],color=Mgiicol[i])
    plt.scatter(Rperp,EW_iline,c=Mgiicol)
    j = np.where(Mgiiind == 1)
    plt.errorbar(Rperp, EW_iline, yerr=EW_iline_sig, ecolor='b', linestyle="None")
    plt.errorbar(Rperp[j], EW_iline[j], yerr=EW_iline_sig[j], ecolor='r', linestyle="None")

    plt.xlabel('R [kpc]')
    plt.ylabel('EW('+iline+') ['+r'$\AA$'+']' )

    plt.show()

    return [Rperp, [EW_iline, EW_iline_sig]]





def stack_plots(ymnx=(-0.35, 1.5), return_fig=True, vlim = None, add_ew=True,nrow=None,outpfolder=None,tight_layout=False):
    """
        Figure EW vs Rperp

    Parameters
    ----------
    LRGsurvey : survey object
    iline: str
      line
    summfile : str
    outpfolder : str
      output folder in which figures are saved.
      None - figures will not be saved

    Returns
    -------
    measured EWs for a survey
    """

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
        #icoords.append(get_coord(name))
        names.append(name)

    #for icoord in icoords:
    for iname in names:
        icoord = get_coord(iname)
        cgmabss = make_cgmabs(icoord, towrite=False, fromfile=True, filename=None)

        ##load guesses and read reliability of components, if = a or b (or c) : add component
        #igm_sys, guesses_file = load_guesses(iname)
        # ...


        abslines = []
        ews = []
        sigews = []

        for icomp in cgmabss.igm_sys._components:
            for iline in icomp._abslines:
                if iline.attrib['EW'] > abs(iline.attrib['sig_EW']):
                    if iline.attrib['EW'] > 0:
                        if iline.analy['do_analysis'] == 1:
                            if iline.analy['flg_eye'] == 0:
                                iline.analy['spec'] = XSpectrum1D.from_file(iline.analy['spec_file'])
                                abslines.append(iline)
                                ews.append(iline.attrib['EW'])
                                sigews.append(iline.attrib['sig_EW'])
                                #print(iline.attrib['EW'], iline.attrib['sig_EW'])

        if vlim == None:
            ivlim = cgmabss.vlim    #  * u.km / u.s
            print('Setting vlim = ', ivlim)
        if nrow == None:
            inrow = round((len(abslines)+1)/2)
        else:
            inrow = nrow

        fig = ltaplots.stack_plot(abslines, vlim=ivlim, ymnx= ymnx, figsz=(14,11), return_fig=return_fig, add_ew=add_ew, nrow=inrow,tight_layout=tight_layout)
        if outpfolder is not None:
            fig.savefig(outpfolder+iname+'.jpg')



'''
def ew_figure(LRGsurvey,iline='HI 1215',summ_file=None):
    """
        Figure EW vs Rperp

    Parameters
    ----------
    LRGsurvey : survey object
    iline: str
      line
    summfile : str

    Returns
    -------
    measured EWs for a survey
    """
    #Rperp and MgII
    summ = load_summ(summ_file)
    Rperp = summ['RP_MPC']*1000.
    Mgiiind = summ['HAVE_MGII']
    Mgiicol = []
    for i in np.arange(len(Mgiiind)):
        if Mgiiind[i] == 1:
            Mgiicol.append('r')
        else:
            Mgiicol.append('b')

    # EWs
    abssystems = LRGsurvey.cgm_abs
    EW_iline = []
    EW_iline_sig = []
    for asys in abssystems:
        HIind = False
        for icomp in asys._components:
            if icomp.name[0:len(iline.split(' ')[0])] == iline.split(' ')[0]:
                for jcomp in icomp._abslines:
                    if jcomp.name[0:len(iline)] == iline:
                        EW_iline.append(jcomp.attrib['EW'].value)
                        EW_iline_sig.append(jcomp.attrib['sig_EW'].value)
                        HIind = True
        if HIind == False:
            EW_iline.append(0)
            EW_iline_sig.append(0)

    # figure
    for i in np.arange(len(Rperp)):
        plt.plot([Rperp[i],Rperp[i]],[EW_iline[i]-EW_iline_sig[i], EW_iline[i]+EW_iline_sig[i] ],color=Mgiicol[i])
        plt.scatter([Rperp[i]],[EW_iline[i]],c=Mgiicol[i])
    plt.xlabel('R [kpc]')
    plt.ylabel('EW('+iline+') ['+r'$\AA$'+']' )

    plt.show()

    return [Rperp, [EW_iline, EW_iline_sig]]

'''
