def get_data(datafile,path_to_files):
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

