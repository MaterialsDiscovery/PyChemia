"""
Several routines to read and write CIF files
"""

import os as _os


def cif_expand(path, dirname=None, verbose=False):
    """
    Split a multi-structure CIF file and return a directory
    with all the CIFs in separated files

    :param path: Path to the file that contain multi-structure

    :param dirname: Path to the directory where the extracted
        CIFs will be located, by default it will be the same
        location of the original file and the name of the directory
        will be the name of the 'filename' without extension

    :param verbose: Print some extra info

    :rtype : str
    """

    if not _os.path.isfile(path):
        raise ValueError('Could not find such filename')
    else:
        rf = open(path, 'r')

    if dirname is None:
        dirname = _os.path.splitext(path)[0]

    if not _os.path.lexists(dirname):
        _os.mkdir(dirname)
    elif not _os.path.isdir(dirname):
        dirname += '_DIR'
        if not _os.path.lexists(dirname):
            _os.mkdir(dirname)
        else:
            raise ValueError('Could not create a suitable directory')

    ndata = 0
    cif = ""
    cifname = None
    for line in rf.readlines():
        cif += line
        if line[:5] == 'data_':
            cifname = line.strip()[5:]
            if verbose:
                print(cifname)
            ndata += 1
        elif line[:13] == '#End of data_':
            if line.strip()[13:] != cifname:
                raise ValueError('Beginning and end of CIF are not consistent')
            wf = open(dirname + '/' + cifname + '.cif', 'w')
            wf.write(cif)
            wf.close()
            cif = ""
    rf.close()
    if verbose:
        print('Number of structures found: ', ndata)
    return dirname


def is_multistructure(path, verbose=False):
    """
    Read a filename in path and returns True if the CIF contains
    several structures
    """
    retval = False
    if not _os.path.isfile(path):
        raise ValueError('Could not find such filename')
    else:
        rf = open(path, 'r')

    ndata = 0
    for line in rf.readlines():
        if line[:5] == 'data_':
            ndata += 1
    if ndata > 1:
        retval = True

    if verbose:
        print('%4d structures in %s' % (ndata, path))

    return retval


def get_singlecifs(dirname, verbose=False):
    """
    Recursively explores a directory searching
    for .cif files, determines if they are single
    or multi-structure, expand those multi and
    return a list of single-strcuture cif files
    """

    single_cifs = []
    lst = [x for x in _os.listdir(dirname) if x[-3:] == 'cif']
    if verbose:
        print('Found '+str(len(lst))+' cifs')

    for cif in lst:
        path = dirname + '/' + cif
        if is_multistructure(path):
            if verbose:
                print(path+' is multistructure')
            cifdir = cif_expand(path)
            sublst = [x for x in _os.listdir(cifdir) if x[-3:] == 'cif']
            for subcif in sublst:
                subpath = cifdir + '/' + subcif
                single_cifs.append(subpath)
        else:
            single_cifs.append(path)
    return single_cifs
