"""
Several routines to read and write CIF files
"""

import os as _os
import re
from collections import OrderedDict

from pychemia.crystal import Lattice
from pychemia.utils.computing import read_file, only_ascii


class CIF:
    def __init__(self, data):
        self.data = data
        self.splits = OrderedDict()
        self.split_data()

    @staticmethod
    def from_string(string):
        return CIF(string)

    @staticmethod
    def from_file(filename):
        return CIF(read_file(filename))

    def split_data(self):
        splits = re.split('^\s*data_', 'CIF\n' + self.data, flags=re.MULTILINE)
        if len(splits) > 1:
            for x in splits[1:]:
                blk = self.block_from_string('data_' + x)
                self._process_block(blk)
                self.splits[blk['title']] = blk

    def block_from_string(self, string):

        ret = OrderedDict()
        string = self._process_string(string)
        n = 0
        in_loop = False
        for x in string.splitlines():
            x = x.strip()
            if x.startswith('data_'):
                ret['title'] = x[5:]
            elif x.startswith('_'):
                y = x.split()
                if in_loop:
                    ret['loop_' + str(n)][y[0]] = " ".join(z for z in y[1:])
                else:
                    ret[y[0]] = " ".join(z for z in y[1:])
            elif x.startswith('loop_'):
                in_loop = True
                n += 1
                ret['loop_' + str(n)] = OrderedDict()
            else:
                if in_loop:
                    y = x.split()
                    ret['loop_' + str(n)][y[0]] = " ".join(z for z in y[1:])
        return ret

    @classmethod
    def _process_block(cls, blk):

        for x in blk.keys():

            if x in ['_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha',
                     '_cell_angle_beta', '_cell_angle_gamma', '_cell_volume']:
                blk[x] = float(blk[x])

            elif x in ['_cell_formula_units_Z', '_symmetry_Int_Tables_number']:
                blk[x] = int(blk[x])

            elif x.startswith('loop'):
                if blk[x].keys()[0] == '_symmetry_equiv_pos_site_id':
                    for y in blk[x].keys()[2:]:
                        blk[x][int(y)] = blk[x][y][1:-1]
                        blk[x].pop(y)
                elif blk[x].keys()[0] == '_atom_site_type_symbol':
                    for y in blk[x].keys():
                        tokens = blk[x][y].split()
                        if len(tokens) == 6:
                            blk[x][y] = (tokens[0], int(tokens[1]), float(tokens[2]),
                                         float(tokens[3]), float(tokens[4]), int(tokens[5]))
        return blk

    @classmethod
    def _process_string(cls, string):
        # Remove comments
        string = re.sub("(\s|^)#.*$", "", string, flags=re.MULTILINE)
        # Remove empty lines
        string = re.sub("^\s*\n", "", string, flags=re.MULTILINE)
        # Remove non_ascii characters
        string = only_ascii(string)
        return string


class CIFStructure:
    def __init__(self, cif_blocks, title=None):

        if title is None:
            title = cif_blocks.keys()[0]
        self.data = cif_blocks[title]

    def get_lattice(self):
        """
        Generate the lattice from the provided lattice parameters. In
        the absence of all six lattice parameters, the crystal system
        and necessary parameters are parsed
        """
        length_strings = ("a", "b", "c")
        angle_strings = ("alpha", "beta", "gamma")

        try:
            lengths = [self.data["_cell_length_" + i] for i in length_strings]
            angles = [self.data["_cell_angle_" + i] for i in angle_strings]
            return Lattice.from_parameters_to_cell(*tuple(lengths + angles))
        except:
            raise ValueError()

    def get_symmetry_operations(self):

        ret = []
        for x in self.data.keys():
            if x.startswith('loop'):

                for symmetry_label in ["_symmetry_equiv_pos_as_xyz",
                                       "_symmetry_equiv_pos_as_xyz_",
                                       "_space_group_symop_operation_xyz",
                                       "_space_group_symop_operation_xyz_"]:
                    if symmetry_label in self.data[x].keys():
                        ret = [self.data[x][y] for y in self.data[x].keys() if self.data[x][y] != '']
                        break

        if not ret:
            print(self.data.keys())
            for symmetry_label in ["_space_group_IT_number",
                                   "_space_group_IT_number_",
                                   "_symmetry_Int_Tables_number",
                                   "_symmetry_Int_Tables_number_"]:

                if symmetry_label in self.data.keys():
                    spg = self.data[symmetry_label]
                    print(spg)
                    break
        if not ret:
            print(self.data.keys())
            for symmetry_label in ["_symmetry_space_group_name_H-M",
                                   "_symmetry_space_group_name_H_M",
                                   "_symmetry_space_group_name_H-M_",
                                   "_symmetry_space_group_name_H_M_",
                                   "_space_group_name_Hall",
                                   "_space_group_name_Hall_",
                                   "_space_group_name_H-M_alt",
                                   "_space_group_name_H-M_alt_",
                                   "_symmetry_space_group_name_hall",
                                   "_symmetry_space_group_name_hall_",
                                   "_symmetry_space_group_name_h-m",
                                   "_symmetry_space_group_name_h-m_"]:

                if symmetry_label in self.data.keys():
                    spg = self.data[symmetry_label]
                    print(spg)
                    break

        return ret


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

    :param dirname: (str) Directory to read CIFS
    :param verbose: (bool) Verbosity of routine
    :return:
    """

    single_cifs = []
    lst = [x for x in _os.listdir(dirname) if x[-3:] == 'cif']
    if verbose:
        print('Found ' + str(len(lst)) + ' cifs')

    for cif in lst:
        path = dirname + '/' + cif
        if is_multistructure(path):
            if verbose:
                print(path + ' is multistructure')
            cifdir = cif_expand(path)
            sublst = [x for x in _os.listdir(cifdir) if x[-3:] == 'cif']
            for subcif in sublst:
                subpath = cifdir + '/' + subcif
                single_cifs.append(subpath)
        else:
            single_cifs.append(path)
    return single_cifs


def get_space_group_symop_operation_xyz(filename):
    rf = open(filename)
    data = rf.read()
    rf.close()

    re.findall('_space_group_symop_operation_xyz\n([\s\w,]+) \n', data)
