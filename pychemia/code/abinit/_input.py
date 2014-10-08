#!/usr/bin/env python

"""
Definition of the class input to read
ABINIT input files and store their information
as a python dictionary called 'variables'
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2012"
__version__ = "1.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "guillermo.avendano@uclouvain.be"
__status__ = "Development"
__date__ = "Aug 27, 2012"

import os as _os
from scipy.io import netcdf_file as _netcdf_file
import numpy as _np

from pychemia.utils.periodic import atomic_symbol, covalent_radius, atomic_number
from pychemia.utils.constants import bohr_angstrom, angstrom_bohr
import pychemia.code.abinit
from pychemia.core import Structure
from pychemia.utils.mathematics import unit_vectors, length_vectors


class InputVariables:
    """
    An input object contains:

    data:
          variables = Dictionary whose keys are ABINIT variable names
                      and contains the values as numpy arrays

    methods:
            write = Write the input into as a text file that ABINIT
                    can use as an input file

            get_value = Get the value of a particular variable
            set_value = Set the value of a particular variable

            get_atomic_structure =
    """

    variables = {}

    def __init__(self, *args, **kwargs):
        """
        Creates a new input object, the input object
        contains a dictionary whose keys are the abinit
        input variables and the values are always numpy
        arrays.
        """
        filename = ''
        if len(args) == 1:
            x = args[0]
            if isinstance(x, pychemia.code.abinit.AbiFiles):
                filename = x.basedir + "/" + x.files["in"]
            elif _os.path.isfile(x):
                filename = x

        if 'filename' in kwargs:
            filename = kwargs['filename']

        if _os.path.isfile(filename):
            if filename[-3:] == '.in':
                self.__import_input(filename)
            elif filename[-6:] == '.files':
                abifile = pychemia.code.abinit.AbiFiles(filename)
                filename = abifile.get_input_filename()
                self.__import_input(filename)
            elif filename[-6:] == 'OUT.nc':
                self.variables = netcdf2dict(filename)
            else:
                try:
                    self.__import_input(filename)
                except ValueError:
                    print('File format not identified')

    def __import_input(self, filename):
        """
        Read an ABINIT input file and return a python dictionary
        with the input variables readed from that. The keys are
        the fullname input variables (acell,xcart21,etc). The
        values are numpy arrays of numbers except for
        the value '*[NUMBER]' that is keeped as string, and
        the string associated to the variable xyzfile

        Args:
            filename:
                ABINIT input filename
        """
        ans = pychemia.code.abinit.parser(filename)
        if ans is not None:
            self.variables = ans

    def write(self, filename):
        """
        Write an input object into a text
        file that ABINIT can use as an input
        file

        Args:
            filename:
                The 'abinit.in' filename that will be written
        """
        wf = open(filename, 'w')
        wf.write(self.__str__())
        wf.close()

    def __str__(self):
        """
        String representation of the object
        """
        ret = ''
        thekeys = self.variables.keys()
        varnames = [x for x in thekeys if not x[-1].isdigit()]

        if 'ndtset' in varnames:
            ret = ret + "#" + 60 * "-" + "\n#" + " MULTI DATASET\n#" + 60 * "-" + "\n\n"
            ret += self.write_key('ndtset')
            varnames.remove('ndtset')
            if 'jdtset' in varnames:
                ret += self.write_key('jdtset')
                varnames.remove('jdtset')
            if 'udtset' in varnames:
                ret += self.write_key('udtset')
                varnames.remove('udtset')
            ret += '\n'

        seqvarnames = [x for x in varnames if
                       (x[-1] == ':' or x[-1] == "+" or x[-1] == "?" or x[-2] == ':' or x[-2] == "+" or x[-2] == "?")]
        if len(seqvarnames) > 0:
            ret = ret + "#" + 60 * "-" + "\n#" + " SEQUENCE\n#" + 60 * "-" + "\n\n"
            seqvarnames.sort()
            for i in seqvarnames:
                ret += self.write_key(i)
                varnames.remove(i)
            ret += '\n'

        if len(varnames) > 0:
            varnames.sort()
            ret = ret + "#" + 60 * "-" + "\n#" + " ALL DATASETS\n#" + 60 * "-" + "\n\n"
            for i in varnames:
                ret += self.write_key(i)
            ret += '\n'

        for dtset in range(1, 100):
            varnames = [x for x in thekeys if
                        (x[-len(str(dtset)):] == str(dtset) and not x[-len(str(dtset)) - 1:].isdigit())]
            if len(varnames) > 0:
                varnames.sort()
                ret = ret + "#" + 60 * "-" + "\n#" + " DATASET " + str(dtset) + "\n#" + 60 * "-" + "\n\n"
                for i in varnames:
                    ret += self.write_key(i)
                ret += '\n'
        return ret

    def write_key(self, varname):
        """
        Receives an input variable and write their contents
        properly according with their kind and length

        Args:
            varname:
                The name of the input variable
            wf:
                The file object where the 'abinit.in' is been written
        """
        ret = ''
        if len(self.variables[varname]) == 0:
            print("[ERROR] input variable: '%s' contains no elements" % varname)
            return

        # Assume that the variables are integer and test if such assumption
        # is true
        integer = True
        real = False
        string = False
        compact = True

        # Get the general kind of values for the input variable
        for j in self.variables[varname]:

            try:
                if not float(j).is_integer():
                    # This is the case of non integer values
                    integer = False
                    real = True
                    string = False
                    if len(str(float(j))) > 7:
                        compact = False

            except ValueError:
                # This is the case of '*1' that could not
                # be converted because we don't know the size
                # of the array
                integer = False
                real = False
                string = True

        ret = ret + (varname.rjust(15)) + "  "

        for j in range(len(self.variables[varname])):

            if real:
                if compact:
                    ret += ("%g" % self.variables[varname][j]).rjust(8)
                else:
                    ret += "%17.10e" % self.variables[varname][j]
            elif integer:
                ret += "%d" % self.variables[varname][j]
            elif string:
                ret += "%s" % self.variables[varname][j]

            # Conditions to jump to a new line
            if ((j + 1) % 3) == 0 and real and j < len(self.variables[varname]) - 1:
                ret += "\n"
                ret += 17 * " "
            elif j < len(self.variables[varname]) - 1:
                ret += " "
        ret += "\n"
        return ret

    def get_value(self, varname, idtset=None, full=False):
        """
        Get the value of the input variable 'varname'
        associated with the dataset 'idtset'
        If 'idtset' is not given will asume that the
        value is not dataset dependent
        If the value is a single number the return
        will be a single number instead of a numpy array,
        or set full=True to get a numpy array always
        """
        name = ''
        fact = 1
        delta = 0

        # Get the right key for the abinit variable
        if idtset is None:
            if varname in self.variables:
                name = varname
        else:
            if (varname + str(idtset)) in self.variables:
                name = varname + str(idtset)
            elif idtset > 10:
                if (varname + '?' + (str(idtset)[1])) in self.variables:
                    name = varname + '?' + (str(idtset)[1])
                elif (varname + (str(idtset)[0]) + '?') in self.variables:
                    name = varname + (str(idtset)[0]) + '?'
                elif (varname + '+?') in self.variables and (varname + ':?') in self.variables:
                    name = varname + ':?'
                    fact = int(str(idtset)[0]) - 1
                    delta = self.variables[varname + '+?']
                elif (varname + '?+') in self.variables and (varname + '?:') in self.variables:
                    name = varname + '?:'
                    fact = int(str(idtset)[1]) - 1
                    delta = self.variables[varname + '?+']
            if name == '' and varname in self.variables:
                name = varname

        #print 'varname=',varname,'name=',name
        # Get the value of the abinit variable
        if name != '':
            npvalue = self.variables[name]
            if not isinstance(delta, int):
                npvalue += fact * delta

        elif (varname + ":") in self.variables and (varname + "+") in self.variables:
            npvalue = self.variables[varname + ":"] + (idtset - 1) * self.variables[varname + "+"]
        else:
            npvalue = None
        # Remove the array if it is one one element
        if npvalue is not None:
            if len(npvalue) == 1 and not full:
                value = npvalue[0]
            else:
                value = npvalue
        else:
            value = None
        # Return the value
        return value

    def set_value(self, varname, value, idtset=''):
        """
        Set the value 'value' into the dictionary
        input with key 'varname'+str(idtset)
        The value could be an integer, real, list
        or numpy array, the internal representation
        is always a numpy array.
        """
        if isinstance(value, (int, float)):
            npvalue = _np.array([value])
        else:
            npvalue = _np.array(value)
        if idtset == '':
            self.variables[varname] = npvalue
        else:
            self.variables[varname + str(idtset)] = npvalue

    def get_crystal(self, idtset='', units='bohr'):
        """
        Return the atomic structure from the input object
        for a given dataset (no dataset by default)
        """
        # NATOM
        natom = self.get_value('natom', idtset)
        # SYMBOLS
        ntypat = self.get_value('ntypat', idtset)
        symbols = []
        for i in range(natom):
            # NOTE: znucl is a real number in OUT.nc
            # Alchemic mixing is not allow here
            znucl = self.get_value('znucl', idtset)
            typat = self.get_value('typat', idtset)
            if ntypat == 1:
                symbols.append(atomic_symbol(int(znucl)))
            else:
                symbols.append(atomic_symbol(int(znucl[typat[i] - 1])))
        # POSITIONS
        xangst = _np.array(self.get_value('xangst', idtset))
        xcart = _np.array(self.get_value('xcart', idtset))
        xred = _np.array(self.get_value('xred', idtset))
        rprim = self.get_value('rprim', idtset)
        acell = _np.array(self.get_value('acell', idtset))

        # Set rprimd and acell using the default values
        # if not found
        if rprim is None:
            rprim = _np.identity(3)
        else:
            rprim = rprim.reshape(3, 3)
        if acell is None:
            acell = _np.ones(3)

        rprimd = _np.zeros((3, 3))
        rprimd[0] = rprim[0] * acell
        rprimd[1] = rprim[1] * acell
        rprimd[2] = rprim[2] * acell

        if xangst is not None:
            positions = xangst.reshape((natom, 3))
            if units == 'bohr':
                positions = positions * angstrom_bohr
        elif xcart is not None:
            positions = xcart.reshape((natom, 3))
            if units == 'angstrom':
                positions = positions * bohr_angstrom
        elif xred is not None:
            xred = xred.reshape(3, natom)
            xcart = _np.zeros(3, natom)

            for i in range(natom):
                xcart[i] = xred[0, i] * rprimd[0] + xred[1, i] * rprimd[1] + xred[2, i] * rprimd[2]

            positions = xcart
            if units == 'angstrom':
                positions = positions * bohr_angstrom

        else:
            positions = None

        if units == 'angstrom':
            rprimd = rprimd * bohr_angstrom

        # Create an object atomic_structure
        crystal = Structure(natom=natom, symbols=symbols, positions=positions, cell=rprimd)

        return crystal

    def from_structure(self, structure):
        """
        Set input variables for a given structure

        :return:
        """
        natom = structure.natom
        ntypat = len(structure.species)
        znucl = atomic_number(structure.species)
        typat_dict = {}
        index = 1
        for ispec in structure.species:
            typat_dict[ispec] = index
            index += 1
        typat = [typat_dict[i] for i in structure.symbols]
        xcart = pychemia.utils.constants.angstrom_bohr * structure.positions.flatten()
        acell = pychemia.utils.constants.angstrom_bohr * length_vectors(structure.cell)
        rprim = unit_vectors(structure.cell).T.flatten()
        for i in ['natom', 'ntypat', 'znucl', 'typat', 'xcart', 'acell', 'rprim']:
            self.set_value(i, eval(i))

    def get_dtsets_keys(self):
        """
        Get the list of dtset suffix acording to
        the values given by ndtset, jdtset and udtset
        """
        ret = None
        if 'jdtset' in self.variables and 'udtset' in self.variables:
            print('ERROR: Both udtset and jdtset could not be used ')
            return None

        elif 'ndtset' in self.variables:
            ndtset = self.get_value('ndtset')

            if ndtset != 0:
                if 'jdtset' in self.variables:
                    ret = list(self.variables['jdtset'])
                elif 'udtset' in self.variables:
                    ret = []
                    udtset = self.get_value('udtset')
                    for i in range(1, udtset[0] + 1):
                        for j in range(1, udtset[1] + 1):
                            ret.append(str(i) + str(j))
                else:
                    ret = range(1, ndtset + 1)
        else:
            ret = ['']

        return ret

    def atom_name(self, iatom, idtset=''):
        """
        Return the name of the atom and the position in the
        list of atoms such as H3, N4, etc
        """
        atomnumber = self.get_value('znucl', idtset=idtset, full=True)[
            self.get_value('typat', idtset=idtset, full=True)[iatom] - 1]
        return atomic_symbol(atomnumber) + str(iatom + 1)

    def view_projections(self):
        """
        Show the 3 projections of the molecule in a single
        figure
        """
        import matplotlib.patches as mpatches
        from matplotlib.collections import PatchCollection
        from matplotlib.pylab import subplots

        fig, ax = subplots(nrows=1, ncols=3)
        fig.set_size_inches(15, 4)
        color = ['r', 'g', 'b']
        j = 0

        cris = self.get_crystal()

        for i in cris.cell:
            ax[0].plot([0, i[0]], [0, i[1]], color[j] + '-', lw=3)
            ax[1].plot([0, i[1]], [0, i[2]], color[j] + '-', lw=3)
            ax[2].plot([0, i[2]], [0, i[0]], color[j] + '-', lw=3)
            j += 1

        proj = [[0, 1], [1, 2], [2, 0]]
        labels = [['x', 'y'], ['y', 'z'], ['z', 'x']]

        for j in range(3):

            patches = []
            for i in range(cris.natom):
                radius = 0.5 * covalent_radius(atomic_number(cris.symbols[i]))
                pos = cris.positions[i]
                art = mpatches.Circle((pos[proj[j][0]], pos[proj[j][1]]), radius, fc='g', ec='g')
                patches.append(art)

            collection = PatchCollection(patches, color='k', alpha=0.5)

            col = ax[j].add_collection(collection)
            ax[j].set_xlim(min(cris.positions[:, proj[j][0]]) - 1, max(cris.positions[:, proj[j][0]]) + 1)
            ax[j].set_ylim(min(cris.positions[:, proj[j][1]]) - 1, max(cris.positions[:, proj[j][1]]) + 1)
            ax[j].set_aspect('equal', adjustable='datalim')
            ax[j].set_xlabel(labels[j][0])
            ax[j].set_ylabel(labels[j][1])

        return fig, ax


def netcdf2dict(filename):
    """
    Read a NetCDF file and create a python dictionary with
    numpy arrays for variables

    Args:
        filename:
            NetCDF filename
    """
    if not _os.path.isfile(filename):
        print('ERROR: No such file: ', filename)
        return None
    output = {}
    netcdfile = _netcdf_file(filename, 'r', mmap=False)
    for ii in netcdfile.variables.keys():
        output[ii] = netcdfile.variables[ii][:]
    netcdfile.close()
    return output
