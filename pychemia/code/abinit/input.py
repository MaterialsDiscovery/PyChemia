
import os
import numpy as np
from pychemia.utils.periodic import atomic_symbol, covalent_radius, atomic_number
from pychemia.utils.constants import bohr_angstrom, angstrom_bohr
from pychemia.utils.mathematics import unit_vectors
from .parser import parser
from pychemia.utils.netcdf import netcdf2dict
from pychemia.core import Structure
from ..codes import CodeInput

"""
Definition of the class input to read
ABINIT input files and store their information
as a python dictionary called 'variables'
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2016"
__version__ = "1.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "guillermo.avendano@uclouvain.be"
__status__ = "Development"
__date__ = "May 13, 2016"


class AbinitInput(CodeInput):

    def __init__(self, input_filename=None):
        CodeInput.__init__(self)
        if input_filename is None:
            self.input_file = 'abinit.in'
        else:
            self.input_file = input_filename
        if os.path.isfile(self.input_file):
            if self.input_file[-6:] == 'OUT.nc':
                self.variables = netcdf2dict(self.input_file)
            else:
                self.read()

    def __import_input(self, filename):
        """
        Read an ABINIT input file and return a python dictionary
        with the input variables readed from that. The keys are
        the fullname input variables (acell,xcart21,etc). The
        values are numbers or lists except for
        the value '*[NUMBER]' that is keeped as string, and
        the string associated to the variable xyzfile

        Args:
            filename:
                ABINIT input filename
        """
        ans = parser(filename)
        if ans is not None:
            self.variables = ans

    def read(self):
        if not os.path.isfile(self.input_file):
            raise ValueError("ERROR: Could not read %s" % self.input_file)

        ans = parser(self.input_file)
        if ans is not None:
            self.variables = ans

    def check(self):
        if self.get_number_variables > 0:
            print("ABINIT input is readable and has %d variables" % self.get_number_variables)

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
                       (x[-1] == ':' or x[-1] == "+" or x[-1] == "?" or x[-2] == ':' or x[-2] == "+" or x[
                           -2] == "?")]
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
                if i == 'dmatpawu' and 'lpawu' in self.variables:
                    if 2 in self.variables['lpawu']:
                        ret += self.write_key(i, ncolumns=5)
                    elif 3 in self.variables['lpawu']:
                        ret += self.write_key(i, ncolumns=7)
                else:
                    ret += self.write_key(i)
            ret += '\n'

        for dtset in range(1, 100):
            varnames = [x for x in thekeys if
                        (x[-len(str(dtset)):] == str(dtset) and not x[-len(str(dtset)) - 1:].isdigit())]
            if len(varnames) > 0:
                varnames.sort()
                ret = ret + "#" + 60 * "-" + "\n#" + " DATASET " + str(dtset) + "\n#" + 60 * "-" + "\n\n"
                for i in varnames:
                    if i == 'dmatpawu' and 'lpawu' in self.variables:
                        if 2 in self.variables['lpawu']:
                            ret += self.write_key(i, ncolumns=5)
                        elif 3 in self.variables['lpawu']:
                            ret += self.write_key(i, ncolumns=7)
                    ret += self.write_key(i)
                ret += '\n'
        return ret

    def write_key(self, varname, ncolumns=None, debug=False):
        """
        Receives an input variable and write their contents
        properly according with their kind and length

        Args:
            varname:
                The name of the input variable
            ncolumns:
                Number of columns for the input variable
            debug:
                Shows contents of variable before creating string out of it
        """

        ret = ''
        if varname not in self.variables:
            print("[ERROR] input variable: '%s' is not defined" % varname)
            return

        # Assume that the variables are integer and test if such assumption
        # is true
        integer = True
        real = False
        string = False
        compact = True

        if isinstance(self.variables[varname], (int, float)):
            varlist = [self.variables[varname]]
        elif isinstance(self.variables[varname], str):
            varlist = [self.variables[varname]]
        else:
            varlist = self.variables[varname]

        if debug:
            print('varlist: %s' % varlist)
            
        # Get the general kind of values for the input variable
        for j in varlist:
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

        known_variables = {'xred': [3], 'acell': [3]}

        if varname in known_variables:
            for i in known_variables[varname]:
                if len(varlist) % i == 0:
                    for j in range(int(len(varlist) / i)):
                        if j == 0:
                            ret += (i * '%17.10E ' + '\n') % tuple(varlist[j * i:j * i + i])
                        else:
                            ret += (17 * ' ' + i * '%17.10E ' + '\n') % tuple(varlist[j * i:j * i + i])
        elif ncolumns is not None:
            i = ncolumns
            for j in range(int(len(varlist) / i)):
                if j == 0:
                    ret += (i * '%17.10E ' + '\n') % tuple(varlist[j * i:j * i + i])
                else:
                    ret += (17 * ' ' + i * '%17.10E ' + '\n') % tuple(varlist[j * i:j * i + i])
        else:
            if debug:
                print("real: %s  integer: %s  string: %s" % (real, integer, string))

            for j in range(len(varlist)):

                if real:
                    if compact:
                        ret += ("%g" % varlist[j]).rjust(8)
                    else:
                        ret += "%17.10e" % varlist[j]
                elif integer:
                    ret += "%d" % varlist[j]
                elif string:
                    ret += "%s" % varlist[j]

                # Conditions to jump to a new line
                if ((j + 1) % 3) == 0 and real and j < len(varlist) - 1:
                    ret += "\n"
                    ret += 17 * " "
                elif j < len(varlist) - 1:
                    ret += " "
            ret += "\n"
        return ret

    def get_structure(self, idtset=None, units='bohr'):
        """
        Return the atomic structure from the input object
        for a given dataset (no dataset by default)
        """
        # NATOM
        natom = self.get_value('natom', idtset)
        # SYMBOLS
        ntypat = self.get_value('ntypat', idtset)
        symbols = []
        znucl = self.get_value('znucl', idtset)
        typat = self.get_value('typat', idtset)
        for i in range(natom):
            # NOTE: znucl is a real number in OUT.nc
            # Alchemic mixing is not allow here
            if ntypat == 1:
                symbols.append(atomic_symbol(int(znucl)))
            else:
                symbols.append(atomic_symbol(int(znucl[typat[i] - 1])))

        # POSITIONS
        xangst = self.get_value('xangst', idtset)
        xcart = self.get_value('xcart', idtset)
        xred = self.get_value('xred', idtset)

        rprim = self.get_value('rprim', idtset)
        acell = np.array(self.get_value('acell', idtset))

        # Set rprimd and acell using the default values
        # if not found
        if rprim is None:
            rprim = np.identity(3)
        else:
            rprim = np.array(rprim).reshape((3, 3))
        if acell is None:
            acell = np.ones(3)

        rprimd = np.zeros((3, 3))
        rprimd[0] = rprim[0] * acell
        rprimd[1] = rprim[1] * acell
        rprimd[2] = rprim[2] * acell

        if xangst is not None:
            xangst = np.array(xangst)
            positions = xangst.reshape((natom, 3))
            if units == 'bohr':
                positions = positions * angstrom_bohr
        elif xcart is not None:
            xcart = np.array(xcart)
            positions = xcart.reshape((natom, 3))
            if units == 'angstrom':
                positions = positions * bohr_angstrom
        elif xred is not None:
            xred = np.array(xred)
            xred = xred.reshape((natom, 3))
            xcart = np.zeros((natom, 3))

            # print rprimd
            # print xred

            for i in range(natom):
                xcart[i] = xred[i, 0] * rprimd[0] + xred[i, 1] * rprimd[1] + xred[i, 2] * rprimd[2]

            positions = xcart
            if units == 'angstrom':
                positions = positions * bohr_angstrom

        else:
            positions = None

        if units == 'angstrom':
            rprimd = rprimd * bohr_angstrom

        # Create an object atomic_structure
        structure = Structure(natom=natom, symbols=symbols, positions=positions, cell=rprimd)

        return structure

    def from_structure(self, structure):
        """
        Set input variables for a given structure

        :param structure: (pychemia.Structure) Structure to set ABINIT input variables
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
        xcart = angstrom_bohr * structure.positions.flatten()
        acell = angstrom_bohr * np.array(structure.lattice.lengths)
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
                            print(ret)
                            print(str(i) + str(j))
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
        atomnumber = self.get_value('znucl', idtset=idtset)
        print("atomnumber=%s" % atomnumber)
        if isinstance(atomnumber, list):
            atomnumber = atomnumber[self.get_value('typat', idtset=idtset)[iatom] - 1]
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

        structure = self.get_structure()

        for i in structure.cell:
            ax[0].plot([0, i[0]], [0, i[1]], color[j] + '-', lw=3)
            ax[1].plot([0, i[1]], [0, i[2]], color[j] + '-', lw=3)
            ax[2].plot([0, i[2]], [0, i[0]], color[j] + '-', lw=3)
            j += 1

        proj = [[0, 1], [1, 2], [2, 0]]
        labels = [['x', 'y'], ['y', 'z'], ['z', 'x']]

        for j in range(3):

            patches = []
            for i in range(structure.natom):
                radius = 0.5 * covalent_radius(atomic_number(structure.symbols[i]))
                pos = structure.positions[i]
                art = mpatches.Circle((pos[proj[j][0]], pos[proj[j][1]]), radius, fc='g', ec='g')
                patches.append(art)

            collection = PatchCollection(patches, color='k', alpha=0.5)

            col = ax[j].add_collection(collection)
            ax[j].set_xlim(min(structure.positions[:, proj[j][0]]) - 1, max(structure.positions[:, proj[j][0]]) + 1)
            ax[j].set_ylim(min(structure.positions[:, proj[j][1]]) - 1, max(structure.positions[:, proj[j][1]]) + 1)
            ax[j].set_aspect('equal', adjustable='datalim')
            ax[j].set_xlabel(labels[j][0])
            ax[j].set_ylabel(labels[j][1])

        return fig, ax

    def clean(self):
        self.variables = {}

    def has_variable(self, varname, section=None):
        return varname in self.variables

    def get_value(self, varname, idtset=None, return_iterable=False):
        """
        Get the value of the input variable 'varname'
        associated with the dataset 'idtset'
        If 'idtset' is not given will asume that the
        value is not dataset dependent
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

        # print 'varname=',varname,'name=',name
        #  Get the value of the abinit variable
        if name != '':
            if isinstance(self.variables[name], list):
                npvalue = list(np.array(self.variables[name]) + fact * np.array(delta))
            else:
                npvalue = self.variables[name] + fact * delta

        elif (varname + ":") in self.variables and (varname + "+") in self.variables:
            if isinstance(self.variables[varname + ":"], list):
                npvalue = list(np.array(self.variables[varname + ":"]) +
                               (idtset - 1) * np.array(self.variables[varname + "+"]))
            else:
                npvalue = self.variables[varname + ":"] + (idtset - 1) * self.variables[varname + "+"]
        else:
            npvalue = None

        if isinstance(npvalue, (int, float)) and return_iterable:
            npvalue = [npvalue]

        return npvalue

    def set_value(self, varname, value, idtset=''):
        """
        Set the value 'value' into the dictionary
        input with key 'varname'+str(idtset)
        The value could be an integer, real, list
        or numpy array, the internal representation
        is always serializable.
        """
        if isinstance(value, (int, float)):
            npvalue = value
        elif isinstance(value, np.ndarray):
            if value[0].dtype == np.dtype('>f8'):
                npvalue = [round(x, 11) for x in value.flatten()]
            elif value[0].dtype == np.dtype('>i4'):
                npvalue = [int(x) for x in value.flatten()]
            else:
                npvalue = list(value)
        else:
            npvalue = list(value)

        if idtset == '':
            self.variables[varname] = npvalue
        else:
            self.variables[varname + str(idtset)] = npvalue


# class InputVariables(collections.MutableMapping):
#     """
#     An input object contains:

#     data:
#           variables = Dictionary whose keys are ABINIT variable names

#     methods:
#             write = Write the input into as a text file that ABINIT
#                     can use as an input file

#             get_value = Get the value of a particular variable
#             set_value = Set the value of a particular variable

#             get_atomic_structure =
#     """


#     def __init__(self, *args, **kwargs):
#         """
#         Creates a new input object, the input object
#         contains a dictionary whose keys are the abinit
#         input variables and the values are always serializable
#         """
#         self.variables = {}
#         filename = ''
#         if len(args) == 1:
#             x = args[0]
#             if isinstance(x, AbiFiles):
#                 filename = x.basedir + "/" + x.files["in"]
#             elif os.path.isfile(x):
#                 filename = x

#         if 'filename' in kwargs:
#             filename = kwargs['filename']

#         if os.path.isfile(filename):
#             if filename[-3:] == '.in':
#                 self.__import_input(filename)
#             elif filename[-6:] == '.files':
#                 abifile = AbiFiles(filename)
#                 filename = abifile.get_input_filename()
#                 self.__import_input(filename)
#             elif filename[-6:] == 'OUT.nc':
#                 self.variables = netcdf2dict(filename)
#             else:
#                 try:
#                     self.__import_input(filename)
#                 except ValueError:
#                     print('File format not identified')

#     def __delitem__(self, key):
#         return self.variables.__delitem__(key)

#     def __setitem__(self, key, value):
#         return self.variables.__setitem__(key, value)

#     def __getitem__(self, key):
#         return self.variables.__getitem__(key)

#     def __iter__(self):
#         return self.variables.__iter__()

#     def __len__(self):
#         return self.variables.__len__()

#     def __import_input(self, filename):
#         """
#         Read an ABINIT input file and return a python dictionary
#         with the input variables readed from that. The keys are
#         the fullname input variables (acell,xcart21,etc). The
#         values are numbers or lists except for
#         the value '*[NUMBER]' that is keeped as string, and
#         the string associated to the variable xyzfile

#         Args:
#             filename:
#                 ABINIT input filename
#         """
#         ans = parser(filename)
#         if ans is not None:
#             self.variables = ans

#     def write(self, filename):
#         """
#         Write an input object into a text
#         file that ABINIT can use as an input
#         file

#         Args:
#             filename:
#                 The 'abinit.in' filename that will be written
#         """
#         wf = open(filename, 'w')
#         wf.write(self.__str__())
#         wf.close()

#     def __str__(self):
#         """
#         String representation of the object
#         """
#         ret = ''
#         thekeys = self.variables.keys()
#         varnames = [x for x in thekeys if not x[-1].isdigit()]

#         if 'ndtset' in varnames:
#             ret = ret + "#" + 60 * "-" + "\n#" + " MULTI DATASET\n#" + 60 * "-" + "\n\n"
#             ret += self.write_key('ndtset')
#             varnames.remove('ndtset')
#             if 'jdtset' in varnames:
#                 ret += self.write_key('jdtset')
#                 varnames.remove('jdtset')
#             if 'udtset' in varnames:
#                 ret += self.write_key('udtset')
#                 varnames.remove('udtset')
#             ret += '\n'

#         seqvarnames = [x for x in varnames if
#                        (x[-1] == ':' or x[-1] == "+" or x[-1] == "?" or x[-2] == ':' or x[-2] == "+" or x[-2] == "?")]
#         if len(seqvarnames) > 0:
#             ret = ret + "#" + 60 * "-" + "\n#" + " SEQUENCE\n#" + 60 * "-" + "\n\n"
#             seqvarnames.sort()
#             for i in seqvarnames:
#                 ret += self.write_key(i)
#                 varnames.remove(i)
#             ret += '\n'

#         if len(varnames) > 0:
#             varnames.sort()
#             ret = ret + "#" + 60 * "-" + "\n#" + " ALL DATASETS\n#" + 60 * "-" + "\n\n"
#             for i in varnames:
#                 ret += self.write_key(i)
#             ret += '\n'

#         for dtset in range(1, 100):
#             varnames = [x for x in thekeys if
#                         (x[-len(str(dtset)):] == str(dtset) and not x[-len(str(dtset)) - 1:].isdigit())]
#             if len(varnames) > 0:
#                 varnames.sort()
#                 ret = ret + "#" + 60 * "-" + "\n#" + " DATASET " + str(dtset) + "\n#" + 60 * "-" + "\n\n"
#                 for i in varnames:
#                     if i == 'dmatpawu':
#                         if 2 in self['lpawu']:
#                             ret += self.write_key(i, ncolumns=5)
#                         elif 3 in self['lpawu']:
#                             ret += self.write_key(i, ncolumns=7)
#                     ret += self.write_key(i)
#                 ret += '\n'
#         return ret

#     def write_key(self, varname, ncolumns=None):
#         """
#         Receives an input variable and write their contents
#         properly according with their kind and length

#         Args:
#             varname:
#                 The name of the input variable
#         """

#         ret = ''
#         if varname not in self.variables:
#             print("[ERROR] input variable: '%s' is not defined" % varname)
#             return

#         # Assume that the variables are integer and test if such assumption
#         # is true
#         integer = True
#         real = False
#         string = False
#         compact = True

#         if isinstance(self.variables[varname], (int, float)):
#             varlist = [self.variables[varname]]
#         else:
#             varlist = self.variables[varname]

#         # Get the general kind of values for the input variable
#         for j in varlist:
#             try:
#                 if not float(j).is_integer():
#                     # This is the case of non integer values
#                     integer = False
#                     real = True
#                     string = False
#                     if len(str(float(j))) > 7:
#                         compact = False

#             except ValueError:
#                 # This is the case of '*1' that could not
#                 # be converted because we don't know the size
#                 # of the array
#                 integer = False
#                 real = False
#                 string = True

#         ret = ret + (varname.rjust(15)) + "  "

#         known_variables = {'xred': [3], 'acell': [3]}

#         if varname in known_variables:
#             for i in known_variables[varname]:
#                 if len(varlist) % i == 0:
#                     for j in range(int(len(varlist) / i)):
#                         if j == 0:
#                             ret += (i * '%17.10E ' + '\n') % tuple(varlist[j * i:j * i + i])
#                         else:
#                             ret += (17 * ' ' + i * '%17.10E ' + '\n') % tuple(varlist[j * i:j * i + i])
#         elif ncolumns is not None:
#             for i in ncolumns:
#                 for j in range(int(len(varlist) / i)):
#                     if j == 0:
#                         ret += (i * '%17.10E ' + '\n') % tuple(varlist[j * i:j * i + i])
#                     else:
#                         ret += (17 * ' ' + i * '%17.10E ' + '\n') % tuple(varlist[j * i:j * i + i])
#         else:
#             for j in range(len(varlist)):

#                 if real:
#                     if compact:
#                         ret += ("%g" % varlist[j]).rjust(8)
#                     else:
#                         ret += "%17.10e" % varlist[j]
#                 elif integer:
#                     ret += "%d" % varlist[j]
#                 elif string:
#                     ret += "%s" % varlist[j]

#                 # Conditions to jump to a new line
#                 if ((j + 1) % 3) == 0 and real and j < len(varlist) - 1:
#                     ret += "\n"
#                     ret += 17 * " "
#                 elif j < len(varlist) - 1:
#                     ret += " "
#             ret += "\n"
#         return ret


def merge(abi_into, abi_from, filename=None):
    """
    From the variables present in the file 'abi_into'
    will try to recover the value from 'abi_from'
    if the value has changed, the value will be updated.

    The new variables will be save in the given filename
    If filename==None the new values will overwrite
    abi_into

    Example:

       merge('abinit.in','abinit_xo_OUT.nc')

    It will update the abinit input with the values of the output

    :param abi_into: (str) Merge abinit variables using this destination
    :param abi_from:  (str) and this source
    :param filename: (str) Storing the final input variables on file

    """
    abinit_into = AbinitInput(abi_into)
    abinit_from = AbinitInput(abi_from)

    for i in abinit_into.variables.keys():
        if i in abinit_from.variables.keys():
            abinit_into.variables[i] = abinit_from.variables[i]

    if filename is None:
        filename = abi_into

    abinit_into.write(filename)


def xyz2input(filename):
    """
    Reads a .xyz and return an ABINIT input
    as a python dictionary
    """

    abiinput = AbinitInput()
    atomdict = atomic_symbol()
    rf = open(filename, 'r')

    natom = int(rf.readline())
    typat = []
    znucl = []
    xangst = []

    ntypat = 0
    rf.readline()
    data = rf.readlines()
    for i in range(natom):
        atom = data[i].split()
        atomnumber = atomdict[atom[0]]
        if atomnumber not in znucl:
            ntypat += 1
            znucl.append(atomnumber)
        typat.append(znucl.index(atomnumber) + 1)
        xangst += [float(atom[1]), float(atom[2]), float(atom[3])]

    abiinput.variables['natom'] = np.array([natom])
    abiinput.variables['znucl'] = np.array(znucl)
    abiinput.variables['ntypat'] = np.array([ntypat])
    abiinput.variables['typat'] = np.array(typat)
    abiinput.variables['xcart'] = angstrom_bohr * np.array(xangst)

    return abiinput
