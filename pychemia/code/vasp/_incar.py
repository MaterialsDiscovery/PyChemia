"""
Routines to read and write INCAR file
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "March 16, 2014"

import os as _os
import math
import numpy as _np
from numbers import Number


def read_incar(path):
    """
    Load the file INCAR in the directory 'path' or
    read directly the file 'path' and return an object
    'inputvars' for pychemia
    """
    if _os.path.isfile(path):
        filename = path
    elif _os.path.isdir(path) and _os.path.isfile(path + '/INCAR'):
        filename = path + '/INCAR'
    else:
        print('INCAR path not found on ', path)
        return

    iv = InputVariables(filename=filename)
    return iv


def write_incar(iv, filepath='INCAR'):
    """
    Takes an object inputvars from pychemia and
    save the file INCAR in the directory 'path' or
    save the file 'path' as a VASP INCAR file
    """

    if _os.path.isdir(filepath):
        filename = filepath + '/INCAR'
    else:
        filename = filepath

    iv.write(filename)


class InputVariables:
    """
    VASP INCAR object

    It contains:

    data:
          variables = Dictionary whose keys are ABINIT variable names
                      and contains the values as numpy arrays

    methods:
            write = Write the input into as a text file that ABINIT
                    can use as an input file

            get_value = Get the value of a particular variable
            set_value = Set the value of a particular variable
    """

    variables = {}

    def __init__(self, filename=None, variables=None):

        if variables is None:
            self.variables = {}
        else:
            self.variables = variables

        if filename is not None and _os.path.isfile(filename):
            try:
                self.__import_input(filename)
            except ValueError:
                print('File format not identified')

    def _clean_variables(self):
        for i in self.variables:
            value = self.variables[i]
            if isinstance(value, _np.ndarray):
                if len(value) == 1:
                    self.variables[i]=value[0]
                else:
                    self.variables[i] = list(value)

    def __import_input(self, filename):
        rf = open(filename, 'r')
        for line in rf.readlines():
            if '=' in line:
                varname = line.split('=')[0].strip().upper()
                value = line.split('=')[1].strip()
                if value[-1] == ';':
                    value = value[:-1]
                if value.upper() == '.FALSE.':
                    self.variables[varname] = False
                elif value.upper() == '.TRUE.':
                    self.variables[varname] = True
                else:
                    try:
                        self.variables[varname] = _np.array([int(value)])
                    except ValueError:
                        try:
                            self.variables[varname] = _np.array([float(value)])
                        except ValueError:
                            self.variables[varname] = _np.array([value])
        rf.close()
        self._clean_variables()

    def write(self, filename='INCAR'):
        """
        Write an inputvars object into a text
        file that VASP can use as an INCAR
        file

        Args:
            filename:
                The path to 'INCAR' filename that will be written
        """
        wf = open(filename, 'w')
        wf.write(self.__str__())
        wf.close()

    def __str__(self):

        ret = ''
        thekeys = self.variables.keys()

        for i in sorted(thekeys):
            ret += self.write_key(i)

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
        ret = (varname.ljust(15)) + " =  "
        if varname not in self.variables:
            print("[ERROR] input variable: '%s' contains no elements" % varname)
            return

        value = self.variables[varname]
        if isinstance(value, bool):
            if value:
                ret += '.TRUE.'
            else:
                ret += '.FALSE.'
        elif isinstance(value, Number):
            ret += str(value)
        elif isinstance(value, basestring):
            ret += value
        else:

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
                    # be converted because we do not know the size
                    # of the array
                    integer = False
                    real = False
                    string = True

            for j in range(len(self.variables[varname])):

                if real:
                    if compact:
                        ret += ("%g" % self.variables[varname][j]).rjust(8)
                    else:
                        ret += ("%17.10e" % self.variables[varname][j])
                elif integer:
                    ret += ("%d" % self.variables[varname][j])
                elif string:
                    ret += ("%s" % self.variables[varname][j])

                # Conditions to jump to a new line
                if ((j + 1) % 3) == 0 and real and j < len(self.variables[varname]) - 1:
                    ret += ";\n"
                    ret += 17 * " "
                elif j < len(self.variables[varname]) - 1:
                    ret += " "

        ret += ";\n"
        return ret

    def set_minimum(self, PREC='Normal', ISPIN=2,  LREAL=False, ISMEAR=0, LORBIT=11):
        self.variables['PREC'] = PREC
        self.variables['LREAL'] = LREAL
        self.variables['ISMEAR'] = ISMEAR
        self.variables['ISPIN'] = ISPIN
        self.variables['LORBIT'] = LORBIT

    def set_encut(self, ENCUT=300, POTCAR=None):
        self.variables['ENCUT'] = ENCUT
        if POTCAR is not None and ENCUT < 10:
            maxvalue = 0
            if not _os.path.isfile(POTCAR):
                raise ValueError('Not such file', POTCAR)
            rf = open(POTCAR)
            for line in rf.readlines():
                if 'ENMAX' in line:
                    list4line = line.split()
                    assert(list4line[0].strip() == 'ENMAX')
                    value = list4line[2].strip()
                    if value[-1] == ';':
                        value = value[:-1]
                    value = float(value)
                    if value > maxvalue:
                        maxvalue = value
            rf.close()
            if ENCUT < 10:
                self.variables['ENCUT'] = int(math.ceil(ENCUT*maxvalue))
            else:
                self.variables['ENCUT'] = maxvalue
        print self.variables['ENCUT']

    def set_ion_relax(self, NSW=50, ISIF=2):
        self.variables['IBRION'] = 2
        self.variables['NSW'] = NSW
        self.variables['ISIF'] = ISIF

    def set_break_conditions(self, EDIFF='1E-4', EDIFFG='-1E-3'):
        self.variables['EDIFF'] = EDIFF
        self.variables['EDIFFG'] = EDIFFG

    def set_rough_relaxation(self):
        self.variables['NELMIN'] = 5
        self.variables['EDIFF'] = 1E-4
        self.variables['EDIFFG'] = -1E-3
        self.variables['IBRION'] = 2
        self.variables['IALGO'] = 48
