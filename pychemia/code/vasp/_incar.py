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
import numpy as _np
from numbers import Number


def load_INCAR(path):
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


def save_INCAR(iv, path):
    """
    Takes an object inputvars from pychemia and
    save the file INCAR in the directory 'path' or
    save the file 'path' as a VASP INCAR file
    """

    if _os.path.isdir(path):
        filename = path + '/INCAR'
    else:
        filename = path

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

        if filename is not None and _os.path.isfile(filename):
            try:
                self.__import_input(filename)
            except ValueError:
                print('File format not identified')

        if variables is not None:
            self.variables = variables

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
