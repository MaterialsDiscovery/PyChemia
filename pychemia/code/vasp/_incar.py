from __future__ import unicode_literals
from builtins import str
import math
import os
import re
from numbers import Number
from pychemia.utils.computing import deep_unicode

import numpy as np

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "May 13, 2016"


def read_incar(filename='INCAR'):
    """
    Load the file INCAR in the directory 'path' or
    read directly the file 'path' and return an object
    'inputvars' for pychemia

    :param filename: (str) Filename of a INCAR file format
    :return:
    """
    if os.path.isfile(filename):
        filename = filename
    elif os.path.isdir(filename) and os.path.isfile(filename + '/INCAR'):
        filename += '/INCAR'
    else:
        print('INCAR path not found on ', filename)
        return

    iv = InputVariables(filename=filename)
    return iv


def write_incar(iv, filepath='INCAR'):
    """
    Takes an object inputvars from pychemia and
    save the file INCAR in the directory 'path' or
    save the file 'path' as a VASP INCAR file

    :param iv: (InputVariables) VASP Input variables
    :param filepath: (str) File path to write the INCAR file
    """

    if os.path.isdir(filepath):
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

        if filename is not None and os.path.isfile(filename):
            try:
                self.__import_input(filename)
            except ValueError:
                print('File format not identified')

    def _clean_variables(self):
        for i in self.variables:
            value = self.variables[i]
            if isinstance(value, list):
                if len(value) == 1:
                    self.variables[i] = value[0]

    def __import_input(self, filename):
        rf = open(filename, 'r')
        for line in rf.readlines():
            line = line.partition('#')[0]
            line = line.rstrip()
            if '=' in line:
                varname = line.split('=')[0].strip().upper()
                value = line.split('=')[1].strip()
                if value[-1] == ';':
                    value = value[:-1]
                if varname == '$SYSTEM':
                    self.variables[varname] = value
                elif value.upper() == '.FALSE.':
                    self.variables[varname] = False
                elif value.upper() == '.TRUE.':
                    self.variables[varname] = True
                else:
                    try:
                        self.variables[varname] = int(value)
                    except ValueError:
                        try:
                            self.variables[varname] = round(float(value), 10)
                        except ValueError:
                            self.variables[varname] = self._deep_parsing(value)
        rf.close()
        self._clean_variables()

    @staticmethod
    def _deep_parsing(value):

        # Try splitting
        val_splt = value.split()
        ret = []
        if len(val_splt) == 1:
            ret = value
        else:
            for i in range(len(val_splt)):
                try:
                    newval = [int(val_splt[i])]
                except ValueError:
                    try:
                        newval = [round(float(val_splt[i]), 10)]
                    except ValueError:
                        if '*' in val_splt[i] and len(val_splt[i].split('*')) == 2:
                            # print 'Trying this', val_splt[i]
                            newval = int(val_splt[i].split('*')[0]) * [float(val_splt[i].split('*')[1])]
                        else:
                            newval = [val_splt[i]]
                ret += newval
        return ret

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
        """
        ret = (varname.ljust(15)) + " =  "
        if varname not in self.variables:
            print("[ERROR] input variable: '%s' contains no elements" % varname)
            return

        value = self.variables[varname]
        value = deep_unicode(value)
        if isinstance(value, bool):
            if value:
                ret += '.TRUE.'
            else:
                ret += '.FALSE.'
        elif isinstance(value, Number):
            ret += str(value)
        elif isinstance(value, str):
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

            if len(self.variables[varname]) > 1:

                tmp = [(self.variables[varname][0], 1)]
                prev = self.variables[varname][0]
                for i in self.variables[varname][1:]:
                    if i == prev:
                        tmp[-1] = (i, tmp[-1][1] + 1)
                    else:
                        tmp.append((i, 1))
                    prev = i

                counter = 0
                for j in tmp:

                    if j[1] > 3:
                        if real:
                            if compact:
                                ret += (" %d*%g" % (j[1] - j[1] % 3, j[0])).rjust(8)
                            else:
                                ret += " %d*%g" % (j[1] - j[1] % 3, j[0])
                        elif integer:
                            ret += " %d*%d" % (j[1] - j[1] % 3, j[0])
                        else:
                            ret += " %d*%s" % (j[1] - j[1] % 3, j[0])

                        if j[1] % 3 != 0:
                            for i in range(j[1] % 3):
                                if real:
                                    if compact:
                                        ret += (" %g" % j[0]).rjust(8)
                                    else:
                                        ret += " %17.10e" % j[0]
                                elif integer:
                                    ret += " %d" % j[0]
                                else:
                                    ret += " %s" % j[0]

                        counter += j[1]
                    else:
                        for i in range(j[1]):
                            if real:
                                if compact:
                                    ret += (" %g" % j[0]).rjust(8)
                                else:
                                    ret += " %17.10e" % j[0]
                            elif integer:
                                ret += " %d" % j[0]
                            elif string:
                                ret += " %s" % j[0]

                            counter += 1

        ret += ";\n"
        return ret

    def set_encut(self, ENCUT=300, POTCAR=None):
        self.variables['ENCUT'] = ENCUT
        if POTCAR is not None and ENCUT < 10:
            maxvalue = 0
            if not os.path.isfile(POTCAR):
                raise ValueError('Not such file', POTCAR)
            rf = open(POTCAR)
            for line in rf.readlines():
                if 'ENMAX' in line:
                    list4line = line.split()
                    assert (list4line[0].strip() == 'ENMAX')
                    value = list4line[2].strip()
                    if value[-1] == ';':
                        value = value[:-1]
                    value = float(value)
                    if value > maxvalue:
                        maxvalue = value
            rf.close()
            if ENCUT < 10:
                self.variables['ENCUT'] = int(math.ceil(ENCUT * maxvalue))
            else:
                self.variables['ENCUT'] = maxvalue
                # pcm_log.debug('ENCUT: %7.3f' % self.variables['ENCUT'])

    def set_ismear(self, kpoints):
        if np.prod(kpoints.grid) > 27:
            self.variables['ISMEAR'] = -5
        else:
            self.variables['ISMEAR'] = 0

    def set_ion_relax(self, NSW=50, ISIF=2, IBRION=2, EDIFFG=-1E-3):
        self.variables['IBRION'] = IBRION
        self.variables['NSW'] = NSW
        self.variables['ISIF'] = ISIF
        self.variables['EDIFFG'] = EDIFFG

    def set_electron_scf(self, NELM=60, NELMIN=2, EDIFF=1E-4, IALGO=48):
        self.variables['NELMIN'] = NELMIN
        self.variables['NELM'] = NELM
        self.variables['EDIFF'] = EDIFF
        self.variables['IALGO'] = IALGO

    def set_rough_relaxation(self):
        self.set_minimum(PREC='Normal', ISPIN=1, LREAL=False, ISMEAR=0, LORBIT=11)
        self.set_electron_scf(NELM=60, NELMIN=2, EDIFF=1E-4, IALGO=48)
        self.set_ion_relax(NSW=50, ISIF=2, IBRION=2, EDIFFG=-1E-2)
        self.variables['NPAR'] = 2

    def set_mit_settings(self):
        self.set_minimum(PREC='Accurate', ISPIN=2, LREAL=False, ISMEAR=-5, LORBIT=11)
        self.set_electron_scf(NELM=100, NELMIN=6, EDIFF=5E-5, IALGO=48)
        self.set_ion_relax(NSW=99, ISIF=3, IBRION=2, EDIFFG=-1E-3)
        self.variables['LWAVE'] = False
        self.variables['SIGMA'] = 0.05
        self.variables['LDAU'] = True
        self.variables['LDAUTYPE'] = 2
        self.variables['ICHARG'] = 1

    def set_minimum(self, PREC='Normal', ISPIN=2, LREAL=False, ISMEAR=0, LORBIT=11):
        self.variables['PREC'] = PREC
        if LREAL is not None:
            self.variables['LREAL'] = LREAL
        else:
            self.variables['LREAL'] = 'Auto'
        self.variables['ISMEAR'] = ISMEAR
        self.variables['ISPIN'] = ISPIN
        self.variables['LORBIT'] = LORBIT
        self.variables['NPAR'] = 2
        self.variables['LASPH'] = True

    def set_density_for_restart(self):
        self.variables['ICHARG'] = 1

    def get_value(self, varname, return_iterable=False):

        if varname not in self.variables:
            raise ValueError('Input variables does not contain value for "%s"' % varname)

        value = self.variables[varname]
        value = deep_unicode(value)
        if return_iterable and isinstance(value, (int, float, str)):
            value = [value]

        return value


def get_potcar_info(filename='POTCAR'):
    rf = open(filename)
    data = rf.read()
    ret = re.findall('([\w ]*)=([.\d ]*)', data)
    return ret
