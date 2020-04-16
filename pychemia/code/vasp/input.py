import os
import math
from numbers import Number
import numpy as np
from pychemia.utils.computing import deep_unicode
from ..codes import CodeInput


class VaspInput(CodeInput):
    """
    VASP INCAR object

    It contains:

    data:
          variables = Dictionary whose keys are VASP variable names
                      and contains the values as numpy arrays

    methods:
            write = Write the input into as a text file that VASP
                    can use as an input file

            get_value = Get the value of a particular variable
            set_value = Set the value of a particular variable
    """
    def __init__(self, filename=None, variables=None):

        # CodeInput.__init__(self)
        if variables is not None:
            for i in variables:
                self.__dict__[i] = variables[i]

        if filename is not None and os.path.isfile(filename):
            try:
                self.__import_input(filename)
            except ValueError:
                print('File format not identified')

    def _clean_variables(self):
        for i in self.__dict__:
            value = self.__dict__[i]
            if isinstance(value, list):
                if len(value) == 1:
                    self.__dict__[i] = value[0]

    @property
    def variables(self):
        return self.__dict__

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
                    self.__dict__[varname] = value
                elif value.upper() == '.FALSE.':
                    self.__dict__[varname] = False
                elif value.upper() == '.TRUE.':
                    self.__dict__[varname] = True
                else:
                    try:
                        self.__dict__[varname] = int(value)
                    except ValueError:
                        try:
                            self.__dict__[varname] = round(float(value), 10)
                        except ValueError:
                            self.__dict__[varname] = self._deep_parsing(value)
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

    def read(self, filename='INCAR'):
        pass

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

    def write_key(self, varname):
        """
        Receives an input variable and write their contents
        properly according with their kind and length

        Args:
            varname:
                The name of the input variable
        """
        ret = (varname.ljust(15)) + " =  "
        if varname not in self.__dict__:
            raise ValueError("[ERROR] input variable: '%s' is not declared" % varname)

        value = self.__dict__[varname]
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
            for j in self.__dict__[varname]:

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

            if len(self.__dict__[varname]) == 1:
                ret += "%s" % self.__dict__[varname][0]
            if len(self.__dict__[varname]) > 1:
                tmp = [(self.__dict__[varname][0], 1)]
                prev = self.__dict__[varname][0]
                for i in self.__dict__[varname][1:]:
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
        self.__dict__['ENCUT'] = ENCUT
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
                self.__dict__['ENCUT'] = int(math.ceil(ENCUT * maxvalue))
            else:
                self.__dict__['ENCUT'] = maxvalue
                # pcm_log.debug('ENCUT: %7.3f' % self.variables['ENCUT'])

    def set_ismear(self, kpoints):
        if np.prod(kpoints.grid) > 27:
            self.__dict__['ISMEAR'] = -5
        else:
            self.__dict__['ISMEAR'] = 0

    def set_ion_relax(self, NSW=50, ISIF=2, IBRION=2, EDIFFG=-1E-3):
        self.__dict__['IBRION'] = IBRION
        self.__dict__['NSW'] = NSW
        self.__dict__['ISIF'] = ISIF
        self.__dict__['EDIFFG'] = EDIFFG

    def set_electron_scf(self, NELM=60, NELMIN=2, EDIFF=1E-4, IALGO=48):
        self.__dict__['NELMIN'] = NELMIN
        self.__dict__['NELM'] = NELM
        self.__dict__['EDIFF'] = EDIFF
        self.__dict__['IALGO'] = IALGO

    def set_rough_relaxation(self):
        self.set_minimum(PREC='Normal', ISPIN=1, LREAL=False, ISMEAR=0, LORBIT=11)
        self.set_electron_scf(NELM=60, NELMIN=2, EDIFF=1E-4, IALGO=48)
        self.set_ion_relax(NSW=50, ISIF=2, IBRION=2, EDIFFG=-1E-2)
        # self.__dict__['NPAR'] = 2

    def set_mit_settings(self):
        self.set_minimum(PREC='Accurate', ISPIN=2, LREAL=False, ISMEAR=-5, LORBIT=11)
        self.set_electron_scf(NELM=100, NELMIN=6, EDIFF=5E-5, IALGO=48)
        self.set_ion_relax(NSW=99, ISIF=3, IBRION=2, EDIFFG=-1E-3)
        self.__dict__['LWAVE'] = False
        self.__dict__['SIGMA'] = 0.05
        self.__dict__['LDAU'] = True
        self.__dict__['LDAUTYPE'] = 2
        self.__dict__['ICHARG'] = 1

    def set_minimum(self, PREC='Normal', ISPIN=2, LREAL=False, ISMEAR=0, LORBIT=11):
        self.__dict__['PREC'] = PREC
        if LREAL is not None:
            self.__dict__['LREAL'] = LREAL
        else:
            self.__dict__['LREAL'] = 'Auto'
        self.__dict__['ISMEAR'] = ISMEAR
        self.__dict__['ISPIN'] = ISPIN
        self.__dict__['LORBIT'] = LORBIT
        # self.__dict__['NPAR'] = 2
        self.__dict__['LASPH'] = True

    def set_density_for_restart(self):
        self.__dict__['ICHARG'] = 1

    def get_defined_variables(self):
        return list(self.keys())

    # The next five methods are requirements of the ABC.
    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __delitem__(self, key):
        del self.__dict__[key]

    def __iter__(self):
        return iter(self.__dict__)

    def __len__(self):
        return len(self.__dict__)

    # The final two methods aren't required, but nice for demo purposes:
    def __str__(self):
        """
        returns simple dict representation of the mapping
        """
        ret = ''
        for i in sorted(self.__dict__.keys()):
            ret += self.write_key(i)
        return ret

    def __repr__(self):
        """
        echoes class, id, & reproducible representation in the REPL
        """
        return '{}, {}(variables={})'.format(super(self.__class__, self).__repr__(),
                                             self.__class__.__name__,
                                             self.__dict__)



