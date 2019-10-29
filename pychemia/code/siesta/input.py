from ..codes import CodeInput
from pychemia.utils.periodic import atomic_symbol
from pychemia import Structure
import os
import numpy as np


class SiestaInput(CodeInput):

    def __init__(self, inputfile=None):
        CodeInput.__init__(self)
        if inputfile is not None:
            if os.path.exists(inputfile):
                self.read(inputfile)

    def read(self, filename):
        if not os.path.isfile(filename):
            raise ValueError('ERROR: File not found: %s' % filename)

        try:
            rf = open(filename)
            data = rf.readlines()
        except UnicodeDecodeError:
            rf = open(filename, encoding='ISO-8859-1')
            data = rf.readlines()

        inblock = False
        ret = {}

        for i in range(len(data)):
#            print(data[i])
            # Ignoring blanks
            if data[i].strip() == '':
                continue
            # Storing a comment line
            elif data[i].strip().startswith('#'):
                comment = data[i].strip()
            # Reading a block
            elif data[i].strip().lower().startswith('%block'):
                blockname = data[i][6:].strip()
                inblock = True
                blockdata = []
            elif data[i].strip().lower().startswith('%endblock'):
                # Store the block of data
                ret['block+'+blockname] = list(blockdata)
                inblock = False
            elif inblock:
                lblock = self.process_line(data[i])
                blockdata.append(lblock)
            else:
                lblock = self.process_line(data[i])
                #print('==> %s' % lblock)
                varname = lblock[0]
                if '_' in varname:
                    varname = varname.replace('_', '')
                if len(lblock) == 2 and isinstance(lblock[1], str):
                    if lblock[1].lower() in ['t', 'true', '.true.', 'yes']:
                        value = True
                    elif lblock[1].lower() in ['f', 'false', '.false.', 'no']:
                        value = False
                    else:
                        value = lblock[1]
                    ret[varname] = value
                elif len(lblock) == 2:
                    ret[varname] = lblock[1]
                elif len(lblock) == 1:
                    ret[varname] = True
                elif len(lblock) > 2 and all(isinstance(x, str) for x in lblock):
                    ret[varname] = data[i].strip()[len(varname):].strip()
                else:
                    ret[varname] = lblock[1:]
        self.variables = ret

    def process_line(self, line):
        # Clean from comments
        if '#' in line:
            data = line.split('#')[0]
        else:
            data = line
        # Line block of data
        lblock = []
        for jtoken in data.split():
            try:
                token = int(jtoken)
                lblock.append(token)
            except ValueError:
                try: 
                    token = float(jtoken)
                    lblock.append(token)
                except ValueError:
                    lblock.append(jtoken)
        return lblock

    def __str__(self):
        ret = ''
        
        if len(self.variables) == 0:
            return ret

        keys = list(self.variables.keys())
        maxlength = max([len(i) for i in keys if not i.startswith('block+')])

        partial_keys = [x for x in keys if x.lower().startswith('system')]
        if len(partial_keys) > 0:
            ret += '# System Variables\n\n'
            # System variables
            for i in partial_keys:
                ret += self.write_one(i, maxlength)
                keys.remove(i)
            ret += '\n'

        partial_keys = [x for x in keys if x.lower().startswith('number')]
        if len(partial_keys) > 0:
            ret += '# Number Variables\n\n'
            # Number variables
            for i in partial_keys:
                ret += self.write_one(i, maxlength)
                keys.remove(i)
            ret += '\n'

        partial_keys = [x for x in keys if x.lower().startswith('save')]
        if len(partial_keys) > 0:
            ret += '# Save Variables\n\n'
            # Save variables
            for i in [x for x in keys if x.lower().startswith('save')]:
                ret += self.write_one(i, maxlength)
                keys.remove(i)
            ret += '\n'

        partial_keys = [x for x in keys if not x.lower().startswith('block+')]
        if len(partial_keys) > 0:
            ret += '# Other variables\n\n'
            # Other variables
            for i in partial_keys:
                ret += self.write_one(i, maxlength)
                keys.remove(i)
            ret += '\n'

        partial_keys = [x for x in keys if x.lower().startswith('block+')]
        if len(partial_keys) > 0:
            ret += '# Block Sections\n\n'
            # Block Sections
            for i in partial_keys:
                ret += self.write_one(i, maxlength) + '\n'
                keys.remove(i)
            ret += '\n'

        return ret

    def write_one(self, key, maxlength):

        if isinstance(self.variables[key], str):
            ret = "%s %s" % (key.ljust(maxlength), self.variables[key])
        elif isinstance(self.variables[key], bool):
            if self.variables[key]:
                ret = "%s %s" % (key.ljust(maxlength), 'T')
            else:
                ret = "%s %s" % (key.ljust(maxlength), 'F')
        elif isinstance(self.variables[key], int):
            ret = "%s %d" % (key.ljust(maxlength), self.variables[key])
        elif isinstance(self.variables[key], float):
            ret = "%s %f" % (key.ljust(maxlength), self.variables[key])
        elif key.startswith('block+'):
            ret = '%' + ('block %s\n' % key[6:])
            for iline in self.variables[key]:
                for jtoken in iline:
                    if isinstance(jtoken, float):
                        ret += " %9.5f" % jtoken
                    elif isinstance(jtoken, int):
                        ret += " %d" % jtoken
                    else:
                        ret += " %s" % jtoken
                ret += '\n'
            ret += '%' + ('endblock %s' % key[6:])
        elif isinstance(self.variables[key], list):
            ret = "%s" % key.ljust(maxlength)
            for i in self.variables[key]:
                ret += " %s" % i
        else:
            raise ValueError("%s %s" % (key, self.variables[key]))
        return ret + '\n'

    def get_structure(self):

        natoms = self.get('NumberOfAtoms')
        nspecies = self.get('NumberOfSpecies')
        chemspecies = self.get('block+ChemicalSpeciesLabel')

        #get lattice constant
        lattice_constant = self.get('LatticeConstant')[0]

        #dot with lattice constant to scale
        cell = np.array(self.get('block+LatticeVectors'),dtype=float)*float(lattice_constant)

        positions =[] #np.zeros(natoms)
        symbols = natoms * [None]

        if self.has_variable('AtomicCoordinatesFormat'):
            atomic_format = self.get('AtomicCoordinatesFormat')
        else:
            atomic_format = 'NotScaledCartesianBohr'

        atomic_coords = self.get('block+AtomicCoordinatesAndAtomicSpecies')

        # index = 0
        # for iline in atomic_coords:
        #     positions[index] = np.array(iline[:3])
        #     for i in chemspecies:
        #         if i[0] == iline[-1]:
        #             symbol = atomic_symbol(i[1])
        #     symbols[index] = symbol
        # index += 1

        #modified a bit
        index = 0
        for iline in atomic_coords:
            positions.append(np.array(iline[:3])) #since we can't set an array element with a sequence in np
            symbols[index] = iline[-1]
            index +=1 

        if atomic_format in ['Fractional', 'ScaledByLatticeVectors']:
            return Structure(symbols=symbols, reduced=positions, cell=cell, periodicity=True)
        elif atomic_format in ['Bohr', 'NotScaledCartesianBohr']:
            return Structure(symbols=symbols, positions=positions, cell=cell, periodicity=True)
        elif atomic_format in ['Ang', 'NotScaledCartesianAng']:
            return Structure(symbols=symbols, positions=positions, cell=cell, periodicity=True)
