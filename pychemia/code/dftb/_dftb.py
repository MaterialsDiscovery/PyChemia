__author__ = 'Guillermo Avendano Franco'

import os
import subprocess
import numpy as np
import numbers
from pychemia.code._codes import Codes
from pychemia.dft._kpoints import KPoints
from pychemia import Structure

class DFTBplus(Codes):

    def __init__(self):
        self.workdir = None
        self.geometry = {}
        self.driver = {}
        self.hamiltonian = {}
        self.options = {}
        self.analysis = {}
        self.parser_options = {}
        self.slater_koster = None
        self.structure = None

    def initialize(self, workdir, structure, kpoints):
        assert structure.is_crystal
        assert structure.is_perfect
        self.structure = structure
        self.get_geometry()
        self.workdir = workdir
        if not os.path.lexists(workdir):
            os.mkdir(workdir)
        self.kpoints = kpoints

    def set_inputs(self):
        self.print_input(filename=self.workdir+os.sep+'dftb_in.hsd')
        self.print_slater_koster()

    def get_ouputs(self):
        pass

    def run(self):
        os.chdir(self.workdir)
        subprocess.call(["dftb+"])

    def finalize(self):
        pass

    @staticmethod
    def read_input(filepath):
        """
        Read an Input for DFTB+
        Support a very restricted subset of
        Human-readable Structured Data (HSD)

        :return: (dict)
        """
        rf = open(filepath, 'r')
        ret = []
        gen_format = False
        value = None
        leftside = None
        rightside = None
        level = 0
        current_container = ret
        tree_pos = [ret]
        for line in rf.readlines():
            print 'line--->', line
            print 'return->', ret
            #print 'container',current_container
            #print 'tree_pos',tree_pos
            print 80*'-'
            if line.strip() == '':
                pass
            elif '=' in line and line.strip()[-1] == '{':
                leftside = line.split('=')[0].strip()
                rightside = line.split('=')[1].split()[0].strip()
                ret.append({'_left': leftside, '_rigth': rightside})
                current_container = ret[-1]
                tree_pos.append(current_container)
                if rightside == 'GenFormat':
                    gen_format = True
                level += 1
                print 'container', current_container
            elif gen_format and line.strip() != '}':
                if value is None:
                    value = line
                else:
                    value += line
            elif line.strip() == '}':
                if gen_format:
                    current_container['value'] = value
                    leftside = None
                    rightside = None
                    value = None
                    gen_format = False
                tree_pos = tree_pos[:-1]
                current_container = tree_pos[-1]
                level -= 1
            elif '=' in line and level > 0:
                current_container[line.split('=')[0].strip()] = line.split('=')[1]
            elif ('=' not in line) and (line.strip()[-1] == '{'):
                if level > 0:
                    current_container[line.split()[0].strip()] = {}
                    current_container = current_container[line.split()[0].strip()]
                    tree_pos.append(current_container)
                    level += 1
                else:
                    ret.append({line.strip().split()[0]: {}})
                    current_container = ret[-1]
                    tree_pos.append(current_container)
                    level += 1
            if line.strip()[-2:] == '{}':
                if level == 0:
                    ret. append({line.split()[0]: '{}'})
                else:
                    raise ValueError('Line not parsed correctly: '+line)
            else:
                raise ValueError('Line not parsed correctly: '+line)
        return ret

    def get_geometry(self):
        self.geometry['TypeNames'] = list(self.structure.species)
        self.geometry['TypesAndCoordinates'] = {'units': 'Angstrom'}
        indices = []
        for isite in self.structure:
            indices.append(self.geometry['TypeNames'].index(isite.symbols[0]) + 1)
        self.geometry['TypesAndCoordinates']['types'] = np.array(indices)
        self.geometry['TypesAndCoordinates']['positions'] = self.structure.positions
        self.geometry['Periodic'] = True
        self.geometry['LatticeVectors'] = {'units': 'Angstrom', 'coordinates': self.structure.cell}

    def basic_driver(self, lattice_optimization=True, maxforce=1E-2, maxsteps=50):

        self.driver['method'] = 'ConjugateGradient'
        self.driver['MaxForceComponent'] = maxforce
        self.driver['MaxSteps'] = maxsteps
        self.driver['LatticeOpt'] = lattice_optimization
        self.driver['MovedAtoms'] = '1:-1'

    def basic_options(self):
        self.options = {'WriteDetailedXML': True}

    def basic_analysis(self):
        self.analysis = {'ProjectStates': []}

    def basic_parser_options(self):
        self.parser_options = {'WriteXMLInput': True}

    def print_block(self, name, block):

        if block == {}:
            return '\n' + name + ' = {}\n'

        ret = '\n' + name + ' = '
        if 'method' in block:
            ret += block['method'] + ' {\n'
        else:
            ret += ' {\n'
        for ivar in sorted([x for x in block if x != 'method']):
            if isinstance(block[ivar], basestring):
                ret += ' ' + ivar + ' = ' + DFTBplus._write_string(block[ivar])
            elif isinstance(block[ivar], bool):
                ret += ' ' + ivar + ' = ' + DFTBplus._write_bool(block[ivar])
            elif isinstance(block[ivar], numbers.Number):
                ret += ' ' + ivar + ' = ' + DFTBplus._write_number(block[ivar])
            elif isinstance(block[ivar], dict):
                ret += ' ' + ivar + DFTBplus._write_dict(block[ivar])
            elif isinstance(block[ivar], list):
                ret += ' ' + ivar + ' = ' + DFTBplus._write_list(block[ivar])
        ret += '}\n'
        return ret

    def basic_hamiltonian(self, scc=True, slater='auto'):
        self.hamiltonian['method'] = 'DFTB'
        self.hamiltonian['SCC'] = scc
        if slater == 'auto':
            self.hamiltonian['SlaterKosterFiles'] = {}
            for i in self.structure.species:
                for j in self.structure.species:
                    self.hamiltonian['SlaterKosterFiles'][i+'-'+j] = '"' + i + '-' + j + '.skf' + '"'

        self.hamiltonian['MaxAngularMomentum'] = {}
        max_angular_momentum = self.get_max_angular_momentum()
        for i in max_angular_momentum:
            self.hamiltonian['MaxAngularMomentum'][i] = '"' + max_angular_momentum[i] + '"'

        self.hamiltonian['KPointsAndWeights'] = self.set_kpoints()

    @staticmethod
    def _write_dict(the_dict):
        ret = ' '
        if 'name' in the_dict:
            ret += ' = '
            ret += the_dict['name'] + ' {'
            if isinstance(the_dict['value'], np.ndarray):
                ret += DFTBplus._write_ndarray(the_dict['value'])
        elif 'units' in the_dict:
            ret += '[' + the_dict['units'] + '] = {'
            if 'coordinates' in the_dict:
                if isinstance(the_dict['coordinates'], np.ndarray):
                    ret += DFTBplus._write_ndarray(the_dict['coordinates'])
            elif 'types' in the_dict:
                ret += '\n'
                for i in range(len(the_dict['types'])):
                    ret += " %3d %10.7f %10.7f %10.7f\n" % ((the_dict['types'][i],) + tuple(the_dict['positions'][i]))
        else:
            ret += ' = {\n'
            for ivar in the_dict:
                ret += ' ' + ivar + ' = ' + the_dict[ivar] + '\n'
        ret += ' }\n'
        return ret

    @staticmethod
    def _write_ndarray(the_ndarray):
        ret = ''
        if the_ndarray.ndim == 1:
            ret += '\n'
            for i in the_ndarray:
                ret += ' %10.7f ' % i
        elif the_ndarray.ndim == 2:
            ret += '\n'
            for i in the_ndarray:
                for j in i:
                    ret += ' %10.7f ' % j
                ret += '\n'
        return ret

    @staticmethod
    def _write_bool(value):
        if value:
            return ' Yes\n'
        else:
            return ' No\n'

    @staticmethod
    def _write_number(value):
        return str(value) + '\n'

    @staticmethod
    def _write_list(value):
        ret = ' {'
        for i in value:
            ret += ' ' + i
        ret += ' }\n'
        return ret

    @staticmethod
    def _write_string(value):
        return value + '\n'

    def basic_input(self):
        self.basic_driver()
        self.basic_hamiltonian()
        self.basic_options()
        self.basic_analysis()
        self.basic_parser_options()

    def print_all(self):
        ret = self.print_block('Geometry', self.geometry)
        ret += self.print_block('Driver', self.driver)
        ret += self.print_block('Hamiltonian', self.hamiltonian)
        ret += self.print_block('Options', self.options)
        ret += self.print_block('Analysis', self.analysis)
        ret += self.print_block('ParserOptions', self.parser_options)
        return ret

    def print_input(self, filename='dftb_in.hsd'):

        wf = open(filename, 'w')
        wf.write(self.print_all())
        wf.close()

    def set_slater_koster(self, search_paths=[]):
        self.slater_koster = []
        for ispecie in self.structure.species:
            for jspecie in self.structure.species:
                pair = ispecie+'-'+jspecie
                pair_found =False
                for ipath in search_paths:
                    path = ipath+os.sep+pair+'.skf'
                    if os.path.exists(ipath+os.sep+ispecie+'-'+jspecie+'.skf'):
                        self.slater_koster.append(path)
                        pair_found = True
                if pair_found:
                    print 'Slater_Koster for ' + pair + ' found on ' + path
                else:
                    print 'ERROR: Slater_Koster for ' + pair + ' not found'

    def print_slater_koster(self):
        for i in self.slater_koster:
            filename = self.workdir+os.sep+os.path.basename(i)
            if not os.path.lexists(filename):
                os.symlink(i, filename)

    def get_shells(self, slater_koster):

        ret = {}
        rf = open(slater_koster, 'r')
        xmldata=''
        isxml =False
        for iline in rf.readlines():
            if iline.strip() == '<Documentation>':
                isxml = True
            elif iline.strip() == '</Documentation>':
                isxml = False
                xmldata += iline
            if isxml:
                xmldata += iline

        import xml.etree.ElementTree as ET
        root = ET.fromstring(xmldata)
        specie1 = root.findall('./General/Element1')[0].text
        specie2 = None
        if root.findall('./General/Element2'):
            specie2 = root.findall('./General/Element2')[0].text

        for ibasis in root.findall('./SK_table/Basis'):
            if ibasis.attrib['atom'] == '1':
                specie = specie1
            elif ibasis.attrib['atom'] == '2':
                specie = specie2
            else:
                specie = None
            ret[specie] = ibasis.find('./Shells').text.strip()
        return ret

    def get_all_shells(self):
        assert (self.slater_koster is not None)
        ret = {}
        for islater in self.slater_koster:
            shells = self.get_shells(islater)
            for jspecie in shells:
                if jspecie in ret:
                    if not  shells[jspecie] in ret[jspecie]:
                        ret[jspecie].append(shells[jspecie])
                else:
                    ret[jspecie] = [shells[jspecie]]
        return ret

    def get_max_angular_momentum(self):

        ret = {}
        maxang = None
        shells = self.get_all_shells()
        for ispecie in shells:
            assert (len(shells[ispecie]) == 1)
            if ispecie in self.structure.species:
                for ishell in shells[ispecie]:
                    # Take the last orbital from a shell like '2s 2p' => p
                    maxang = ishell.split()[-1][-1]
                ret[ispecie] = maxang
        return ret

    def set_kpoints(self):
        assert isinstance(self.kpoints, KPoints)
        ret = None
        if self.kpoints.kmode in ['cartesian', 'reciprocal']:
            ret = np.zeros((self.kpoints.nkpt, 4))
            ret[0:3,0:3] = self.kpoints.kpoints_list
            ret[:,3] = self.kpoints.weights
        elif self.kpoints.kmode in ['gamma', 'monkhorst-pack']:
            ret = {'name': 'SupercellFolding', 'value': np.zeros((4,3)) }
            ret['value'][:3,:3] = np.diag(self.kpoints.grid)
            if self.kpoints.kmode == 'gamma':
                ret['value'][3,:] = np.zeros(3)
            elif self.kpoints.kmode == 'monkhorst-pack':
                for i in range(3):
                    # Test if grid is odd
                    if self.kpoints.grid[i] & 1:
                        # ODD
                        ret['value'][3,i] = 0
                    else:
                        # EVEN
                        ret['value'][3,i] = 0.5
        return ret


def read_geometry_gen(filename):
    rf = open(filename, 'r')
    line = rf.readline()
    natom = int(line.split()[0])
    kind = line.split()[1].strip()
    if kind in ['S', 'F']:
        periodic = True
    elif kind == 'C':
        periodic = False
    else:
        raise ValueError('Wrong type of geometry' + kind + ". Should be ('F', 'S' or 'C')")
    species = rf.readline().split()
    coords = np.zeros((natom, 3))
    symbols = []
    for i in range(natom):
        line = rf.readline()
        coords[i, :] = np.array([float(x) for x in line.split()[2:]])
        # The indices for atoms on a gen file are in the second col [1] and start in 1
        symbols.append(species[int(line.split()[1])-1])
    if periodic:
        line = rf.readline()
        coords_origin = np.array(line.split())
        cell = np.zeros((3, 3))
        for i in range(3):
            line = rf.readline()
            cell[i, :] = np.array([float(x) for x in line.split()])

    if kind == 'F':
        return Structure(symbols=symbols, cell=cell, reduced=coords, periodicity=True)
    elif kind == 'S':
        return Structure(symbols=symbols, cell=cell, positions=coords, periodicity=True)
    elif kind == 'C':
        return Structure(symbols=symbols, positions=coords, periodicity=False)
