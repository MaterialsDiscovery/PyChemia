import os
import re
import subprocess
import numpy as np
import numbers
from ..codes import CodeRun
from pychemia.crystal import KPoints
from pychemia import Structure, pcm_log
from pychemia.utils.serializer import generic_serializer


class DFTBplus(CodeRun):
    def __init__(self, workdir='.'):

        if not os.path.lexists(workdir):
            os.mkdir(workdir)

        CodeRun.__init__(self, executable='dftb+', workdir=workdir)
        self.geometry = {}
        self.driver = {}
        self.hamiltonian = {}
        self.options = {}
        self.analysis = {}
        self.parser_options = {}
        self.slater_koster = []
        self.structure = Structure()
        self.runner = None
        self.kpoints = None
        self.stdout_file = None
        self.output = None
        self.kp_density = None

    def initialize(self, structure, kpoints=None, kp_density=10000):
        assert structure.is_crystal
        assert structure.is_perfect
        self.structure = structure
        self.get_geometry()
        self.kp_density = kp_density
        if kpoints is None:
            kpoints = KPoints.optimized_grid(self.structure.lattice, kp_density=self.kp_density, force_odd=True)
        self.kpoints = kpoints

    def set_inputs(self):
        self.write_input(filename=self.workdir + os.sep + 'dftb_in.hsd')
        self.link_slater_koster()

    def get_outputs(self):
        self.output = read_detailed_out(filename=self.workdir + os.sep + 'detailed.out')

    def run_status(self):
        if self.runner is None:
            pcm_log.info('DFTB+ not finish')
            filename = self.workdir + os.sep + 'dftb_stdout.log'
            if os.path.exists(filename):
                booleans, geom_optimization, stats = read_dftb_stdout(filename=filename)
                pcm_log.debug(str(booleans))
                pcm_log.debug(str(geom_optimization))
                pcm_log.debug(str(stats))
            return
        if self.runner.poll() == 0:
            pcm_log.info('DFTB+ complete normally')

    def finalize(self):
        self.stdout_file.close()

    @staticmethod
    def read_input(filepath='dftb_in.hsd'):
        """
        Read an Input for DFTB+
        Support a very restricted subset of
        Human-readable Structured Data (HSD)

        :param filepath: (str) Filename to read as DFTB+ input
        :return: (dict)
        """
        rf = open(filepath, 'r')
        ret = []
        gen_format = False
        value = None
        level = 0
        current_container = ret
        tree_pos = [ret]
        for line in rf.readlines():
            line = line.partition('#')[0]
            line = line.rstrip()
            print('line--->', line)
            print('return->', ret)
            # print 'container',current_container
            # print 'tree_pos',tree_pos
            print(80 * '-')
            if line.strip() == '':
                continue
            elif '=' in line and line.strip()[-1] == '{':
                leftside = line.split('=')[0].strip()
                rightside = line.split('=')[1].split()[0].strip()
                ret.append({'_left': leftside, '_rigth': rightside})
                current_container = ret[-1]
                tree_pos.append(current_container)
                if rightside == 'GenFormat':
                    gen_format = True
                level += 1
                print('container', current_container)
            elif gen_format and line.strip() != '}':
                if value is None:
                    value = line
                else:
                    value += line
            elif line.strip() == '}':
                if gen_format:
                    current_container['value'] = value
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
                    ret.append({line.split()[0]: '{}'})
                else:
                    raise ValueError('Line not parsed correctly: ' + line)
            else:
                raise ValueError('Line not parsed correctly: ' + line)
        return ret

    def get_geometry(self):
        assert (self.structure is not None)

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

    def set_static(self):
        self.driver = {}

    def basic_options(self):
        self.options = {'WriteDetailedXML': True}

    def basic_analysis(self):
        self.analysis = {'ProjectStates': []}

    def basic_parser_options(self):
        self.parser_options = {'WriteXMLInput': True}

    @staticmethod
    def print_block(name, block):

        if block == {}:
            return '\n' + name + ' = {}\n'

        ret = '\n' + name + ' = '
        if 'method' in block:
            ret += block['method'] + ' {\n'
        else:
            ret += ' {\n'
        for ivar in sorted([x for x in block if x != 'method']):
            if isinstance(block[ivar], str):
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
                    self.hamiltonian['SlaterKosterFiles'][i + '-' + j] = '"' + i + '-' + j + '.skf' + '"'

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
            if 'value' in the_dict and isinstance(the_dict['value'], np.ndarray):
                ret += DFTBplus._write_ndarray(the_dict['value'])
        elif 'units' in the_dict:
            ret += '[' + the_dict['units'] + '] = {'
            if 'coordinates' in the_dict:
                if isinstance(the_dict['coordinates'], np.ndarray):
                    ret += DFTBplus._write_ndarray(the_dict['coordinates'])
            elif 'types' in the_dict:
                ret += '\n'
                for i in range(len(the_dict['types'])):
                    ret += " %3d %10.7f %10.7f %10.7f\n" % (the_dict['types'][i],
                                                            the_dict['positions'][i][0],
                                                            the_dict['positions'][i][1],
                                                            the_dict['positions'][i][2])
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
        if not self.slater_koster:
            pcm_log.debug('The Slater-Koster files were not selected')
            pcm_log.debug('The Hamiltonian could not be setup')
        else:
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

    def write_input(self, filename='dftb_in.hsd'):

        wf = open(filename, 'w')
        wf.write(self.print_all())
        wf.close()

    def set_slater_koster(self, search_paths):
        ret = True
        if isinstance(search_paths, str):
            search_paths = [search_paths]
        elif not isinstance(search_paths, list):
            raise ValueError('search_path is not an string or list')
        self.slater_koster = []
        for ispecie in self.structure.species:
            for jspecie in self.structure.species:
                pair = ispecie + '-' + jspecie
                pair_found = False
                path = self.workdir + os.sep + pair + '.skf'
                if os.path.exists(path):
                    self.slater_koster.append(path)
                    pair_found = True
                else:
                    for ipath in search_paths:
                        path = ipath + os.sep + pair + '.skf'
                        if os.path.exists(ipath + os.sep + ispecie + '-' + jspecie + '.skf'):
                            self.slater_koster.append(path)
                            pair_found = True
                if pair_found:
                    # print 'Slater-Koster %7s : %s' % (pair, path)
                    pcm_log.debug('Slater-Koster %7s : %s' % (pair, path))
                else:
                    pcm_log.debug('ERROR: Slater_Koster for ' + pair + ' not found')
                    ret = False
        return ret

    def link_slater_koster(self):
        """
        Create symbolic links to Slater-Koster files not already on the workdir
        """
        for i in self.slater_koster:
            filename = self.workdir + os.sep + os.path.basename(i)
            if os.path.exists(filename):
                pass
            elif os.path.lexists(filename):
                os.remove(filename)
                os.symlink(os.path.abspath(i), filename)
            else:
                os.symlink(os.path.abspath(i), filename)

    @staticmethod
    def get_shells(slater_koster):

        ret = {}
        rf = open(slater_koster, 'r')
        xmldata = ''
        isxml = False
        for iline in rf.readlines():
            if iline.strip() == '<Documentation>':
                isxml = True
            elif iline.strip() == '</Documentation>':
                isxml = False
                xmldata += iline
            if isxml:
                xmldata += iline

        import xml.etree.ElementTree
        root = xml.etree.ElementTree.fromstring(xmldata)
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
        ret = {}
        for islater in self.slater_koster:
            shells = self.get_shells(islater)
            for jspecie in shells:
                if jspecie in ret:
                    if not shells[jspecie] in ret[jspecie]:
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
            ret[0:3, 0:3] = self.kpoints.kpoints_list
            ret[:, 3] = self.kpoints.weights
        elif self.kpoints.kmode in ['gamma', 'monkhorst-pack']:
            ret = {'name': 'SupercellFolding', 'value': np.zeros((4, 3))}
            ret['value'][:3, :3] = np.diag(self.kpoints.grid)
            if self.kpoints.kmode == 'gamma':
                ret['value'][3, :] = np.zeros(3)
            elif self.kpoints.kmode == 'monkhorst-pack':
                for i in range(3):
                    # Test if grid is odd
                    if self.kpoints.grid[i] & 1:
                        # ODD
                        ret['value'][3, i] = 0
                    else:
                        # EVEN
                        ret['value'][3, i] = 0.5
        return ret

    def roll_outputs(self, value):
        for ifile in ['detailed.out', 'detailed.xml', 'dftb_stdout.log', 'geo_end.gen', 'geo_end.xyz', 'dftb_in.hsd']:
            if os.path.exists(self.workdir + os.sep + ifile):
                os.rename(self.workdir + os.sep + ifile, self.workdir + os.sep + ('%03d_%s' % (value, ifile)))


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
        symbols.append(species[int(line.split()[1]) - 1])
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


def read_detailed_out(filename='detailed.out'):
    """
    Read a 'detailed.out' file and extract the values
    of forces, stress and total energy
    If those any of those three properties are not present
    its value will be set to None

    :param filename: The filename in the format of DFTB+ 'detailed.out'
    :return: tuple (forces, stress, total_energy)
    """
    ret = {}

    rf = open(filename, 'r')
    data = rf.read()

    # Extracting the data from 'detailed.out'
    forces = re.findall(r'Total\s*Forces\s*([.\-+E\d\s]+)\n', data)
    stress = re.findall(r'Total\s*stress\s*tensor\s*([.\-+E\d\s]+)\n', data)
    total_energy = re.findall(r'Total\s*energy:\s*[.\-+E\d\s]*\s*H\s*([.\-\d\sE]+)\s*eV', data)

    max_force = re.findall(r'Max force for moved atoms::\s*([.\-+E\d\s]+)\s*au', data)
    if len(max_force) > 0:
        max_force = float(max_force[-1])
    else:
        max_force = None

    max_deriv = re.findall(r'Maximal derivative component:\s*([.\-+E\d\s]+)\s*au', data)
    if len(max_deriv) > 0:
        max_deriv = float(max_deriv[-1])
    else:
        max_deriv = None

    if forces:
        forces = np.array(forces[0].split(), dtype=float).reshape((-1, 3))
        # pcm_log.debug('Forces :\n'+str(forces))
    else:
        forces = None

    if stress:
        stress = np.array(stress[0].split(), dtype=float).reshape((3, 3))
        # pcm_log.debug('Stress :\n' + str(stress))
    else:
        stress = None

    if total_energy:
        total_energy = float(total_energy[0])
        # pcm_log.debug('Energy :' + str(total_energy))
    else:
        total_energy = None

    scc = re.findall('iSCC Total electronic\s+Diff electronic\s+SCC error\s+([\s\dE+-.]*)', data)
    if len(scc) > 0:
        ret['SCC'] = {}
        ret['SCC']['iSCC'] = int(scc[0].split()[0])
        ret['SCC']['Total_electronic'] = float(scc[0].split()[1])
        ret['SCC']['Diff_electronic'] = float(scc[0].split()[2])
        ret['SCC']['SCC_error'] = float(scc[0].split()[3])

    ret['forces'] = forces
    ret['stress'] = stress
    ret['total_energy'] = total_energy
    if max_force is not None:
        ret['max_force'] = max_force
    if max_deriv is not None:
        ret['max_deriv'] = max_deriv
    return ret


def read_dftb_stdout(filename='dftb_stdout.log'):
    """
    Read the standard output stored in a file
    (by default 'dftb_stdout.log') and extract
    relevant information useful for a relaxation
    procedure

    :param filename: The standard output stored in a file
    :return:
    """

    rf = open(filename, 'r')
    data = rf.read()

    ret = {}

    start_ini = re.findall(r'Starting initialization...[\s-]+([.\s\d\w\-:(),+]+)---', data)
    if len(start_ini) == 1:
        start_ini = start_ini[0].replace('\n  ', '  ').split('\n')
        ret['Starting_Initialization'] = {}
        for line in start_ini:
            if len(line.split(':')) == 2:
                left, right = line.split(':')
            elif line[:3] != '---':
                left = line[:29].replace(':', '')
                right = line[29:]
            ret['Starting_Initialization'][left.strip()] = parse_value(right.strip())

            if left.strip() == 'K-points and weights':
                right = right.replace(':', '')
                right = np.array(right.split(), dtype=float).reshape((-1, 5))
                ret['Starting_Initialization'][left.strip()] = {'kpoints': generic_serializer(right[:, 1:-1]),
                                                                'weights': generic_serializer(right[:, -1])}

    geom_blocks = re.findall(r'Geometry step:([\^*.\s\d\w\-:(),+]+)Pa', data)

    ret['Geometry_Steps'] = []
    for iblock in geom_blocks:
        ib = iblock.split('\n')
        if 'Lattice' in ib[0]:
            tmp = {'iter': int(ib[0].split(',')[0]), 'scc': {}}
        else:
            tmp = {'iter': int(ib[0]), 'scc': {}}
        tmp['scc']['tot_electronic'] = []
        tmp['scc']['diff_electronic'] = []
        tmp['scc']['scc_error'] = []
        for iline in ib[4:]:
            fields = iline.split()
            if len(fields) > 0 and fields[0].strip().isdigit():
                tmp['scc']['tot_electronic'].append(float(fields[1]))
                tmp['scc']['diff_electronic'].append(float(fields[2]))
                tmp['scc']['scc_error'].append(float(fields[3]))
            elif len(iline.split(':')) == 2:
                left, right = iline.split(':')
                tmp[left.strip()] = {'value': float(right.split()[0]), 'units': right.split()[1]}

        ret['Geometry_Steps'].append(tmp)

    max_force = re.findall(r'Maximal force component:\s*([\s\dE+-.]+)\s*\n', data)
    if len(max_force) > 0:
        ret['max_force'] = float(max_force[-1])

    if re.findall('Geometry did NOT converge', data):
        pcm_log.debug('Convergence not achieved!')
        ret['ion_convergence'] = False
    else:
        ret['ion_convergence'] = True

    if re.findall('Geometry converged', data):
        pcm_log.debug('Convergence achieved!')
        ret['ion_convergence'] = True
    else:
        ret['ion_convergence'] = False
    return ret


def parse_value(value):
    try:
        ret = int(value)
    except ValueError:
        try:
            ret = float(value)
        except ValueError:
            if value == 'Yes':
                ret = True
            elif value == 'No':
                ret = False
            else:
                ret = value
    return ret
