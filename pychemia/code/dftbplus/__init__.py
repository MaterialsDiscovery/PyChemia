__author__ = 'Guillermo Avendano Franco'

import os

from pychemia.code import Codes


class DFTBplus(Codes):

    def __init__(self):
        self.workingdirectory = None
        self.geometry = {}
        self.driver = {}
        self.hamiltonian = {}
        self.options = {}
        self.analysis = {}
        self.parser_options = {}
        self.slater_koster = {}

    def initialize(self, dirpath):
        self.workingdirectory = dirpath
        if not os.path.lexists(dirpath):
            os.mkdir(dirpath)

    def set_inputs(self):
        pass

    def get_ouputs(self):
        pass

    def run(self):
        pass

    def finalize(self):
        pass

    @property
    def dirpath(self):
        return self.workingdirectory

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
            elif (not '=' in line) and (line.strip()[-1] == '{'):
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

    def get_geometry(self, structure):
        self.geometry['TypeNames'] = list(structure.species)
        self.geometry['TypesAndCoordinates']={'units':'Angstrom', 'coordinates':[]}
        for i in structure
