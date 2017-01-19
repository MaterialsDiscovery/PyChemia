import subprocess
import numpy as np
from pychemia.code.vasp import write_poscar


class PhonopyJob:

    def __init__(self, structure, dim=(2, 2, 2)):

        self.structure = structure
        self.dim = np.array(dim)

    def create_poscars(self):

        write_poscar(self.structure)
        subprocess.call(['phonopy', '-d', r'--dim="%d' % self.dim[0], '%d' % self.dim[1],  r'%d"' % self.dim[2], '-c',
                         'POSCAR'])

