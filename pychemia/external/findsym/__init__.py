import re
import subprocess


class FindSym:
    def __init__(self, structure, tolerance=0.1, executable='findsym'):

        self.structure = structure
        self.tolerance = tolerance
        self.executable = executable
        self.data = None

    def create_input(self, filename='findsym.inp'):

        wf = open(filename, 'w')

        wf.write('Input for FINDSYM, from filename\n')
        wf.write('%0.3E Tolerance\n' % self.tolerance)
        wf.write('1    form of lattice parameters: to be entered as vectors\n')
        counter = 1
        for i in self.structure.cell:
            wf.write('%12.8f %12.8f %12.8f Lattice vector %d\n' % (i[0], i[1], i[2], counter))
            counter += 1
        wf.write('P    unknown centering\n')
        wf.write('%d    number of atoms\n' % self.structure.natom)
        for i in self.structure.symbols:
            wf.write('%d ' % i)
        wf.write('Atom kinds\n')
        counter = 1
        for i in self.structure.reduced:
            wf.write('%24.16f %24.16f %24.16f Reduced Coordinates %d' % (i[0], i[1], i[2], counter))
            counter += 1

    def run(self):
        self.create_input()
        rf = open('findsym.inp')
        self.data = subprocess.check_output([self.executable], stdin=rf)
        rf.close()

    def write_cif(self, filename='structure.cif'):

        if self.data is None:
            self.run()

        findre = re.findall('# CIF file[\s\w.#,-]*', self.data)

        if len(findre) != 1:
            print('ERROR: Not CIF data was found on FindSym output')
            return
        else:
            cifdata = findre[0]

        wf = open(filename, 'w')
        wf.write(cifdata)
        wf.close()
