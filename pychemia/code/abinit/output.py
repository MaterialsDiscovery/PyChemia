import re
import os
import numpy as np
from pychemia.utils.netcdf import file2dict
from ..codes import CodeOutput

class AbinitOutput(CodeOutput):

    def __init__(self, filename='abinit.out'):

        CodeOutput.__init__(self)
        self.filename = filename
        rf = open(self.filename)
        self.data = rf.read()
        self.read()

    def reload(self):
        rf = open(self.filename)
        self.data = rf.read()

    def get_energetics(self):

        ret = re.findall('ETOT\s*([\d]+)\s*([.E\d\-+]+)\s*([.E\d\-+]+)\s*([.E\d\-+]+)\s*([.E\d\-+]+)\n',
                         self.data)
        if not ret:
            return None
        else:
            ret = [[float(x) for x in y] for y in ret]
            ret = np.array(ret)
            energetics = {'iter': list(ret[:, 0]),
                          'etot': list(ret[:, 1]),
                          'deltaEh': list(ret[:, 2]),
                          'residm': list(ret[:, 3]),
                          'nres2': list(ret[:, 4])}
            return energetics

    def read(self):
        self.output_values = self.get_energetics()

    @property
    def is_finished(self):
        return "Calculation completed." in self.data

    def get_occupation_matrix(self):

        occ_block = re.findall(r"=== For Atom[\s\w\d\-.,=>:]*\n \n \n", self.data)
        if not occ_block:
            return None
        occs = re.findall('Occupation matrix for spin[\s\d\w]*([\s\d\-.]*)', occ_block[0])
        atoms = [int(x) for x in re.findall('For Atom([\s\w\d]*)', occ_block[0])]
        spins = [int(x) for x in re.findall('Occupation matrix for spin([\s\w\d]*\n)', occ_block[0])]

        ret = []
        n = len(occs)
        if 2 * len(atoms) == len(spins) and len(spins) == n:
            for i in range(n):
                atom = atoms[(i + 1) / 2 + (i + 1) % 2 - 1]
                spin = spins[i]
                occ_matrix = [float(x) for x in occs[i].strip().split()]
                ret.append({'atom': atom, 'spin': spin, 'occ_matrix': occ_matrix})
        return ret

    def get_dmatpawu(self):
        om = self.get_occupation_matrix()
        min_atom = min([x['atom'] for x in om])
        max_atom = max([x['atom'] for x in om])
        ret = []
        for i in range(min_atom, max_atom + 1):
            for iom in om:
                if iom['atom'] == i and iom['spin'] == 1:
                    ret.append(iom['occ_matrix'])
        return np.array(ret)

    def get_integrated_electronic_densities(self):
        """
        Return the integrated electronic densities in atomic spheres from the output
        Data is returned serialized as python dictionary with keys
        'radius', 'up_density', 'dn_density', 'total_up+dn' and 'diff_up-dn'
        If the method is unable to parse the output file None is returned.
        
        :return: (dict) Integrated electronic densities around an atomic radius, up, down, sum and diff
        """

        ref = re.findall("Integrated electronic and magnetization densities in atomic spheres:\s*"
                         "---------------------------------------------------------------------\s*"
                         "Note: Diff\(up-dn\) is a rough approximation of local magnetic moment\s*"
                         "Atom\s*Radius\s*up_density\s*dn_density\s*Total\(up\+dn\)\s*Diff\(up\-dn\)\s*"
                         "([\s\d\w.-]*)\n -+\n", self.data)
        # Return None if could not parse the magnetization block
        if len(ref) == 0:
            return None
        else:
            magmoms = np.fromstring(ref[0], dtype=float, sep=' ')
            magmoms.reshape((-1, 6))
            return {'radius': list(magmoms[:, 1]),
                    'up_density': list(magmoms[:, 2]),
                    'dn_density': list(magmoms[:, 3]),
                    'total_up+dn': list(magmoms[:, 4]),
                    'diff_up-dn': list(magmoms[:, 5])}

    @staticmethod
    def read_output_netcdf(filename):
        return file2dict(filename)
