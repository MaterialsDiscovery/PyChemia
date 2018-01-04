import re
import os
import numpy as np
from math import sqrt
from pychemia.utils.netcdf import file2dict
from ..codes import CodeOutput


class AbinitOutput(CodeOutput):

    def __init__(self, filename='abinit.out'):

        CodeOutput.__init__(self)
        self.filename = None
        self.data = ''
        if os.path.isfile(filename):
            self.filename = filename
        elif os.path.isdir(filename) and os.path.isfile(filename + os.sep + 'abinit.out'):
            self.filename = filename + os.sep + 'abinit.out'
        if self.filename is not None:
            rf = open(self.filename)
            self.data = rf.read()
            if "Calculation completed." in self.data:
                self.read()

    def reload(self):
        rf = open(self.filename)
        self.data = rf.read()

    def get_energetics(self):

        ret = re.findall('ETOT\s*([\d]+)\s*([.E\d\-+]+)\s*([.E\d\-+]+)\s*([.E\d\-+]+)\s*([.E\d\-+]+)',
                         self.data)
        if not ret:
            raise RuntimeError("Could not extract energetic lines from: %s" % self.filename)

        newret = []
        for y in ret:
            try:
                etotline = [float(x) for x in y]
                newret.append(etotline)
            except ValueError:
                print("This line could not be parsed correctly: %s" % y)

        ret = np.array(newret)
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

    def get_occupation_matrices(self):

        occ_block = re.findall(r"=== For Atom[\s\w\d\-.,=>:]*\n", self.data)
        if not occ_block:
            raise RuntimeError("ERROR: Could not get occupation matrices blocks: %s" % self.filename)
        if len(occ_block) != 1:
            raise RuntimeError("ERROR: Occupation matrices should be in one block from the output")

        occs = re.findall('Occupation matrix for spin[\s\d\w]*([\s\d\-.]*)', occ_block[0])
        atoms = [int(x) for x in re.findall('For Atom([\s\w\d]*)', occ_block[0])]
        spins = [int(x) for x in re.findall('Occupation matrix for spin([\s\w\d]*\n)', occ_block[0])]

        # print(occs)
        # print(atoms)
        # print(spins)

        # print("Number of atoms: %d" % len(atoms))
        # print("Number of spins blocks: %d" % len(spins))
        # print("Number of occ matrices: %d" % len(occs))

        if 2 * len(atoms) != len(spins) or len(spins) != len(occs):
            raise RuntimeError("Number of matrices is not consistent with number of atoms and spin")

        ret = {}
        for i in range(len(occs)):
            atom = atoms[int((i + 1) / 2) + (i + 1) % 2 - 1]
            spin = spins[i]
            occ_matrix = [float(x) for x in occs[i].strip().split()]
            if len(occ_matrix) not in [25, 49]:
                raise RuntimeError("Occupation matrices must be 5x5 or 7x7")
            ret[(atom, spin)] = np.array(occ_matrix).reshape((-1, int(sqrt(len(occ_matrix))),
                                                              int(sqrt(len(occ_matrix)))))
        return ret

    def get_dmatpawu(self, both_spins=True):
        om = self.get_occupation_matrices()
        min_atom = min([x[0] for x in om.keys()])
        max_atom = max([x[0] for x in om.keys()])
        min_spin = min([x[1] for x in om.keys()])
        max_spin = max([x[1] for x in om.keys()])
        if not both_spins:
            max_spin = min_spin
        ret = None
        for i in range(min_atom, max_atom + 1):
            for j in range(min_spin, max_spin + 1):
                if ret is None:
                    ret = om[(i, j)]
                else:
                    ret = np.vstack((ret, om[(i, j)]))
        return ret

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
                         "Atom\s*Radius\s*up_density\s*dn_density\s*Total\(up\+dn\)\s*Diff\(up-dn\)\s*"
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

# DEPRECATED CODE, all ABINIT specific operations from orbitaldftu
# is moved here. 
#

# def get_final_correlation_matrices_from_output(filename):
#     rf = open(filename)
#     data = rf.read()
#     mainblock = re.findall('LDA\+U DATA[\s\w\d\-.=,>:]*\n', data)
#     if len(mainblock) != 1:
#         for isection in mainblock:
#             print(isection)
#         raise ValueError("Occupations for correlated orbitals could not be detected right")


#     pattern = "For Atom\s*(\d+), occupations for correlated orbitals. lpawu =\s*([\d]+)\s*Atom\s*[\d]+\s*. Occ. " \
#               "for lpawu and for spin\s*\d+\s*=\s*([\d\.]+)\s*Atom\s*[\d]+\s*. Occ. for lpawu and for " \
#               "spin\s*\d+\s*=\s*([\d\.]+)\s*=> On atom\s*\d+\s*,  local Mag. for lpawu is[\s\d\w\.\-]*== Occupation " \
#               "matrix for correlated orbitals:\s*Occupation matrix for spin  1\s*([\d\.\-\s]*)Occupation matrix " \
#               "for spin  2\s*([\d\.\-\s]*)"
#     ans = re.findall(pattern, mainblock[0])
#     # print(ans)

#     ret = []
#     for i in ans:
#         atom_data = {'atom number': int(i[0]),
#                      'orbital': int(i[1]),
#                      'occ spin 1': float(i[2]),
#                      'occ spin 2': float(i[3])}
#         matrix = [float(x) for x in i[4].split()]
#         atom_data['matrix spin 1'] = list(matrix)
#         matrix = [float(x) for x in i[5].split()]
#         atom_data['matrix spin 2'] = list(matrix)
#         ret.append(atom_data)
#     return ret


# def get_final_dmatpawu(filename):
#     ret = get_final_correlation_matrices_from_output(filename)
#     # Looking for input file
#     dirname = os.path.dirname(os.path.abspath(filename))
#     outnc = [ x for x in os.listdir(dirname) if x[-6:]=='OUT.nc']
#     if len(outnc) < 1:
#         raise ValueError("ERROR: Could not find *.OUT.nc needed to determine nspol, nspden and nspinor")
#     variables = netcdf2dict(dirname + os.sep + outnc[0])

#     if 'nspinor' in variables:
#         nspinor = variables['nspinor']
#     else:
#         nspinor = 1

#     if 'nsppol' in variables:
#         nsppol =  variables['nsppol']
#     else:
#         nsppol = 1

#     if 'nspden' in variables:
#         nspden = variables['nspden']
#     else:
#         nspden = nsppol

#     # number of matrices per atom:
#     nmpa = get_num_matrices_per_atom(nsppol, nspden, nspinor)

#     dmatpawu = []
#     for i in ret:
#         dmatpawu += i['matrix spin 1']
#         if nmpa == 2:
#             dmatpawu += i['matrix spin 2']        
#     return dmatpawu
