import re
import os
import numpy as np
from scipy.io import netcdf_file


class AbinitOutput:
    def __init__(self, filename='abinit.out'):

        self.filename = filename
        rf = open(self.filename)
        self.data = rf.read()

    def reload(self):
        rf = open(self.filename)
        self.data = rf.read()

    def get_energetics(self):

        ret = re.findall('ETOT\s*([\d]+)\s*([\.E\d\-\+]+)\s*([\.E\d\-\+]+)\s*([\.E\d\-\+]+)\s*([\.E\d\-\+]+)\n',
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

    def get_occupation_matrix(self):

        occ_block = re.findall(r"=== For Atom[\s\w\d\-\.,=>:]*\n \n \n", self.data)
        if not occ_block:
            return None
        occs = re.findall('Occupation matrix for spin[\s\d\w]*([\s\d\-\.]*)', occ_block[0])
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

    @staticmethod
    def read_output_netcdf(filename):
        return netcdf2dict(filename)


def netcdf2dict(filename):
    """
    Read a NetCDF file and create a python dictionary with
    numbers or lists for each variable

    Args:
        filename:
            NetCDF filename
    """
    if not os.path.isfile(filename):
        print('ERROR: No such file: ', filename)
        return None
    ret = {}
    netcdfile = netcdf_file(filename, 'r', mmap=False)
    for ii in netcdfile.variables.keys():
        ret[ii] = netcdfile.variables[ii][:]
    netcdfile.close()

    for i in ret:
        if ret[i].dtype == np.dtype('>f8'):
            ret[i] = [round(x, 11) for x in ret[i].flatten()]
        elif ret[i].dtype == np.dtype('>i4'):
            ret[i] = [int(x) for x in ret[i].flatten()]

    for i in ret:
        if len(ret[i]) == 1:
            ret[i] = ret[i][0]

    return ret
