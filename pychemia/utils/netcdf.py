import os
import numpy as np
from scipy.io.netcdf import netcdf_file


def file2dict(filename):

    if not os.path.isfile(filename):
        raise ValueError("ERROR: Could not read %s" % filename)

    nc = netcdf_file(filename, 'r', mmap=False)
    ret = {}
    for ikey in nc.variables.keys():
        data = nc.variables[ikey].data
        if type(data[0]) == np.float64:
            if len(data) == 1:
                data = float(data[0])
            else:
                data = [float(x) for x in data]
        elif type(data[0]) == np.int32:
            if len(data) == 1:
                data = int(data[0])
            else:
                data = [int(x) for x in data]
        else:
            data = list(data)

        ret[ikey] = data
        del data

    nc.close()
    return ret


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
