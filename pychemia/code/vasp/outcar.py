import os
import re
from .output import VaspOutput


def read_vasp_stdout(filename):
    if not os.path.isfile(filename):
        raise ValueError("Could not read: %s" % filename)
    rf = open(filename)
    data = rf.read()
    rf.close()

    flaws = re.findall('[0-9]*\.[0-9]*E[+-][0-9]*-[0-9]*\.[0-9]*E[+-][0-9]*', data)
    flaws_splitted = re.findall('([0-9]*\.[0-9]*E[+-][0-9]*)(-[0-9]*\.[0-9]*E[+-][0-9]*)', data)
    for i in range(len(flaws)):
        data = data.replace(flaws[i], flaws_splitted[i][0] + ' ' + flaws_splitted[i][1])

    re_str = r'\n([\w]{3})\:\s*([\d]+)+\s*([\dE+-.]+)\s*([\dE+-.]+)\s*([\dE+-.]+)\s*([\d]+)\s*([\dE+-.]+)\s*[-\ .\w+]+'

    result = re.findall(re_str, data)

    ret = []
    for i in result:
        ret.append([i[0], int(i[1]), float(i[2]), float(i[3]), float(i[4]), int(i[5]), float(i[6])])

    number_of_scf_per_ionic_iter = []
    final_energy_after_scf = []
    index = 0
    counter = 0
    for i in ret:
        if i[1] > index:
            index = i[1]
        else:
            number_of_scf_per_ionic_iter.append(index)
            final_energy_after_scf.append(ret[counter][2])
            index = 0
        counter += 1

    return {'iterations': number_of_scf_per_ionic_iter, 'energies': final_energy_after_scf, 'data': ret}
