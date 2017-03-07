import re


def read_log(filename):

    rf = open(filename)
    data = rf.read()
    energy_data = re.findall('Total Energy\s[\-]+\s([\s\d\w/\-.=]*)Forces', data)

    energetics = {}
    if len(energy_data) > 0:
        for iline in energy_data[-1].split('\n'):
            if '=' in iline:
                energetics[iline.split('=')[0].strip()] = float(iline.split('=')[1])

    forces_data = re.findall('Total Forces:\s+\s([\s\d\w/\-+.=:()]*)\n \n', data)

    forces = []
    if len(forces) > 0:
        for iline in forces[-1].split('\n'):
            if '=' in iline:
                forces.append([float(x) for x in iline.split('=')[1].split()])

    ret = {'energetics': energetics, 'forces': forces}
    return ret


def write_inp(structure, filename):

    wf = open(filename, 'w')

    if structure.is_crystal:
        periodicity_tag = 0
    else:
        periodicity_tag = 1

    wf.write("%d %d\n" % (structure.natom, periodicity_tag))


