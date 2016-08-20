
def write_inp(structure, filename):

    wf = open(filename, 'w')

    if structure.is_crystal:
        periodicity_tag = 0
    else:
        periodicity_tag = 1

    wf.write("%d %d\n" % (structure.natom, periodicity_tag))
