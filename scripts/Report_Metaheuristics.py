#!/usr/bin/env python

import pychemia
from pychemia.visual.searcher import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


if len(sys.argv) != 3:
    print('Use:')
    print("Report_Metaheuristics.py dbname formula")
    exit(1)

gen_size = 16
name = sys.argv[1]
formula = sys.argv[2]
navg_energy = 6
letterpage = (11, 8.5)
host = "localhost"
user = "guilleaf"
passwd = "zxcvbnm"
ssl = True


if __name__ == '__main__':
    pdf = PdfPages(name + 'multipage.pdf')

    pcdb = pychemia.db.PyChemiaDB(name, host=host, user=user, passwd=passwd, ssl=ssl)
    popu = pychemia.population.RelaxStructures(pcdb)
    print(popu)

    tags = get_all_status(popu)
    for i in ['lock', 'target_forces', 'relaxation']:
        if i in tags:
            tags.remove(i)
    print(tags)

    colors = ['r', 'g', 'm', 'c', 'y', 'k']
    color_tags = {}
    for i in tags:
        if len(colors) == 0:
            colors = ['r', 'g', 'm', 'c', 'y', 'k']
        color_tags[i] = colors.pop(0)

    for itag in tags:
        pcdb = pychemia.db.PyChemiaDB(name, host=host, user=user, passwd=passwd, ssl=ssl)
        popu = pychemia.population.RelaxStructures(pcdb, tag=itag)
        patches = []

        plot_generation_chart(popu, gen_size)

    # Change of energy vs distance
    print('Plot: Change of energy vs distance')
    ret = {}
    color = {'replace_by_other': 'r', 'replace_by_random': 'b'}
    plt.figure(figsize=letterpage)
    for change in ['replace_by_other', 'replace_by_random']:
        print('Searching for: ', change)
        ret[change] = []
        for entry in popu.pcdb.db.generation_changes.find({'change': change}):
            id_from = entry['from']
            id_to = entry['to']
            # print "%s --> %s" % (id_from, id_to)
            x = popu.distance(id_from, id_to)
            if popu.is_evaluated(id_from) and popu.is_evaluated(id_to):
                values = popu.get_values([id_from, id_to])
                delta_energy = values[id_to] - values[id_from]
                y = delta_energy
                ret[change].append([x, y])
        ret[change] = np.array(ret[change])
        if len(ret[change]) > 0:
            plt.plot(ret[change][:, 0], ret[change][:, 1], color[change] + 'o', label=change)

    plt.xlim(min(ret[change][:, 0]) - 0.01, 0.01 + max(ret[change][:, 0]))
    plt.xlabel('Distance between structures (after relaxation)')
    plt.ylabel('Energy change')
    plt.title('Change of energy due to movement of structures')
    plt.legend()
    pdf.savefig()
    plt.close()

    # Getting all known Structures
    pcdb = pychemia.db.PyChemiaDB('ICSD_Pures', host='157.182.27.249', user=None, passwd=None, ssl=False)
    known_structures = {'energy_pa': [], 'spacegroup': []}
    number_known = 0
    for entry in pcdb.entries.find({'structure.formula': formula,
                                    'properties.energy_pa': {'$exists': 1}}, {'properties': 1}):
        known_structures['energy_pa'].append(entry['properties']['energy_pa'])
        known_structures['spacegroup'].append(entry['properties']['spacegroup'])
        number_known += 1

    # Best energies vs Generations
    miny = min(known_structures['energy_pa'])
    maxy = max(known_structures['energy_pa'])
    rangey = (miny - 0.1 * (maxy - miny), maxy + 0.1 * (maxy - miny))

    plt.figure(figsize=letterpage)

    maxx = 0
    for itag in tags:

        pcdb = pychemia.db.PyChemiaDB(name, host=host, user=user, passwd=passwd, ssl=ssl)
        popu = pychemia.population.RelaxStructures(pcdb, tag=itag)

        inigen = fingen = 0
        inigen, fingen = get_generation_limits(popu, gen_size)

        if fingen > maxx:
            maxx = fingen

        # All energies vs Generations
        sys.stdout.write("Plot: All energies vs Generations:")
        for j in range(inigen, fingen):
            sys.stdout.write(' %d' % j)
            thegen = get_generation(popu, j)
            energies = []
            for i in range(len(thegen)):
                struct, properties, status = popu.pcdb.get_dicts(thegen[i])
                energy = properties['energy_pa']
                energies.append(energy)
            if len(energies) > 0:
                for i in range(len(energies)):
                    plt.plot([j, j + 1 - 0.1], [energies[i], energies[i]], color_tags[itag] + '-')
                    plt.text(j + 0.9 * float(i) / len(energies), energies[i], str(i),
                             ha="center", family='sans-serif', size=6)
        sys.stdout.write('\n')

    for i in range(number_known):
        energy_pa = known_structures['energy_pa'][i]
        spacegroup = known_structures['spacegroup'][i]
        plt.plot([inigen, maxx], energy_pa * np.ones(2), 'b--')
        plt.text(fingen, energy_pa, str(spacegroup), ha="left", family='sans-serif', size=6, color='b')

    plt.xticks(range(maxx))
    plt.ylabel('Energy per atom [eV]')
    plt.xlabel('Generations')
    plt.title('Energies')
    pdf.savefig()
    plt.close()

    plt.figure(figsize=letterpage)
    maxx = 0
    for itag in tags:

        pcdb = pychemia.db.PyChemiaDB(name, host=host, user=user, passwd=passwd, ssl=ssl)
        popu = pychemia.population.RelaxStructures(pcdb, tag=itag)

        inigen, fingen = get_generation_limits(popu, gen_size)
        if fingen > maxx:
            maxx = fingen

        sys.stdout.write("Plot: Best energies vs Generations:")
        best = 16
        for j in range(inigen, fingen):
            sys.stdout.write(' %d' % j)
            thegen = get_generation(popu, j)[:best]
            energies = []
            spacegroups = []
            for i in range(len(thegen)):
                struct, properties, status = popu.pcdb.get_dicts(thegen[i])
                energies.append(properties['energy_pa'])
                spacegroups.append(properties['spacegroup'])
            if len(energies) > 0:
                for i in range(len(energies)):
                    if rangey[0] < energies[i] < rangey[1]:
                        plt.plot([j, j + 1 - 0.1], [energies[i], energies[i]], color_tags[itag] + '-')
                        plt.text(j + 0.9 * float(i) / len(energies), energies[i], str(i),
                                 ha="center", family='sans-serif', size=6, color='0.9')
                        plt.text(j + 0.9, energies[i], str(spacegroups[i]),
                                 ha="center", family='sans-serif', size=6, color=color_tags[itag])
        sys.stdout.write('\n')

    for i in range(number_known):
        energy_pa = known_structures['energy_pa'][i]
        spacegroup = known_structures['spacegroup'][i]
        plt.plot([inigen, maxx], energy_pa * np.ones(2), 'b--')
        plt.text(maxx, energy_pa, str(spacegroup), ha="left", family='sans-serif', size=6, color='b')

    plt.ylim(*rangey)
    plt.xticks(range(maxx))
    plt.ylabel('Energy per atom [eV]')
    plt.xlabel('Generations')
    plt.title('Energies')
    pdf.savefig()

    plt.figure(figsize=letterpage)

    for i in range(number_known):
        energy_pa = known_structures['energy_pa'][i]
        spacegroup = known_structures['spacegroup'][i]
        plt.plot([0, 1], energy_pa * np.ones(2), 'b--')
        plt.text(0.5, energy_pa, str(spacegroup), ha="left", family='sans-serif', size=6, color='b')

    index = 1
    for itag in tags:

        pcdb = pychemia.db.PyChemiaDB(name, host=host, user=user, passwd=passwd, ssl=ssl)
        popu = pychemia.population.RelaxStructures(pcdb, tag=itag)
        pre_energy = 1E10

        for entry_id in popu.ids_sorted(popu.evaluated):
            entry = popu.get_entry(entry_id)

            if itag in entry['status']:
                energy_pa = entry['properties']['energy_pa']
                spacegroup = entry['properties']['spacegroup']
                if np.abs(pre_energy - energy_pa) > 1E-3 and rangey[0] < energy_pa < rangey[1]:
                    plt.plot([index, index + 1], energy_pa * np.ones(2), color_tags[itag] + '--')
                    plt.text(index + 0.5, energy_pa, str(spacegroup),
                             ha="left",
                             family='sans-serif',
                             size=6,
                             color=color_tags[itag])
                    pre_energy = energy_pa

        index += 1

    plt.ylim(*rangey)
    plt.xticks(0.5 + np.arange(len(tags) + 1), ['CIFs'] + tags)
    plt.ylabel('Energy per atom [eV]')
    plt.xlabel('Structures found')
    plt.title('Energies')
    pdf.savefig()

    plt.figure(figsize=letterpage)

    for i in range(number_known):
        energy_pa = known_structures['energy_pa'][i]
        spacegroup = known_structures['spacegroup'][i]
        plt.plot([0, 1], energy_pa * np.ones(2), 'b-')
        # plt.text(0.5, energy_pa,str(spacegroup), ha="left", family='sans-serif', size=6, color='b')

    pcdb = pychemia.db.PyChemiaDB(name, host=host, user=user, passwd=passwd, ssl=ssl)
    popu = pychemia.population.RelaxStructures(pcdb, distance_tolerance=0.4, value_tol=0.1)
    pre_energy = 1E10

    selection = popu.evaluated
    n = len(selection)
    while True:
        selection = popu.cleaned_from_duplicates(selection)
        if len(selection) == n:
            break
        n = len(selection)
        print('Selection', n)

    for entry_id in popu.ids_sorted(selection):
        entry = popu.get_entry(entry_id)
        energy_pa = entry['properties']['energy_pa']
        spacegroup = entry['properties']['spacegroup']
        if np.abs(pre_energy - energy_pa) > 1E-3 and rangey[0] < energy_pa < rangey[1]:
            plt.plot([1, 2], energy_pa * np.ones(2), 'r-')
            # plt.text(1.5, energy_pa,str(spacegroup), ha="left", family='sans-serif', size=6, color='r')
            pre_energy = energy_pa

    plt.ylim(*rangey)
    plt.xticks(0.5 + np.arange(2), ['Known Structures', 'New predicted'])
    plt.ylabel('Energy per atom [eV]')
    plt.xlabel('Structures found')
    plt.title('Energies')
    pdf.savefig()

    pdf.close()
