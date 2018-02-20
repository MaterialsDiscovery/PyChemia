#!/usr/bin/env python

import sys
import numpy as np
import pychemia
from pychemia.visual.searcher import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from pychemia.crystal import CrystalSymmetry
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


def get_all_status(population):
    ret = []
    for entry in population.pcdb.entries.find({}, {'status': 1}):
        for x in entry['status']:
            if x not in ret:
                ret.append(x)
    return ret


def spacegroup2poly(spacegroup):
    """
    Return the size of polygon appropiated to represent a given
    crystal system

    :param spacegroup:
    :return:
    """
    if spacegroup <= 2:
        return 3
    elif spacegroup <= 15:
        return 5
    elif spacegroup <= 74:
        return 7
    elif spacegroup <= 142:
        return 9
    elif spacegroup <= 167:
        return 8
    elif spacegroup <= 194:
        return 6
    else:
        return 4


def label(ax, xy, spacegroup, energy):
    """
    Write the energy and spacegroup number on the center of each
    candidate polygon

    :param ax:
    :param xy:
    :param spacegroup:
    :param energy:
    :return:
    """
    y = xy[1] - 0.0  # shift y-value for label so that it's below the artist
    ax.text(xy[0], y, "%d" % spacegroup, ha="center", family='sans-serif', size=8)
    y = xy[1] - 0.2  # shift y-value for label so that it's below the artist
    ax.text(xy[0], y, "%.3f" % energy, ha="center", family='sans-serif', size=8)


def add_structure(ax, patches, position, spacegroup, energy):
    """
    Add one polygon to the patches list

    :param ax:
    :param patches:
    :param position:
    :param spacegroup:
    :param energy:
    :return:
    """
    polygon = mpatches.RegularPolygon(position, spacegroup2poly(spacegroup), 0.5, clip_on=True)
    patches.append(polygon)
    label(ax, position, spacegroup, energy)


def change_symbol(change):
    """
    Return the text that inform about the conditions that create that candidate

    :param change:
    :return:
    """

    if change == 'promoted':
        return 'P'
    elif change == 'modified':
        return 'M'
    elif change == 'replace_by_random':
        return 'RR'
    elif change == 'replace_by_other':
        return 'MV'
    elif change == 'duplicate':
        return '___'
    else:
        print(change)


def get_generation(searcher, tag):
    """
    Return the identifiers of all members of a population associated with
    a given searcher tag

    :param searcher:
    :param tag:
    :return:
    """
    lista = [x['_id'] for x in searcher.get_all_generations(tag) if searcher.population.is_evaluated(x['_id'])]
    return searcher.population.ids_sorted(lista)


def get_generation_limits(searcher, gen_size):
    inigen = 0
    fingen = 0
    i = 0
    while True:
        n = len(get_generation(searcher, i))
        print(' [Generation %d: Number candidates = %d]' % (i, n))
        if n < gen_size:
            fingen = i
            break
        i += 1
    return inigen, fingen


def plot_generation_chart(searcher, gen_size):
    colors = []
    patches = []
    population = searcher.population

    inigen, fingen = get_generation_limits(searcher, gen_size)

    if inigen == fingen:
        return

    fig, ax = plt.subplots(figsize=(gen_size, 2 * (fingen - inigen)))

    avg_energies = np.zeros(fingen)

    # Structures
    sys.stdout.write("Plotting Structures for generations:")
    for j in range(inigen, fingen):
        sys.stdout.write(' %d' % j)
        thegen = get_generation(population, j)[:gen_size]
        plt.text(-0.9, -2 * j, "Gen %d" % j, ha="center", family='sans-serif', size=12)
        energies = []
        for i in range(gen_size):
            struct, properties, status = population.pcdb.get_dicts(thegen[i])
            # print properties['spacegroup'], properties['energy_pa']
            spacegroup = properties['spacegroup']
            energy = properties['energy_pa']
            energies.append(energy)
            add_structure(ax, patches, [i, -2 * j], spacegroup, energy)
            colors.append(energy)
        avg_energies[j] = np.mean(energies[:navg_energy])

    collection = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.3)
    collection.set_array(np.array(colors))
    ax.add_collection(collection)
    sys.stdout.write('\n')

    # Duplicates
    sys.stdout.write("Marking duplicates for generations:")
    for j in range(inigen, fingen):
        sys.stdout.write(' %d' % j)
        thegen = get_generation(population, j)[:gen_size]
        dups = population.get_duplicates(thegen)
        for i in dups:
            x = [list(thegen).index(i), list(thegen).index(dups[i])]
            y = [-2 * j + 0.5, -2 * j + 0.5]
            plt.plot(x, y, 'k-')
    sys.stdout.write('\n')

    # Changes
    sys.stdout.write("Plotting changes for generations:")
    for j in range(inigen, fingen - 1):
        sys.stdout.write(' %d' % j)
        sys.stdout.flush()
        next_generation = get_generation(population, j + 1)
        thegen = get_generation(population, j)[:gen_size]
        for i in range(gen_size):
            # add a line
            for entry in population.pcdb.db.generation_changes.find({'from': thegen[i]}):
                csymbol = change_symbol(entry['change'])
                if csymbol is None:
                    print(entry)
                plt.text(i, -2 * j - 0.6, csymbol, ha="center", family='sans-serif', size=6)
                if 'to' in entry and entry['to'] in next_generation[:gen_size]:
                    next_id = list(next_generation).index(entry['to'])
                    x, y = np.array([[i, next_id], [-2 * j - 0.7, -2 * j - 1.3]])
                    line = mlines.Line2D(x, y, lw=3., alpha=0.3)
                    ax.add_line(line)
                elif entry['change'] == 'promoted' and entry['from'] in next_generation[:gen_size]:
                    next_id = list(next_generation[:gen_size]).index(entry['from'])
                    x, y = np.array([[i, next_id], [-2 * j - 0.7, -2 * j - 1.3]])
                    line = mlines.Line2D(x, y, lw=3., alpha=0.3)
                    ax.add_line(line)
    sys.stdout.write('\n')

    plt.title('%s (tag=%s)' % (population.name, population.tag))
    plt.text(float(gen_size) / 2, 1, '%s (tag=%s)' % (population.name, population.tag),
             ha="center",
             family='sans-serif',
             size=16)
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.axis('equal')
    plt.axis('off')
    pdf.savefig()
    plt.clf()
    plt.close()

    # Average energy of best structures
    print('Plot: Average energy of best structures')
    plt.figure(figsize=letterpage)
    plt.plot(range(fingen), avg_energies, 'ro')
    plt.plot(range(fingen), avg_energies, 'b--')
    plt.xlabel('Generation')
    plt.ylabel('Average energy of %d best structures' % navg_energy)
    plt.title('%s (tag=%s)' % (population.name, population.tag))
    plt.xticks(range(fingen))
    pdf.savefig()
    plt.close()


def plot_evolution_circular(searcher, target_function='energy_pa', tag='spacegroup'):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.cla()
    plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0, wspace=None, hspace=None)
    ax.set_xlim((-1.1, 1.1))
    ax.set_ylim((-1.1, 1.1))
    ax.set_axis_off()
    ax.set_aspect(1.0)

    for gen in range(searcher.current_generation):
        radius = 1.0 / (1.0 + np.exp(-6 * float(gen + 1) / searcher.current_generation))
        radius = (float(gen + 1) / searcher.current_generation) ** 2
        entries = {}
        counter = 0
        nele = len(searcher.lineage)
        for lin in range(nele):
            entry = searcher.population.get_entry(searcher.lineage[str(lin)][gen])
            _id = entry['_id']
            entries[_id] = {}
            entries[_id][target_function] = entry['properties'][target_function]
            if tag is not None:
                if tag in entry['properties']:
                    entries[_id][tag] = entry['properties'][tag]
                elif tag in entry['structure']:
                    entries[_id][tag] = entry['structure'][tag]
                elif tag == 'spacegroup':
                    st = searcher.population.get_structure(_id)
                    ss = CrystalSymmetry(st)
                    entries[_id][tag] = ss.number(1E-3)
                else:
                    raise ValueError('Tag not found')
            entries[_id]['x'] = np.cos(counter * 2 * np.pi / nele)
            entries[_id]['y'] = np.sin(counter * 2 * np.pi / nele)
            counter += 1

        ids = [i for i in entries]
        xx = np.array([entries[i]['x'] for i in ids])
        yy = np.array([entries[i]['y'] for i in ids])
        aa = np.array([4E3 for i in ids])
        cc = np.array([entries[i][target_function] for i in ids])
        best_candidate = entries[ids[cc.argsort()[0]]]
        plt.scatter(radius * xx, radius * yy, s=radius ** 2 * aa, c=cc, alpha=0.5 * radius, edgecolors='none')
        if tag is not None:
            for i in entries:
                if i == best_candidate:
                    ax.text(radius * entries[i]['x'], radius * entries[i]['y'], str(entries[i][tag]),
                            va='center',
                            ha='center', size=int(radius * 24), color='red')
                else:
                    ax.text(radius * entries[i]['x'], radius * entries[i]['y'], str(entries[i][tag]),
                            va='center',
                            ha='center', size=int(radius * 14))
    plt.savefig(searcher.population.name + '_circular.png')


gen_size = 80
letterpage = (11, 8.5)

dbsettings = {'host': 'mongo01.systems.wvu.edu', 'name': 'OrbitalGAF', 'user': 'gufranco', 'passwd': 'zxcvbnm'}
pcdb = pychemia.db.get_database(dbsettings)
popu = pychemia.population.orbitaldftu.OrbitalDFTU(pcdb)
searcher = pychemia.searcher.FireFly(popu, generation_size=80)
searcher.recover()

if __name__ == '__main__':
    pdf = PdfPages(name + 'multipage.pdf')

    print(popu)

    colors = ['r', 'g', 'm', 'c', 'y', 'k']
    color_tags = {}
    for i in tags:
        if len(colors) == 0:
            colors = ['r', 'g', 'm', 'c', 'y', 'k']
        color_tags[i] = colors.pop(0)

    plot_generation_chart(searcher, gen_size)

    sys.exit(1)

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
    popu = pychemia.population.RelaxStructures(pcdb, distance_tol=0.4, value_tol=0.1)
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
