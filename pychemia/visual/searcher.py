import sys
import numpy as np
from pychemia.crystal import CrystalSymmetry
import matplotlib.pyplot as plt
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
        x = np.array([entries[i]['x'] for i in ids])
        y = np.array([entries[i]['y'] for i in ids])
        a = np.array(len(ids)*[4E3])
        c = np.array([entries[i][target_function] for i in ids])
        best_candidate = entries[ids[c.argsort()[0]]]
        plt.scatter(radius * x, radius * y, s=radius ** 2 * a, c=c, alpha=0.5 * radius, edgecolors='none')
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
