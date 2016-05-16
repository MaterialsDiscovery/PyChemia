import numpy as np
import matplotlib.pyplot as plt


def plot_evolution_circular(searcher, target_function='energy_pa', tags='spacegroup'):
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
            if tags is not None:
                entries[_id][tags] = entry['properties'][tags]
            entries[_id]['x'] = np.cos(counter * 2 * np.pi / nele)
            entries[_id]['y'] = np.sin(counter * 2 * np.pi / nele)
            counter += 1

        ids = [i for i in entries]
        X = np.array([entries[i]['x'] for i in ids])
        Y = np.array([entries[i]['y'] for i in ids])
        A = np.array([4E3 for i in ids])
        C = np.array([entries[i][target_function] for i in ids])
        best_candidate = entries[ids[C.argsort()[0]]]
        plt.scatter(radius * X, radius * Y, s=radius ** 2 * A, c=C, alpha=0.5 * radius, edgecolors='none')
        if tags is not None:
            for i in entries:
                if i == best_candidate:
                    ax.text(radius * entries[i]['x'], radius * entries[i]['y'], str(entries[i]['spacegroup']),
                            va='center',
                            ha='center', size=int(radius * 24), color='red')
                else:
                    ax.text(radius * entries[i]['x'], radius * entries[i]['y'], str(entries[i]['spacegroup']),
                            va='center',
                            ha='center', size=int(radius * 14))
    plt.savefig(searcher.population.name + '_circular.png')
