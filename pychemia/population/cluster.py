import numpy as np
from ._population import Population
import pychemia
import math
import scipy.spatial


class LJCluster(Population):

    def __init__(self, name, composition=None, tag='global', target_forces=1E-3, value_tol=1E-2,
                 distance_tol=0.25):
        if composition is not None:
            self.composition = pychemia.Composition(composition)
        else:
            self.composition = None
        self.tag = tag
        self.target_forces = target_forces
        self.value_tol = value_tol
        self.distance_tol = distance_tol
        Population.__init__(self, name, tag)

    def add_random(self):
        """
        Add one random structure to the population
        """
        if self.composition is None:
            raise ValueError('No composition associated to this population')
        comp = self.composition.composition.copy()
        structure = pychemia.Structure.random_cluster(composition=comp)

        return self.new_entry(structure), None

    def check_duplicates(self, ids):
        ret = {}
        selection = self.ids_sorted(ids)
        values = np.array([self.value(i) for i in selection])
        if len(values) == 0:
            return ret
        diffs = np.ediff1d(values)

        for i in range(len(diffs)):
            idiff = diffs[i]
            if idiff < self.value_tol:
                ident1 = selection[i]
                ident2 = selection[i+1]
                pychemia.pcm_log.debug('Testing distances between %s and %s' % (str(ident1), str(ident2)))
                distance = self.distance(ident1, ident2)
                # print 'Distance = ', distance
                if distance < self.distance_tol:
                    pychemia.pcm_log.debug('Distance %7.3f < %7.3f' % (distance, self.distance_tol))
                    ret[ident2] = ident1
        if len(ret) > 0:
            pychemia.pcm_log.debug('Number of duplicates %d' % len(ret))
        return ret

    def cross(self, ids):
        if len(ids) != 2:
            raise ValueError("Crossing only implemented between two clusters")

        entry0 = self.get_entry(ids[0])
        entry1 = self.get_entry(ids[1])

        pos0 = np.array(entry0['structure']['positions']).reshape((-1, 3))
        pos1 = np.array(entry1['structure']['positions']).reshape((-1, 3))

        cut = np.random.randint(1, len(pos0))

        new_pos0 = np.concatenate((pos0[:cut], pos1[cut:]))
        new_pos1 = np.concatenate((pos1[:cut], pos0[cut:]))

        new_structure = pychemia.Structure(positions=new_pos0, symbols=entry0['structure']['symbols'],
                                           periodicity=False)
        entry_id = self.new_entry(structure=new_structure)
        new_structure = pychemia.Structure(positions=new_pos1, symbols=entry0['structure']['symbols'],
                                           periodicity=False)
        entry_jd = self.new_entry(structure=new_structure)

        return entry_id, entry_jd

    def distance(self, entry_id, entry_jd):

        entry_orig = self.get_entry(entry_id)
        entry_dest = self.get_entry(entry_jd)

        pos_orig = np.array(entry_orig['structure']['positions']).reshape((-1, 3))
        pos_dest = np.array(entry_dest['structure']['positions']).reshape((-1, 3))

        dist_cm_orig = np.apply_along_axis(np.linalg.norm, 1, pos_orig)
        dist_cm_dest = np.apply_along_axis(np.linalg.norm, 1, pos_dest)

        return np.linalg.norm(np.abs(dist_cm_dest-dist_cm_orig))

    @property
    def to_dict(self):
        return {'name': self.name,
                'tag': self.tag,
                'target_forces': self.target_forces,
                'value_tol': self.value_tol,
                'distance_tol': self.distance_tol}

    def get_duplicates(self, ids, fast=False):
        dupes_dict = {}
        dupes_list = []
        selection = self.ids_sorted(ids)
        print 'Searching duplicates in %d structures' % len(selection)
        for i in range(len(selection)-1):
            print i, 'of', len(selection)
            entry_id = selection[i]
            value_i = self.value(entry_id)
            for j in range(i+1, len(selection)):
                entry_jd = selection[j]
                if fast and entry_jd in dupes_list:
                    continue
                value_j = self.value(entry_jd)
                if abs(value_i-value_j) < self.value_tol:
                    distance = self.distance(entry_id, entry_jd)
                    if distance < self.distance_tol:
                        if entry_id in dupes_dict:
                            dupes_dict[entry_id].append(entry_jd)
                        else:
                            dupes_dict[entry_id] = [entry_jd]
                        dupes_list.append(entry_jd)
        return dupes_dict, [x for x in selection if x in dupes_list]

    def is_evaluated(self, entry_id):
        entry = self.get_entry(entry_id)
        if entry is not None and entry['properties'] is not None:
            properties = entry['properties']
            if 'forces' not in properties:
                forces = None
            elif properties['forces'] is None:
                forces = None
            else:
                forces = np.max(np.apply_along_axis(np.linalg.norm, 1, np.array(properties['forces']).reshape((-1, 3))))
        else:
            forces = None

        if forces is not None and forces < self.target_forces:
            return True
        else:
            return False

    def from_dict(self, population_dict):
        return LJCluster(name=population_dict['name'],
                         tag=population_dict['tag'],
                         target_forces=population_dict['target_forces'],
                         value_tol=population_dict['value_tol'],
                         distance_tol=population_dict['distance_tol'])

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):

        entry_orig = self.get_entry(entry_id)
        entry_dest = self.get_entry(entry_jd)

        pos_orig = np.array(entry_orig['structure']['positions']).reshape((-1, 3))
        pos_dest = np.array(entry_dest['structure']['positions']).reshape((-1, 3))

        newpos = parametric_move(pos_orig, pos_dest, fraction=f)

        new_structure = pychemia.Structure(positions=newpos,
                                           symbols=symbols, # entry_orig['structure']['symbols'],
                                           periodicity=False)

        if in_place:
            structure, properties, relax = self.evaluate(new_structure)
            return self.pcdb.update(entry_id, structure=structure, properties=properties)
        else:
            return self.new_entry(new_structure, active=False)

    def evaluate(self, structure, gtol=None):

        if gtol is None:
            gtol = self.target_forces
        lj = pychemia.code.LennardJones(structure)
        relax = lj.local_minimization(gtol=gtol)
        structure.set_positions(relax.x.reshape((-1, 3)))
        structure.relocate_to_cm()
        structure.align_inertia_momenta()
        sorted_indices = structure.sort_sites()
        forces = relax.jac.reshape((-1, 3))[sorted_indices]
        properties = {'forces': pychemia.serializer.generic_serializer(forces), 'energy': relax.fun}
        return structure, properties, relax

    def refine(self, entry_id, gtol=None):
        structure = self.get_structure(entry_id)
        structure, properties, relax = self.evaluate(structure, gtol=gtol)
        return self.pcdb.update(entry_id, structure=structure, properties=properties)

    def maxforce(self, entry_id):
        return np.max(pychemia.utils.mathematics.length_vectors(self.get_forces(entry_id)))

    def refine_progressive(self, entry_id):
        gtol = 10**math.ceil(math.log10(self.maxforce(entry_id)))
        while True:
            print 'Local minimization up to ', gtol
            gtol /= 10
            structure = self.get_structure(entry_id)
            structure, properties, relax = self.evaluate(structure, gtol=gtol)
            self.pcdb.update(entry_id, structure=structure, properties=properties)
            if self.maxforce(entry_id) > gtol:
                print 'The relaxation was not successful'
                print relax
                break

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):

        entry = self.get_entry(entry_id)
        pos = np.array(entry['structure']['positions']).reshape((-1, 3))
        # Unit Vectors
        uv = pychemia.utils.mathematics.unit_vectors(2*np.random.rand(*pos.shape)-1)
        new_pos = pychemia.serializer.generic_serializer(pos+factor*uv)

        structure = pychemia.Structure(positions=new_pos,
                                               symbols=entry['structure']['symbols'],
                                               periodicity=False)

        if in_place:
            structure, properties, relax = self.evaluate(structure)
            return self.pcdb.db.pychemia_entries.update_one({'_id': entry_id}, {'$set':
                                                        {'structure': structure.to_dict, 'properties': properties}})
        else:
            structure = pychemia.Structure(positions=new_pos,
                                               symbols=entry['structure']['symbols'],
                                               periodicity=False)
            return self.new_entry(structure, active=False)

    def get_structure(self, entry_id):
        entry = self.get_entry(entry_id)
        if 'structure' not in entry:
            raise ValueError('structure is not present on %s' % entry_id)
        if entry['structure'] is None:
            raise ValueError('structure is None for %s' % entry_id)
        return pychemia.Structure.from_dict(entry['structure'])

    def get_forces(self, entry_id):
        entry = self.pcdb.db.pychemia_entries.find_one({'_id': entry_id}, {'properties.forces': 1})
        forces = np.array(entry['properties']['forces']).reshape((-1, 3))
        return forces

    def str_entry(self, entry_id):
        structure = self.get_structure(entry_id)
        return str(structure)

    def new_entry(self, structure, active=True):

        structure, properties, relax = self.evaluate(structure)
        status = {self.tag: active}
        entry_id = self.pcdb.insert(structure=structure, properties=properties, status=status)
        pychemia.pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def recover(self):
        data = self.pcdb.db.population_info.find_one({'tag': self.tag})
        if data is not None:
            self.distance_tol = data['distance_tol']
            self.value_tol = data['value_tol']
            self.name = data['name']
            self.target_forces = data['target_forces']

    def value(self, entry_id):
        entry = self.get_entry(entry_id)
        structure = self.get_structure(entry_id)
        if 'properties' not in entry:
            pychemia.pcm_log.debug('This entry has no properties %s' % str(entry['_id']))
            return None
        elif entry['properties'] is None:
            return None
        elif 'energy' not in entry['properties']:
            pychemia.pcm_log.debug('This entry has no energy in properties %s' % str(entry['_id']))
            return None
        else:
            return entry['properties']['energy'] / structure.get_composition().gcd


def rotation_move(pos_orig, pos_dest, fraction):

    newpos = np.zeros(pos_orig.shape)
    for i in range(len(pos_orig)):
        newpos[i] = pychemia.utils.mathematics.rotate_towards_axis(pos_orig[i], pos_dest[i], fraction=fraction)
        uv = newpos[i]/np.linalg.norm(newpos[i])
        newpos[i] = fraction * np.linalg.norm(pos_dest[i]) * uv + (1-fraction) * np.linalg.norm(pos_orig[i]) * uv
    return newpos


def direct_move(pos_orig, pos_dest, fraction):
    return pos_orig + fraction*(pos_dest - pos_orig)


def movement_sweep(pos_orig, pos_dest, symbols, figname='figure.pdf'):

    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(ncols=1, nrows=3, sharex=True, figsize=(11, 8.5))
    plt.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.08, hspace=0.08)

    ee = []
    ff = []
    dd = []
    delta = 2E-3
    xx = np.arange(0.0, 1.0+0.9*delta, delta)

    for f in xx:
        newpos = direct_move(pos_orig, pos_dest, fraction=f)

        new_structure = pychemia.Structure(positions=newpos,
                                           symbols=symbols, # entry_orig['structure']['symbols'],
                                           periodicity=False)
        lj = pychemia.code.LennardJones(new_structure)
        ee.append(lj.get_energy())
        ff.append(np.max(lj.get_magnitude_forces()))
        # Distance Matrix
        dm = scipy.spatial.distance_matrix(newpos, newpos)
        # Min distance
        md = np.min((np.array(dm)+100*np.eye(len(pos_orig))).flatten())
        dd.append(md)

    ax[0].plot(xx, ee)
    ax[0].set_ylim(min(ee), 0.1)
    ax[1].semilogy(xx, ff)
    ax[2].plot(xx, dd)

    st = pychemia.Structure(positions=pos_orig, symbols=symbols, periodicity=False)
    lj = pychemia.code.LennardJones(st)
    ax[0].plot(0, lj.get_energy(), 'ro')
    ax[1].semilogy(0, np.max(lj.get_magnitude_forces()), 'ro')
    dm = scipy.spatial.distance_matrix(lj.structure.positions, lj.structure.positions)
    md = np.min((np.array(dm)+100*np.eye(len(pos_orig))).flatten())
    ax[2].plot(0, md, 'ro')

    st = pychemia.Structure(positions=pos_dest, symbols=symbols, periodicity=False)
    lj = pychemia.code.LennardJones(st)
    ax[0].plot(1, lj.get_energy(), 'ro')

    ax[1].semilogy(1, np.max(lj.get_magnitude_forces()), 'ro')
    dm = scipy.spatial.distance_matrix(lj.structure.positions, lj.structure.positions)
    md = np.min((np.array(dm)+100*np.eye(len(pos_orig))).flatten())
    ax[2].plot(1, md, 'ro')

    ax[2].set_xlim(-0.01, 1.01)

    ax[0].set_ylabel('Energy')
    ax[1].set_ylabel('Max Force')
    ax[2].set_ylabel('Minimal interatomic distance')

    plt.savefig(figname)
