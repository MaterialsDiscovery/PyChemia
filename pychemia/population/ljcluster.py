
import math
import sys
import numpy as np
import scipy.spatial
from pychemia import Composition, Structure, pcm_log
from pychemia.analysis import ClusterAnalysis, ClusterMatch
from pychemia.code.lennardjones import lj_compact_evaluate
from pychemia.utils.mathematics import unit_vector, length_vectors, unit_vectors, rotate_towards_axis, length_vector
from pychemia.utils.periodic import covalent_radius, atomic_number
from pychemia.utils.serializer import generic_serializer
from pychemia.code.lennardjones import LennardJones
from pychemia.external.symmol import get_point_group
from ._population import Population
from ._distances import FingerPrints, StructureDistances


class LJCluster(Population):

    def __init__(self, name, composition=None, tag='global', target_forces=1E-3, value_tol=1E-2,
                 distance_tolerance=0.1, minimal_density=70.0, refine=True, direct_evaluation=False):
        if composition is not None:
            self.composition = Composition(composition)
        else:
            self.composition = None
        Population.__init__(self, name=name, tag=tag, direct_evaluation=direct_evaluation,
                            distance_tolerance=distance_tolerance)
        self.tag = tag
        self.target_forces = target_forces
        self.value_tol = value_tol
        self.minimal_density = minimal_density
        self.refine = refine
        self.fingerprinter = FingerPrints(self.pcdb)
        self.distancer = StructureDistances(self.pcdb)

    def add_random(self):
        """
        Add one random structure to the population
        """
        if self.composition is None:
            raise ValueError('No composition associated to this population')
        comp = self.composition.composition.copy()
        structure = Structure.random_cluster(composition=comp)

        return self.new_entry(structure), None

    def get_duplicates(self, ids, tolerance=None, fast=True):
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
                ident2 = selection[i + 1]
                pcm_log.debug('Testing distances between %s and %s' % (str(ident1), str(ident2)))
                distance = self.distance(ident1, ident2)
                if distance < self.distance_tolerance:
                    pcm_log.debug('Distance %7.3f < %7.3f' % (distance, self.distance_tolerance))
                    ret[ident2] = ident1
        if len(ret) > 0:
            pcm_log.debug('Number of duplicates %d' % len(ret))
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

        new_structure = Structure(positions=new_pos0, symbols=entry0['structure']['symbols'],
                                  periodicity=False)
        entry_id = self.new_entry(structure=new_structure)
        new_structure = Structure(positions=new_pos1, symbols=entry0['structure']['symbols'],
                                  periodicity=False)
        entry_jd = self.new_entry(structure=new_structure)

        return entry_id, entry_jd

    def distance(self, entry_id, entry_jd, rcut=50):
        """
        Return a measure of the distance between two clusters by computing
        a n-dimensional vector of the distances between each atom to the
        origin and

        :param rcut:
        :param entry_id: The id of one population entry
        :param entry_jd: The id of another population entry
        :return: (int) The distance between two clusters
        """

        ids_pair = tuple(np.sort([entry_id, entry_jd]))
        distance_entry = self.distancer.get_distance(ids_pair)

        if distance_entry is None:
            fingerprints = {}
            for entry_ijd in [entry_id, entry_jd]:

                if self.fingerprinter.get_fingerprint(entry_ijd) is None:
                    structure = self.get_structure(entry_ijd)
                    analysis = ClusterAnalysis(structure)
                    x, ys = analysis.discrete_radial_distribution_function()
                    fingerprint = {'_id': entry_ijd}
                    for k in ys:
                        atomic_number1 = atomic_number(k[0])
                        atomic_number2 = atomic_number(k[1])
                        pair = '%06d' % min(atomic_number1 * 1000 + atomic_number2,
                                            atomic_number2 * 1000 + atomic_number1)
                        fingerprint[pair] = list(ys[k])

                    if self.fingerprinter.get_fingerprint(entry_ijd) is None:
                        self.fingerprinter.set_fingerprint(fingerprint)
                    else:
                        self.fingerprinter.update(entry_ijd, fingerprint)
                    fingerprints[entry_ijd] = fingerprint
                else:
                    fingerprints[entry_ijd] = self.fingerprinter.get_fingerprint(entry_ijd)

            dij = []
            for pair in fingerprints[entry_id]:
                if pair in fingerprints[entry_jd] and pair != '_id':
                    vect1 = fingerprints[entry_id][pair]
                    vect2 = fingerprints[entry_jd][pair]
                    if len(vect1) < len(vect2):
                        tmp = np.zeros(len(vect2))
                        tmp[:len(vect1)] = vect1
                        vect1 = tmp
                    elif len(vect1) > len(vect2):
                        tmp = np.zeros(len(vect1))
                        tmp[:len(vect2)] = vect2
                        vect2 = tmp
                    uvect1 = unit_vector(vect1)
                    uvect2 = unit_vector(vect2)
                    dij.append(0.5 * (1.0 - np.dot(uvect1, uvect2)))
            distance = float(np.mean(dij))
            self.distancer.set_distance(ids_pair, distance)
        else:
            distance = distance_entry['distance']
        return distance

    @property
    def to_dict(self):
        return {'name': self.name,
                'tag': self.tag,
                'target_forces': self.target_forces,
                'value_tol': self.value_tol,
                'distance_tolerance': self.distance_tolerance,
                'minimal_density': self.minimal_density}

    def get2_duplicates(self, ids, fast=True):
        dupes_dict = {}
        dupes_list = []
        selection = self.ids_sorted(ids)
        pcm_log.debug('Searching duplicates in %d structures' % len(selection))
        for i in range(len(selection) - 1):
            ncomps = 0
            entry_id = selection[i]
            if fast and entry_id in dupes_list:
                continue
            sys.stdout.write(" %5d of %5d: " % (i, len(selection)))
            value_i = self.value(entry_id)
            for j in range(i + 1, len(selection)):
                entry_jd = selection[j]
                if fast and entry_jd in dupes_list:
                    continue
                value_j = self.value(entry_jd)
                if abs(value_i - value_j) < self.value_tol:
                    ncomps += 1
                    distance = self.distance(entry_id, entry_jd)
                    if distance < self.distance_tolerance:
                        if entry_id in dupes_dict:
                            dupes_dict[entry_id].append(entry_jd)
                        else:
                            if fast:
                                dupes_dict[entry_id] = [entry_jd]
                            else:
                                dupes_dict[entry_id] = entry_jd
                        dupes_list.append(entry_jd)
            sys.stdout.write(' comparisons: %d\n' % ncomps)
        return dupes_dict, [x for x in selection if x in dupes_list]

    def is_evaluated(self, entry_id):

        entry = self.get_entry(entry_id)
        if entry is not None and 'properties' not in entry:
            raise ValueError('Anomalous condition for %s' % entry_id)

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
                         distance_tolerance=population_dict['distance_tolerance'],
                         minimal_density=population_dict['minimal_density'])

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):

        st_orig = self.get_structure(entry_id)
        st_dest = self.get_structure(entry_jd)

        cm = ClusterMatch(st_orig, st_dest)
        cm.match()

        # pos_orig = np.array(entry_orig['structure']['positions']).reshape((-1, 3))
        # pos_dest = np.array(entry_dest['structure']['positions']).reshape((-1, 3))
        pos_orig = cm.structure1.positions
        pos_dest = cm.structure2.positions

        # Move to a position with negative energy
        reduc = 1
        new_positions = np.array(pos_orig)
        while True:
            new_positions = rotation_move(pos_orig, pos_dest, fraction=reduc * factor)
            new_structure = Structure(positions=new_positions, symbols=st_orig.symbols, periodicity=False)
            lj = LennardJones(new_structure)
            if lj.get_energy() < 0.0:
                break
            reduc -= 0.05
            pcm_log.debug('Effective factor will be reduced to %7.3f, original factor %7.3f' % (reduc * factor, factor))
            if reduc <= 0.0:
                # print 'No movement effective'
                break

        # Avoid condition with atoms too close
        distance_matrix = scipy.spatial.distance_matrix(new_positions, new_positions)
        tmp = np.max(distance_matrix.flatten())
        # print 'Scaling by', tmp
        minimal_distance = np.min((distance_matrix + tmp * np.eye(len(new_positions))).flatten())

        if minimal_distance < 1E-8:
            pcm_log.debug("Null distance between different atoms, no moving")
            new_positions = pos_orig

        if tmp > 5:
            # print 'Big scaling, better not to move'
            new_positions = pos_orig
        else:
            max_cov = np.max(covalent_radius(st_orig.symbols))
            new_positions *= max_cov / minimal_distance

        new_structure = Structure(positions=new_positions, symbols=st_orig.symbols, periodicity=False)
        # print 'Density of cluster', new_structure.density

        if in_place:
            self.unset_properties(entry_id)
            return self.set_structure(entry_id, new_structure)
        else:
            return self.new_entry(new_structure, active=False)

    def evaluate(self, structure, gtol=None):

        if gtol is None:
            gtol = self.target_forces

        positions, forces, energy = lj_compact_evaluate(structure, gtol, self.minimal_density)

        structure.set_positions(positions)
        structure.relocate_to_cm()
        if structure.natom > 2:
            structure.align_inertia_momenta()
        sorted_indices = structure.sort_sites()
        forces = forces[sorted_indices]
        pg = get_point_group(structure, executable='symmol')
        properties = {'forces': generic_serializer(forces), 'energy': energy, 'point_group': pg}
        return structure, properties

    def evaluate_entry(self, entry_id):
        pcm_log.debug('Evaluating %s target density= %7.3F' % (entry_id, self.minimal_density))
        structure = self.get_structure(entry_id)
        structure, properties = self.evaluate(structure, gtol=self.target_forces)
        self.update_properties(entry_id=entry_id, new_properties=properties)
        return self.set_structure(entry_id, structure)

    def refine(self, entry_id, gtol=None):
        if self.refine:
            pcm_log.debug('Evaluating %s target density= %7.3F' % (entry_id, self.minimal_density))
            structure = self.get_structure(entry_id)
            structure, properties = self.evaluate(structure, gtol=gtol)
            self.update_properties(entry_id=entry_id, new_properties=properties)
            return self.set_structure(entry_id, structure)

    def maxforce(self, entry_id):
        return np.max(length_vectors(self.get_forces(entry_id)))

    def refine_progressive(self, entry_id):

        if self.refine:
            inivalue = self.value(entry_id)
            gtol = 10 ** math.ceil(math.log10(self.maxforce(entry_id)))
            while True:
                pcm_log.debug('Local minimization up to %7.3f ' % gtol)
                gtol /= 10
                pcm_log.debug('Evaluating %s target density= %7.3F' % (entry_id, self.minimal_density))
                structure = self.get_structure(entry_id)
                structure, properties = self.evaluate(structure, gtol=gtol)
                if properties['energy'] / structure.natom < inivalue:
                    self.update_properties(entry_id=entry_id, new_properties=properties)
                    return self.set_structure(entry_id, structure)
                else:
                    pcm_log.debug('Relaxation raise value %7.3f < %7.3f' % (inivalue, properties['energy'] /
                                                                            structure.natom))
                if self.maxforce(entry_id) > gtol:
                    pcm_log.debug('I cannot relax more...')
                    break

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):

        entry = self.get_entry(entry_id)
        pos = np.array(entry['structure']['positions']).reshape((-1, 3))
        # Unit Vectors
        uv = unit_vectors(2 * np.random.rand(*pos.shape) - 1)
        new_pos = generic_serializer(pos + factor * uv)

        structure = Structure(positions=new_pos, symbols=entry['structure']['symbols'], periodicity=False)

        if in_place:
            self.update_properties(entry_id=entry_id, new_properties={})
            return self.set_structure(entry_id, structure)
        else:
            structure = Structure(positions=new_pos, symbols=entry['structure']['symbols'], periodicity=False)
            return self.new_entry(structure, active=False)

    def get_structure(self, entry_id):
        entry = self.get_entry(entry_id)
        if 'structure' not in entry:
            raise ValueError('structure is not present on %s' % entry_id)
        if entry['structure'] is None:
            raise ValueError('structure is None for %s' % entry_id)
        return Structure.from_dict(entry['structure'])

    def get_forces(self, entry_id):
        entry = self.get_entry(entry_id, projection={'properties.forces': 1})
        forces = np.array(entry['properties']['forces']).reshape((-1, 3))
        return forces

    def str_entry(self, entry_id):
        structure = self.get_structure(entry_id)
        entry = self.get_entry(entry_id, projection={'properties': 1})
        msg = 'Cluster: LJ%d  Point group: %s  Energy: %7.3f Forces: %7.1E'
        pg = entry['properties']['point_group']
        return msg % (structure.natom, pg, entry['properties']['energy'], self.maxforce(entry_id))

    def new_entry(self, structure, active=True):

        if active and self.direct_evaluation:
            structure, properties = self.evaluate(structure, gtol=self.target_forces)
        else:
            properties = {}
        status = {self.tag: active}
        entry = {'structure': structure.to_dict, 'properties': properties, 'status': status}
        entry_id = self.set_entry(entry)
        pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def recover(self):
        data = self.get_population_info()
        if data is not None:
            self.distance_tolerance = data['distance_tolerance']
            self.value_tol = data['value_tol']
            self.name = data['name']
            self.target_forces = data['target_forces']
            self.minimal_density = data['minimal_density']

    def value(self, entry_id):
        entry = self.get_entry(entry_id)
        structure = self.get_structure(entry_id)
        if 'properties' not in entry:
            pcm_log.debug('This entry has no properties %s' % str(entry['_id']))
            return None
        elif entry['properties'] is None:
            return None
        elif 'energy' not in entry['properties']:
            pcm_log.debug('This entry has no energy in properties %s' % str(entry['_id']))
            return None
        else:
            return entry['properties']['energy'] / structure.get_composition().gcd


def rotation_move(pos_orig, pos_dest, fraction):
    new_positions = np.zeros(pos_orig.shape)
    for i in range(len(pos_orig)):
        if length_vector(pos_orig[i]) < 1E-5 or length_vector(pos_dest[i]) < 1E-5:
            new_positions[i] = pos_dest[i]
        else:
            new_positions[i] = rotate_towards_axis(pos_orig[i], pos_dest[i], fraction=fraction)
            uv = new_positions[i] / np.linalg.norm(new_positions[i])
            new_positions[i] = fraction * np.linalg.norm(pos_dest[i]) * uv + (1 - fraction) * \
                                np.linalg.norm(pos_orig[i]) * uv
    return new_positions


def direct_move(pos_orig, pos_dest, fraction):
    return pos_orig + fraction * (pos_dest - pos_orig)


def movement_sweep(pos_orig, pos_dest, symbols, figname='figure.pdf'):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(ncols=1, nrows=3, sharex=True, figsize=(11, 8.5))
    plt.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.08, hspace=0.08)

    ee = []
    ff = []
    dd = []
    delta = 2E-3
    xx = np.arange(0.0, 1.0 + 0.9 * delta, delta)

    for f in xx:
        new_positions = direct_move(pos_orig, pos_dest, fraction=f)

        new_structure = Structure(positions=new_positions, symbols=symbols, periodicity=False)
        lj = LennardJones(new_structure)
        ee.append(lj.get_energy())
        ff.append(np.max(lj.get_magnitude_forces()))
        # Distance Matrix
        dm = scipy.spatial.distance_matrix(new_positions, new_positions)
        # Min distance
        md = np.min(np.array(np.array(dm) + 100 * np.eye(len(pos_orig))).flatten())
        dd.append(md)

    ax[0].plot(xx, ee)
    ax[0].set_ylim(min(ee), 0.1)
    ax[1].semilogy(xx, ff)
    ax[2].plot(xx, dd)

    st = Structure(positions=pos_orig, symbols=symbols, periodicity=False)
    lj = LennardJones(st)
    ax[0].plot(0, lj.get_energy(), 'ro')
    ax[1].semilogy(0, np.max(lj.get_magnitude_forces()), 'ro')
    dm = scipy.spatial.distance_matrix(lj.structure.positions, lj.structure.positions)
    md = np.min(np.array(np.array(dm) + 100 * np.eye(len(pos_orig))).flatten())
    ax[2].plot(0, md, 'ro')

    st = Structure(positions=pos_dest, symbols=symbols, periodicity=False)
    lj = LennardJones(st)
    ax[0].plot(1, lj.get_energy(), 'ro')

    ax[1].semilogy(1, np.max(lj.get_magnitude_forces()), 'ro')
    dm = scipy.spatial.distance_matrix(lj.structure.positions, lj.structure.positions)
    md = np.min(np.array(np.array(dm) + 100 * np.eye(len(pos_orig)).flatten()))
    ax[2].plot(1, md, 'ro')

    ax[2].set_xlim(-0.01, 1.01)

    ax[0].set_ylabel('Energy')
    ax[1].set_ylabel('Max Force')
    ax[2].set_ylabel('Minimal inter atomic distance')

    plt.savefig(figname)
