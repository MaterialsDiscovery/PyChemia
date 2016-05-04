import numpy as np
import pychemia
import scipy.optimize
import scipy.spatial


class LennardJones:
    def __init__(self, structure, ljparams=None, cp=0.0):
        self.initial_structure = structure
        self.structure = structure.copy()

        if ljparams is None:
            ljparams = {}
            for i in self.structure.species:
                for j in self.structure.species:
                    ljparams[i + '-' + j] = {'sigma': 1.0, 'epsilon': 1.0}
        self.ljparams = ljparams
        self.cp = cp
        self.sigmas = np.zeros((structure.natom, structure.natom))
        self.epsilons = np.zeros((structure.natom, structure.natom))
        for i in range(structure.natom):
            for j in range(structure.natom):
                self.sigmas[i, j], self.epsilons[i, j] = self.get_sigma_epsilon(i, j)

    def get_forces(self):
        return lj_forces(self.structure.positions, self.sigmas, self.epsilons, self.cp).reshape((-1, 3))

    def get_magnitude_forces(self):
        return pychemia.utils.mathematics.length_vectors(self.get_forces())

    def steepest_descent(self, dt=1, tolerance=1E-3):
        forces = self.get_forces()
        while np.max(self.get_magnitude_forces()) > tolerance:
            for i in range(self.structure.natom):
                mass = pychemia.utils.periodic.mass(self.structure.symbols[i])
                self.structure.positions[i] += 0.5 * forces[i] / mass * dt ** 2

            forces = self.get_forces()
            maxforce = np.max(self.get_magnitude_forces())
            dt = max(100.0 / maxforce, 1.0)
            print('Steepest Descent maxforce= %9.3E   density=%7.4f dt=%7.3f' % (maxforce, self.structure.density, dt))

    def get_sigma_epsilon(self, i, j):
        specie1 = self.structure.symbols[i]
        specie2 = self.structure.symbols[j]
        sort_species = [specie1, specie2]
        sort_species.sort()

        sigma = self.ljparams[sort_species[0] + '-' + sort_species[1]]['sigma']
        epsilon = self.ljparams[sort_species[0] + '-' + sort_species[1]]['epsilon']

        return float(sigma), float(epsilon)

    def get_energy(self):
        return lj_energy(self.structure.positions, self.sigmas, self.epsilons, self.cp)

    def local_minimization(self, method='BFGS', gtol=1E-4, soft_max_ncalls=5, hard_max_ncalls=20):

        if method in ['Nelder-Mead', 'Powell']:
            jac = None
            options = {'maxiter': 100, 'disp': False}
        else:
            jac = lj_gradient
            options = {'gtol': 0.1 * gtol, 'disp': False}

        x0 = self.structure.positions.flatten()

        hard_call = 1
        soft_call = 1
        res = 0.0

        while True:

            res = scipy.optimize.minimize(lj_energy, x0, args=(self.sigmas, self.epsilons, self.cp), method=method,
                                          jac=jac, options=options)

            if method in ['BFGS', 'CG']:
                maxforce = np.max(np.apply_along_axis(np.linalg.norm, 1, np.array(res.jac).reshape((-1, 3))))
                print('[%2d/%2d] MaxForce/gtol= %9.3E/%9.3E  Value= %12.5f' % (soft_call,
                                                                               soft_max_ncalls,
                                                                               maxforce, gtol,
                                                                               res.fun))

                x0 = res.x
                soft_call += 1
                hard_call += 1
                if maxforce < gtol:
                    break
                if soft_call > soft_max_ncalls:
                    print('Anomalous condition, maxforce= %9.3E gtol= %9.3E' % (maxforce, gtol))
                    pos = np.array(x0).reshape((-1, 3))
                    mindis = np.min(np.array(np.array(scipy.spatial.distance_matrix(pos, pos)) +
                                             100 * np.eye(len(pos))).flatten())
                    x0 = np.array(1.0 / mindis * pos).flatten()
                    soft_call = 1
                if hard_call > hard_max_ncalls:
                    print('Could not reach target forces, maxforce= %9.3E gtol= %9.3E' % (maxforce, gtol))
                    print('Hard limit reached, aborting local minimization')
                    break
            else:
                break
        return res

    def move_to_minima(self):
        relax = self.local_minimization()
        self.structure.positions = relax.x.reshape((-1, 3))
        self.structure.relocate_to_cm()


def lj_forces(pos, sigmas, epsilons, cp=0.0):
    pos = np.array(pos).reshape((-1, 3))
    n = len(pos)
    ret = np.zeros((n, 3))
    for i in range(n - 1):
        for j in range(i + 1, n):
            vector = pos[j] - pos[i]
            distance = np.linalg.norm(vector)
            sr6 = (sigmas[i, j] / distance) ** 6
            magnitude = 24 * epsilons[i, j] / distance * (2.0 * sr6 ** 2 - sr6)
            ret[i] -= magnitude * vector / distance
            ret[j] += magnitude * vector / distance
    if cp > 0.0:
        for i in range(n):
            ret[i] -= cp * pos[i]
    return ret.flatten()


def lj_gradient(pos, sigmas, epsilons, cp=0.0):
    return -lj_forces(pos, sigmas, epsilons, cp)


def lj_energy(pos, sigmas, epsilons, cp=0.0):
    pos = np.array(pos).reshape((-1, 3))
    ret = 0
    n = len(pos)
    for i in range(n - 1):
        for j in range(i + 1, n):
            vector = pos[j] - pos[i]
            distance = np.linalg.norm(vector)
            if distance > 1E-10:
                sr6 = (sigmas[i, j] / distance) ** 6
            else:
                sr6 = 0
            ret += 4 * epsilons[i, j] * (sr6 ** 2 - sr6)
    if cp > 0.0:
        for i in range(n):
            ret += 0.5 * cp * np.linalg.norm(pos[i]) ** 2
    return ret


def lj_compact_evaluate(structure, gtol, minimal_density):
    k = 1
    while structure.density < minimal_density:
        iniden = structure.density
        lj = pychemia.code.LennardJones(structure, cp=k)
        relax = lj.local_minimization(gtol=gtol)
        structure.set_positions(relax.x.reshape((-1, 3)))
        finden = structure.density
        print('Compacting Cluster (target_density=%7.3f): Density I= %7.3f   F= %7.3f' % (minimal_density,
                                                                                          iniden, finden))
        k += 1
        if k > 10:
            print('I tried too much...')
            break

    lj = pychemia.code.LennardJones(structure)
    relax = lj.local_minimization(gtol=gtol)
    # Return relaxed positions, forces and energy
    return relax.x.reshape((-1, 3)), relax.jac.reshape((-1, 3)), relax.fun
