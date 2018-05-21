from ..codes import CodeOutput
import numpy as np


class OctopusOutput(CodeOutput):

    @property
    def is_finished(self):
        return False

    def __init__(self, dirname):
        self.dirname = dirname
        CodeOutput.__init__(self)

    def read(self):
        coords = coordinates(self.dirname)
        energs = energy(self.dirname)
        self.output_values = {'energetics': energs, 'coordinates': coords}


def coordinates(dirname):
    """
    Reads the 'coordinates' file and extracts its
    data as the numpy arrays 'positions','velocities'
    and 'forces'.

    Each of those arrays has the following structure:

    Args:
       dirname:
          The path were the octopus data is produced

    Returns:
       natom:
          Number of atoms
       iterations:
          Array with iterations
       times:
          Array with times
       positions:
          positions[0:MaxIteration,0:natom,0:3]
       velocities:
          velocities[0:MaxIteration,0:natom,0:3]
       forces:
          forces[0:MaxIteration,0:natom,0:3]

    """

    geo = np.loadtxt(dirname + '/td.general/coordinates')
    natom = (geo.shape[1] - 2) / 9
    iterations = geo[:, 0]
    times = geo[:, 1]

    positions = geo[:, 2:3 * natom + 2]
    velocities = geo[:, 3 * natom + 2:6 * natom + 2]
    forces = geo[:, 6 * natom + 2:9 * natom + 2]

    positions = positions.reshape([-1, 3])
    velocities = velocities.reshape([-1, 3])
    forces = forces.reshape([-1, 3])

    return natom, iterations, times, positions, velocities, forces


def energy(dirname):
    """
    Reads the 'energy' file and extracts its
    data as the numpy arrays 'Total','Kinetic', etc.

    Args:
       dirname:
          The directory where the octopus data was produced

    Returns:
       iterations:
          numpy array of iterations recorded
       times:
          numpy array with time for each iteration
       energy_dict:
          dictionary of numpy arrays with values of
          Total, Kinetic (ions), Ion-Ion, Electronic,
          Eigenvalues, Hartree, Int[n v_xc], Exchange
          and Correlation

    """

    ene = np.loadtxt(dirname + '/td.general/energy')
    iterations = ene[:, 0]
    times = ene[:, 1]

    energy_dict = {'Total': ene[:, 2], 'Kinetic': ene[:, 3], 'Ion-Ion': ene[:, 4], 'Electronic': ene[:, 5],
                   'Eigenvalues': ene[:, 6], 'Hartree': ene[:, 7], 'Int': ene[:, 8], 'Exchange': ene[:, 9],
                   'Correlation': ene[:, 10]}

    return iterations, times, energy_dict
