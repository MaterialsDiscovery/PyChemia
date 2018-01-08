from abc import ABCMeta, abstractmethod
import numpy as np


class Relaxator:
    __metaclass__ = ABCMeta

    def __init__(self, target_forces):
        self.target_forces = target_forces

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def get_forces_stress_energy(self):
        return None, None, None

    @abstractmethod
    def get_final_geometry(self):
        pass

    def get_max_force_stress(self):

        forces, stress, total_energy = self.get_forces_stress_energy()

        if forces is not None:
            max_force = np.max(np.apply_along_axis(np.linalg.norm, 1, forces))
        else:
            max_force = None

        if stress is not None:
            max_stress = np.max(np.abs(stress.flatten()))
        else:
            max_stress = None

        return max_force, max_stress

    def relaxation_status(self):

        max_force, max_stress = self.get_max_force_stress()

        if max_force is None:
            good_forces = False
        elif max_force < self.target_forces:
            good_forces = True
        else:
            good_forces = False

        if max_stress is None:
            good_stress = False
        elif max_stress < self.target_forces:
            good_stress = True
        else:
            good_stress = False

        return good_forces, good_stress
