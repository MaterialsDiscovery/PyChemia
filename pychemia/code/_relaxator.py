from abc import ABCMeta, abstractmethod
import numpy as np
from pychemia import pcm_log


class Relaxator:
    __metaclass__ = ABCMeta

    def __init__(self, target_forces):
        self.target_forces = target_forces

    @abstractmethod
    def set_params(self, params):
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def get_forces_stress_energy(self):
        pass

    @abstractmethod
    def get_final_geometry(self):
        pass

    def get_max_force_stress(self):

        forces, stress, total_energy = self.get_forces_stress_energy()

        if forces is not None:
            max_force = np.max(np.abs(forces.flatten()))
        else:
            max_force = None

        if stress is not None:
            max_stress = np.max(np.abs(stress.flatten()))
        else:
            max_stress = None

        return max_force, max_stress

    def relaxation_status(self):

        forces, stress, total_energy = self.get_forces_stress_energy()

        if forces is None or stress is None:
            return False, False, False

        # Forces are good if we are at least one order of magnitude higher than target
        if forces is not None and np.max(np.abs(forces.flatten())) < 10 * self.target_forces:
            good_forces = True
        else:
            if forces is not None:
                pcm_log.info('Forces: ' + str(np.max(np.abs(forces.flatten()), 10 * self.target_forces)))
            good_forces = False

        # Stress is good if we are at least one order of magnitude higher than target
        if stress is not None and np.max(np.abs(stress.flatten())) < 10 * self.target_forces:
            good_stress = True
        else:
            if stress is not None:
                pcm_log.info('Stress: ' + str(np.max(np.abs(stress.flatten()), 10 * self.target_forces)))
            good_stress = False
        good_data = True

        return good_forces, good_stress, good_data
