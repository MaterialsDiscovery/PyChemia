import os
import time
import numpy as np
from threading import Thread
from ase.calculators.lj import LennardJones
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms, UnitCellFilter
import pychemia
import pychemia.external.ase
from pychemia.utils.serializer import generic_serializer


class AseObjectiveFunction:
    def __init__(self):
        self.population = None

    def initialize(self, population):
        self.population = population

    def ids_sorted(self, selection):
        values = np.array([self.population.value(i) for i in selection])
        argsort = np.argsort(values)
        return np.array(selection)[argsort]

    def get_values(self, selection):
        ret = {}
        for i in selection:
            ret[i] = self.population.value(i)
        return ret


class AseEvaluator:
    def __init__(self):
        self.process = None
        self.thread = None
        self.population = None

    def initialize(self, population):
        self.population = population

    def evaluate(self, imember):
        entry = self.population.get_entry(imember)
        pcm_structure = pychemia.Structure.from_dict(entry['structure'])

        ase_structure = pychemia.external.ase.pychemia2ase(pcm_structure)
        ase_structure.set_calculator(LennardJones())

        dyn = QuasiNewton(ase_structure)
        dyn.run()

        ase_structure.set_constraint(FixAtoms(mask=[True for atom in ase_structure]))
        ucf = UnitCellFilter(ase_structure)
        qn = QuasiNewton(ucf)
        qn.run()

        new_structure = pychemia.external.ase.ase2pychemia(ase_structure)
        energy = ase_structure.get_potential_energy()
        forces = ase_structure.get_forces()
        stress = ase_structure.get_stress()
        new_properties = {'energy': float(energy), 'forces': generic_serializer(forces),
                          'stress': generic_serializer(stress)}

        self.population.db.update(imember, structure=new_structure, properties=new_properties)

    @property
    def is_running(self):
        if self.thread is not None:
            return self.thread.is_alive()
        else:
            return False

    def run(self):

        def worker(evaluator, population):
            while True:
                for i in population.actives:
                    if not population.is_evaluated(i):
                        evaluator.evaluate(i)
                time.sleep(1)
                if os.path.exists('stop'):
                    os.remove('stop')
                    return

        # self.process = Process(target=worker, args=(self.population,))
        # self.process.start()
        self.thread = Thread(target=worker, args=(self, self.population,))
        self.thread.daemon = True
        self.thread.start()

    def stop(self):
        if self.process is not None and not self.process.is_alive():
            self.process.terminate()
