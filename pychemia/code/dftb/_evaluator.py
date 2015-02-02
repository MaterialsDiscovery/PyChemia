__author__ = 'Guillermo Avendano Franco'

import os
import bson
from pychemia.evaluator import Evaluator
from multiprocessing import Pool, Process
from pychemia.code.dftb import Relaxator, read_detailed_out, read_geometry_gen, read_dftb_stdout
from pychemia.serializer import generic_serializer
import logging

logging.basicConfig(level=logging.DEBUG)

def f(setup):
    logging.info(str(setup['id'])+' [Geometry optimization start] ')
    workdir = setup['workdir']
    population = setup['population']
    structure = population.get_structure(setup['id'])
    slater_path = setup['slater_path']
    target_forces = setup['target_forces']
    relax = Relaxator(workdir=workdir, structure=structure, slater_path=slater_path, target_forces=target_forces)
    relax.run()
    update(workdir, population, setup['id'])
    summary = str(setup['id'])+' [Geometry optimization end]\n'
    filename = workdir+os.sep+'dftb_stdout.log'
    if os.path.isfile(filename):
        booleans, geom_optimization, stats = read_dftb_stdout(filename)
        for i in booleans:
            summary += i.ljust(30) + str(booleans[i]) + '\n'
        for i in geom_optimization:
            summary += i.ljust(30) + str(geom_optimization[i]) + '\n'
        for i in stats:
            summary += i.ljust(30) + str(stats[i]) + '\n'
        ret = stats['ion_convergence']
    else:
        ret = False
    if ret:
        summary += 'SUCCESS\n'
    else:
        summary += 'FAILED\n'
    logging.info(summary)
    return ''


def update(basedir, population, imember):
    workdir = basedir + os.sep + str(imember)
    filename = workdir+os.sep+'detailed.out'
    if os.path.isfile(filename):
        forces, stress, total_energy = read_detailed_out(filename=filename)
        geometry = workdir+os.sep+'geo_end.gen'
        if forces is not None and stress is not None and total_energy is not None and os.path.isfile(geometry):
            new_structure = read_geometry_gen(geometry)
            properties = {'forces': generic_serializer(forces), 'stress': generic_serializer(stress), 'energy': total_energy}
            population.update_entry(imember=imember, structure=new_structure, properties=properties)


class DFTBplusEvaluator(Evaluator):

    def __init__(self, workdir, slater_path, target_forces, nparal):

        self.workdir = workdir
        if not os.path.isdir(workdir):
            os.mkdir(workdir)
        self.slater_path = slater_path
        self.target_forces = target_forces
        self.population = None
        self.nparal = nparal
        self._running = False

    def initialize(self, population):

        self.population = population

    def run(self):
        self._running = True

        actives_to_evaluate = self.population.actives
        setups = []
        for i in actives_to_evaluate:
            setup = {}
            setup['workdir'] = self.workdir + os.sep + str(i)
            setup['population'] = self.population
            setup['slater_path'] = self.slater_path
            setup['target_forces'] = self.target_forces
            setup['id'] = i
            setups.append(setup)

        #pool = Pool(processes=self.nparal)
        #pool.map(f, setups, chunksize=1)

        procs = []
        index = 0
        for i in range(self.nparal):
            p = Process(target=f, args=(setups[index],))
            p.start()
            procs.append(p)
            index += 1

        while index < len(actives_to_evaluate):
            for p in procs:
                if not p.is_alive():
                    procs.remove(p)
                    p = Process(target=f, args=(setups[index],))
                    p.start()
                    procs.append(p)
                    index += 1

        for i in actives_to_evaluate:
            self.update(i)
        self._running = False

    @property
    def is_running(self):
        return self._running

    def stop(self):
        pass

    def update(self, i):
        update(self.workdir, population=self.population, imember=i)

    def update_db_from_files(self):

        ids = os.listdir(self.workdir)
        for i in ids:
            obid = bson.ObjectId(i)
            if obid in self.population.members:
                self.update(obid)
