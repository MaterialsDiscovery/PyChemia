__author__ = 'Guillermo Avendano-Franco'

import random
from _metaheuristics import MetaHeuristics
from pychemia.sc import Population
from pychemia import Composition
from pychemia.analysis import StructureChanger


class HarmonySearch(MetaHeuristics):

    def __init__(self, name, composition=None, size=None, dbmaster=None, nmaxfromdb=None):

        self.name = name
        self.size = size
        self.hmcr = None
        self.par = None
        self.relax_population = None
        self.runner = None
        if composition is not None:
            comp = Composition(composition)
            self.population = Population(name, new=True)

            nindb = dbmaster.entries.find({'formula': comp.formula}).count()
            if nindb < nmaxfromdb:
                nrandstruct = size - nindb
            else:
                nrandstruct = size - nmaxfromdb

            self.population.add_random(composition, size=nrandstruct)
            self.population.add_from_db(composition, dbmaster.name, sizemax=nmaxfromdb)

        else:
            print 'Reading from already present database'
            self.population = Population(name)
            print 'Number of entries: ', len(self)

    def __len__(self):
        return len(self.population)

    def set_hmcr(self, hmcr):
        if 1.0 >= hmcr >= 0.0:
            self.hmcr = hmcr

    def set_par(self, par):
        if 1.0 >= par >= 0.0:
            self.par = par

    def set_run(self, code, runner, basedir):

        if code == 'vasp':
            from pychemia.code.vasp import RelaxPopulation
            self.relax_population = RelaxPopulation(self.population, basedir)

        self.runner = runner

        self.relax_population.create_dirs(clean=True)
        self.relax_population.create_inputs()

    def run(self, nparal):

        entries_ids = self.population.entries_ids

        def worker(workdir):
            self.runner.run(dirpath=workdir)

        def checker(workdir):
            return self.relax_population.update(workdir)

        workdirs = [self.relax_population.basedir+os.sep+i for i in entries_ids]
        self.runner.run_multidirs(workdirs, worker, checker, nconcurrent=nparal)

        self.search_harmony()

    def search_harmony(self):

        for i in self.population.actives:

            rnd = random.random()
            if rnd < self.hmcr:
                rnd = random.random()
                if rnd < self.par:
                    changer = StructureChanger(self.population.structure[i])
                    newstruct = changer.random_change(self.delta)
