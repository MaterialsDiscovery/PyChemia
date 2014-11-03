__author__ = 'Guillermo Avendano Franco'

import os
import shutil
import pychemia


class RelaxPopulation():

    def __init__(self, population, basedir):
        self.population = population
        self.basedir = basedir
        self.AbinitJobs = {}

    def create_dirs(self, clean=False):
        if not os.path.isdir(self.basedir):
            os.makedirs(self.basedir)
        elif clean:
            for i in os.listdir(self.basedir):
                path = self.basedir+os.sep+i
                if os.path.isfile(path):
                    os.remove(path)
                elif os.path.isdir(path):
                    shutil.rmtree(path)
        for i in self.population.pcdb.entries.find():
            name = self.basedir+os.sep+str(i['_id'])
            if not os.path.isdir(name):
                os.mkdir(name)

    def create_inputs(self):
        kpoints = pychemia.dft.KPoints(kmode='gamma', grid=[4, 4, 4])
        for i in self.population.pcdb.entries.find():
            name = str(i['_id'])
            workdir = self.basedir+os.sep+name
            struct = pychemia.Structure().fromdict(i)
            job = pychemia.code.abinit.AbinitJob(struct, workdir)
            job.set_kpoints(kpoints)
            job.ion_relax(tolmxf=1E-4, tolrff=1E-2,)
            job.energetics()
            job.write_all()
            self.AbinitJobs[name] = job

