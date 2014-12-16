__author__ = 'Guillermo Avendano Franco'

import os

from _abifiles import AbiFiles
from _input import InputVariables


class AbinitJob():

    def __init__(self, structure=None, workdir=None, exchange='LDA', kind='FHI'):

        self.structure = structure
        self.workdir = workdir
        self.inp = InputVariables()
        self.inp.from_structure(self.structure)
        self.abifile = AbiFiles(self.workdir)
        self.abifile.set_input(self.inp)
        self.abifile.set_psps(exchange=exchange, kind=kind)

    def _check_workdir(self):

        if self.workdir is None:
            raise ValueError("A proper working directory has not been setup")
        elif not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)

    def write_abiinput(self):
        self.inp.write(self.workdir+os.sep+'abinit.in')

    def write_abifiles(self):
        self.abifile.create()

    def write_all(self):
        self.write_abifiles()
        self.write_abiinput()

    def set_kpoints(self, kp):
        self.inp.set_value('ngkpt', list(kp.grid))

    def ion_relax(self, internal_opt=True, external_opt=True, ionmov=2, nstep=20, ntime=30,
                  tolmxf=1E-7, tolrff=1E-3, dilatmx=1.05):

        if internal_opt and external_opt:
            self.inp.set_value('ndtset', 2)
            self.inp.set_value('optcell', 0, idtset=1)
            self.inp.set_value('optcell', 2, idtset=2)
            self.inp.set_value('ecutsm', 0.5, idtset=2)
        elif external_opt:
            self.inp.set_value('optcell', 2)

        self.inp.set_value('ionmov', ionmov)
        self.inp.set_value('nstep', nstep)
        self.inp.set_value('ntime', ntime)
        self.inp.set_value('tolmxf', tolmxf)
        self.inp.set_value('tolrff', tolrff)
        self.inp.set_value('dilatmx', dilatmx)

    def energetics(self, ecut=50):
        self.inp.set_value('ecut', ecut)
