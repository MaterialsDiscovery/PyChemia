__author__ = 'Guillermo Avendano Franco'


class Tasks():

    def __init__(self, abinput):
        """
        Create a ABINIT tasks object completing information in the input for running
        specific calculations using ABINIT

        :param abinput:
        """
        self.abinput = abinput

    def geometry_opt(self, internal_opt=True, external_opt=True, ionmov=2, nstep=20, ntime=30,
                     tolmxf=1E-7, tolrff=1E-3):

        if internal_opt and external_opt:
            self.abinput.set_value('optcell', 0, idtset=1)
            self.abinput.set_value('optcell', 2, idtset=2)
        elif external_opt:
            self.abinput.set_value('optcell', 2)

        self.abinput.set_value('ionmov', ionmov)
        self.abinput.set_value('nstep', nstep)
        self.abinput.set_value('ntime', ntime)
        self.abinput.set_value('tolmxf', tolmxf)
        self.abinput.set_value('tolrff', tolrff)

    def kpoints(self, ngkpt=(3, 3, 3)):
        self.abinput.set_value('ngkpt', list(ngkpt))

    def energetics(self, ecut=50):
        self.abinput.set_value('ecut', ecut)

