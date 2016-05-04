import numpy as np
from . import procar


class BandStructure:
    def __init__(self, structure, filename='PROCAR'):
        self.structure = structure
        rcell = self.structure.lattice.reciprocal().cell
        print(rcell)
        self.procarFile = procar.ProcarParser()
        self.procarFile.readFile(procar=filename, permissive=False, recLattice=rcell)

        self.data = procar.ProcarSelect(self.procarFile, deepCopy=True)

    def simple_bands(self):
        self.data.selectIspin([0])
        self.data.selectAtoms([-1])
        self.data.selectOrbital([-1])

        fermi = 0

        bands = (self.data.bands.transpose() - np.array(fermi)).transpose()

        plot = procar.ProcarPlot(bands, self.data.spd, self.data.kpoints)

        fig, ax = plot.scatterPlot()
        ax.set_ylabel(r"Energy [eV]")

        return fig, ax
