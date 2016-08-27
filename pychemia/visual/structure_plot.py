from mayavi import mlab
from pychemia.utils.periodic import covalent_radius
import numpy as np


class StructurePlot:

    def __init__(self, structure):
        """

        :param structure: (Structure) A pychemia Structure object
        """
        self.structure = structure

    def plot(self, figname=None, size=(300, 325), view=(30, 30), color=(1.0, 1.0, 1.0)):

        fig = mlab.figure(size=size)
        figure = mlab.gcf()
        fig.scene.disable_render = True
        figure.scene.background = (0.0, 0.0, 0.0)
        mlab.view(0, 90, distance=0.2)
        assert (self.structure.natom > 0)

        x = self.structure.positions[:, 0]
        y = self.structure.positions[:, 1]
        z = self.structure.positions[:, 2]
        cr = covalent_radius(self.structure.symbols)
        s = np.apply_along_axis(np.linalg.norm, 1, self.structure.positions)

        mlab.points3d(x, y, z, s, scale_factor=1.0, resolution=8, opacity=1.0,
                      color=color,
                      scale_mode='none')

        if self.structure.is_crystal:
            frame, line1, line2, line3 = self.structure.get_cell().get_path()

            mlab.plot3d(frame[:, 0], frame[:, 1], frame[:, 2], tube_radius=.02, color=(1, 1, 1))
            mlab.plot3d(line1[:, 0], line1[:, 1], line1[:, 2], tube_radius=.02, color=(1, 1, 1))
            mlab.plot3d(line2[:, 0], line2[:, 1], line2[:, 2], tube_radius=.02, color=(1, 1, 1))
            mlab.plot3d(line3[:, 0], line3[:, 1], line3[:, 2], tube_radius=.02, color=(1, 1, 1))
        else:
            for i in range(self.structure.natom - 1):
                for j in range(i + 1, self.structure.natom):
                    vector = self.structure.positions[i] - self.structure.positions[j]
                    mvector = np.linalg.norm(vector)
                    uvector = 1.0 / mvector * vector
                    if 2 * mvector < covalent_radius(self.structure.symbols[i]) + \
                            covalent_radius(self.structure.symbols[j]):
                        pair = np.concatenate(
                            (self.structure.positions[i] - 0.1 * uvector,
                             self.structure.positions[j] + 0.1 * uvector)).reshape((-1, 3))
                        mlab.plot3d(pair[:, 0], pair[:, 1], pair[:, 2], tube_radius=0.15, opacity=1.0, color=(1, 1, 1))

        mlab.view(distance=12.0)
        fig.scene.disable_render = False
        if figname is not None:
            mlab.savefig(figname)
        return figure
