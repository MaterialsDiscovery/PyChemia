import numpy as np
import itertools
from tvtk.api import tvtk
from pychemia import pcm_log
from mayavi import mlab


class LatticePlot:
    def __init__(self, lattice):
        self.lattice = lattice

    def plot(self, points=None):

        tube_rad = max(self.lattice.lengths) / 100.0

        frame, line1, line2, line3 = self.lattice.get_path()
        for i, j, k in [[frame[:, 0], frame[:, 1], frame[:, 2]],
                        [line1[:, 0], line1[:, 1], line1[:, 2]],
                        [line2[:, 0], line2[:, 1], line2[:, 2]],
                        [line3[:, 0], line3[:, 1], line3[:, 2]]]:
            mlab.plot3d(i, j, k, tube_radius=tube_rad, color=(1, 1, 1), tube_sides=24, transparent=True, opacity=0.5)

        if points is not None:
            ip = np.array(points)
            mlab.points3d(ip[:, 0], ip[:, 1], ip[:, 2], tube_rad * np.ones(len(ip)), scale_factor=1)

        return mlab.pipeline

    def plot_wigner_seitz(self, scale=1):

        ws = np.array(self.lattice.get_wigner_seitz())
        points = np.array([])
        for i in ws:
            for j in i:
                points = np.concatenate((points, j))
        points = scale * (points.reshape(-1, 3))

        index = 0
        # triangles = index + np.array(list(itertools.combinations(range(len(ws[0])), 3)))
        # scalars = _np.ones(len(ws[0]))
        # for i in ws[1:]:
        #     index += len(i)
        #     triangles = np.concatenate((triangles, index + _np.array(list(itertools.combinations(range(len(i)), 3)))))
        #     scalars = np.concatenate((scalars, _np.random.random() * _np.ones(len(i))))
        # scalars = _np.ones(len(ws[0]))
        scalars = None
        triangles = None
        for i in ws:
            pcm_log.debug(i)
            iscalars = np.ones(len(i))
            itriangles = index + np.array(list(itertools.combinations(range(len(i)), 3)))
            pcm_log.debug(iscalars)
            pcm_log.debug(itriangles)
            if triangles is None:
                triangles = itriangles
            else:
                triangles = np.concatenate((triangles, itriangles))
            if scalars is None:
                scalars = iscalars
            else:
                scalars = np.concatenate((scalars, iscalars))
            index += len(i)

        # print triangles
        # print scalars
        # The TVTK dataset.
        mesh = tvtk.PolyData(points=points, polys=triangles)
        mesh.point_data.scalars = scalars
        mesh.point_data.scalars.name = 'scalars'

        pipeline = self.plot()
        pipeline.surface(mesh, color=(0.9, 0.1, 0.1), opacity=0.2)
        return pipeline
