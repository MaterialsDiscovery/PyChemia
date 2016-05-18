#!/usr/bin/env python

import os
import scipy.io
from mayavi import mlab
import numpy as np
import pychemia
import pychemia.io.xyz


class OctopusVisualization:
    def __init__(self, fig, fs, gs, bs):
        """
        Creates a container for the
        sources of a visualization
        """
        self.fig = fig  # Figure
        self.fs = fs  # Field sources
        self.gs = gs  # Atomic Geometry sources
        self.bs = bs  # Cell box sources


def _compute_box(b):
    box = np.array([[0, 0, 0], [b[0], 0, 0], [b[0], b[1], 0],
                    [0, b[1], 0], [0, 0, 0], [0, 0, b[2]],
                    [0, b[1], b[2]], [b[0], b[1], b[2]],
                    [b[0], 0, b[2]], [0, 0, b[2]]])

    x0 = box[:, 0]
    y0 = box[:, 1]
    z0 = box[:, 2]

    x1 = np.array([b[0], b[0]])
    y1 = np.array([b[1], b[1]])
    z1 = np.array([0, b[2]])

    return x0, y0, z0, x1, y1, z1


def _read_netcdf(filename):
    if not os.path.isfile(filename):
        print('ERROR: Non existing file: ' + filename)
        exit(1)
    data = scipy.io.netcdf_file(filename)
    rdata = data.variables['rdata'][:]
    pos = data.variables['pos'][:]
    data.close()
    return np.abs(rdata), pos


def oct_visualize(path='static',
                  spin=None,
                  spin_pol=None,
                  size=(640, 480),
                  bgcolor=(0, 0, 0),
                  structure=True,
                  box=True,
                  figname='figure.png',
                  view=None,
                  visual=None):
    if spin is None:
        spin = [True, True]

    if visual is None:
        fig = mlab.figure(size=size, bgcolor=bgcolor)
    else:
        fig = visual.fig

    if spin_pol is None:
        inp = pychemia.code.octopus.InputVariables('inp')
        if inp.variables['SpinComponents'] == 'spin_polarized':
            spin_pol = True
        else:
            spin_pol = False

    rdata = []
    pos = []
    if not spin_pol:
        r, p = _read_netcdf(path + '/density.ncdf')
        rdata.append(r)
        pos.append(p)
        b = list(np.array(rdata).shape)
    else:
        if spin[0]:
            r, p = _read_netcdf(path + '/density-sp1.ncdf')
            rdata.append(r)
            pos.append(p)
            b = list(r.shape)
            pos = p
        if spin[1]:
            r, p = _read_netcdf(path + '/density-sp2.ncdf')
            rdata.append(r)
            pos.append(p)
            b = list(r.shape)
            pos = p

        if spin[0] is True and spin[1] is True and not np.all(pos[0] == pos[1]):
            print('ERROR box not consistent on', path)
            exit()

    if visual is None:
        fs = []
        if not spin_pol:
            fs1 = mlab.pipeline.scalar_field(np.log(rdata))
            mlab.pipeline.volume(fs1, vmin=-4, vmax=-1)
            fs.append(fs1)
        else:
            if spin[0]:
                fs1 = mlab.pipeline.scalar_field(np.log(rdata1))
                mlab.pipeline.volume(fs1, vmin=-4, vmax=-1, color=(1, 0, 0))
                fs.append(fs1)
            if spin[1]:
                fs2 = mlab.pipeline.scalar_field(np.log(rdata2))
                mlab.pipeline.volume(fs2, vmin=-4, vmax=-1, color=(0, 0, 1))
                fs.append(fs2)
    else:
        count = 0
        if not spin_pol:
            visual.fs[0].mlab_source.scalars = np.log(rdata)
        else:
            if spin[0]:
                visual.fs[count].mlab_source.scalars = np.log(rdata1)
                count += 1
            if spin[1]:
                visual.fs[count].mlab_source.scalars = np.log(rdata2)

    if structure:
        geometry_file = pychemia.io.xyz.load(path + '/geo.xyz')

        orig = pos[:, 0]
        delta = pos[:, 1]

        x = (geometry_file.positions[:, 0] - orig[0]) / delta[0]
        y = (geometry_file.positions[:, 1] - orig[1]) / delta[1]
        z = (geometry_file.positions[:, 2] - orig[2]) / delta[2]
        cr = pychemia.utils.periodic.covalent_radius(geometry_file.symbols)

        if visual is None:
            gs = mlab.points3d(x, y, z, cr, scale_factor=5)
        else:
            visual.gs.mlab_source.set(x=x, y=y, z=z)
    else:
        gs = None

    if visual is None:
        bs = []
        if box:
            x0, y0, z0, x1, y1, z1 = _compute_box(b)
            ze = np.zeros(2)

            bb1 = mlab.plot3d(x0, y0, z0, tube_radius=.15, color=(1, 1, 1))
            bb2 = mlab.plot3d(ze, y1, z1, tube_radius=.15, color=(1, 1, 1))
            bb3 = mlab.plot3d(x1, y1, z1, tube_radius=.15, color=(1, 1, 1))
            bb4 = mlab.plot3d(x1, ze, z1, tube_radius=.15, color=(1, 1, 1))

            bs = [bb1, bb2, bb3, bb4]
    else:
        if box:
            x0, y0, z0, x1, y1, z1 = _compute_box(b)
            ze = np.zeros(2)
            visual.bs[0].start()
            visual.bs[1].start()
            visual.bs[2].start()
            visual.bs[3].start()

            visual.bs[0].mlab_source.set(x=x0, y=y0, z=z0, tube_radius=.15, color=(1, 1, 1))
            visual.bs[1].mlab_source.set(x=ze, y=y1, z=z1, tube_radius=.15, color=(1, 1, 1))
            visual.bs[2].mlab_source.set(x=x1, y=y1, z=z1, tube_radius=.15, color=(1, 1, 1))
            visual.bs[3].mlab_source.set(x=x1, y=ze, z=z1, tube_radius=.15, color=(1, 1, 1))

        else:
            try:
                visual.bs[0].stop()
                visual.bs[1].stop()
                visual.bs[2].stop()
                visual.bs[3].stop()
            except ValueError:
                pass

    if visual is None:
        visual = OctopusVisualization(fig, fs, gs, bs)

    if view is None:
        mlab.view(45, 55, 500, np.array([b[0] / 2, b[1] / 2, b[2] / 2]))
    else:
        mlab.view(view[0], view[1], view[2], view[3])

    mlab.savefig(figname)

    return visual
