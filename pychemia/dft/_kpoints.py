"""
Definition of the mesh in reciprocal space
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2012"
__version__ = "1.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "March 16, 2014"

import numpy as _np


class KPoints():
    """
    Defines an object that contains information
    about kpoints mesh
    """

    def __init__(self, kmode, comment=None, kpts=None, wgts=None, shifts=None, grid=None):
        """
        Creates a new object kpoints
        if no arguments are provided
        the default is to create one single
        kpoint in gamma.
        The arguments evaluated are

        kpts : Numpy array of kpoint positions
        wgts : Weights for all kpoints
        """
        if comment is None:
            self.comment = ''
        else:
            self.comment = comment.strip()

        self.kmode = kmode.lower()
        self.nkpt = 0

        self.shifts = None
        self.kpts = None
        self.wgts = None
        self.grid = None

        if kpts is None and grid is None:
            self.nkpt = 1
            self.kpts = _np.array([[0, 0, 0]])
            self.wgts = _np.array([1])
        elif kpts is not None and grid is not None:
            raise ValueError("Both 'kpts' and 'grid' cannot be entered simultaneously")
        elif kpts is not None:
            self.kpts = _np.array(kpts)
            self.nkpt = len(self.kpts)
            if wgts is not None and len(wgts) == self.nkpt:
                self.wgts = _np.array(wgts)
            else:
                self.wgts = _np.ones(self.nkpt)
        elif grid is not None:
            self.grid = _np.array(grid)
            if shifts is not None and len(shifts) == 3:
                self.shifts = _np.array(shifts)
            else:
                self.shifts = _np.zeros(3)

    def add_kpt(self, pos, wgt):

        if self.nkpt == 0:
            self.kpts = _np.array(pos).reshape([-1, 3])
            self.wgts = _np.array([wgt]).reshape([-1])
        else:
            self.kpts = _np.append(self.kpts, pos).reshape([-1, 3])
            self.wgts = _np.append(self.wgts, wgt)

        self.nkpt += 1

    def set_grid(self, grid, shifts=None):
        self.grid = _np.array(grid)
        self.nkpt = 0
        if shifts is not None and len(shifts) == 3:
            self.shifts = _np.array(shifts)
        else:
            self.shifts = _np.zeros(3)

    def __str__(self):
        """
        String representation of the kpoints object
        """
        kp = ' Mode  : '+self.kmode+'\n'
        if self.kmode == 'gamma' or self.kmode == 'monkhorst-pack':
            kp += ' Grid  : '+str(self.grid)+'\n'
            kp += ' Shift : '+str(self.shifts)+'\n'
        elif self.kmode == 'cartesian' or self.kmode == 'reciprocal':
            kp = str(self.nkpt) + '\n\n'
            for i in range(self.nkpt):
                kp += (" %15.7f %15.7f %15.7f %20.7f\n"
                       % (self.kpts[i, 0],
                       self.kpts[i, 1],
                       self.kpts[i, 2],
                       self.wgts[i]))
        return kp
