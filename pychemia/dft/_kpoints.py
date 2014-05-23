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

    def __init__(self, **kwargs):
        """
        Creates a new object kpoints
        if no arguments are provided
        the default is to create one single
        kpoint in gamma.
        The arguments evaluated are

        kpts : Numpy array of kpoint positions
        wgts : Weights for all kpoints
        """
        self.comment = ''

        # Default to a single kpoint at gamma
        if len(kwargs) == 0:
            self.nkpt = 1
            self.kpts = _np.array([[0, 0, 0]])
            self.wgts = _np.array([1])

        # Setting the options provided by kwargs
        if 'nkpt' in kwargs:
            if kwargs['nkpt'] == 0 or kwargs['nkpt'] is None:
                # Create an empty kpoints object
                self.nkpt = 0
                self.kpts = _np.array([[]])
                self.wgts = _np.array([])
            else:
                self.nkpt = kwargs['nkpt']

        if 'kpts' in kwargs:
            self.kpts = _np.array(kwargs['kpts'])

        if 'wgts' in kwargs:
            self.wgts = _np.array(kwargs['wgts'])

        # Autocompleting missing data
        if self.nkpt != len(self.kpts):
            self.nkpt = len(self.kpts)

        if len(self.wgts) != len(self.kpts):
            self.wgts = _np.ones(self.nkpt)

    def add_kpt(self, pos, wgt):

        if self.nkpt == 0:
            self.kpts = _np.array(pos).reshape([-1, 3])
            self.wgts = _np.array([wgt]).reshape([-1])
        else:
            self.kpts = _np.append(self.kpts, pos).reshape([-1, 3])
            self.wgts = _np.append(self.wgts, wgt)

        self.nkpt += 1

    def __str__(self):
        """
        String representation of the kpoints object
        """
        if self.nkpt > 0:
            kp = str(self.nkpt) + '\n\n'
            for i in range(self.nkpt):
                kp += (" %15.7f %15.7f %15.7f %20.7f\n"
                       % (self.kpts[i, 0],
                       self.kpts[i, 1],
                       self.kpts[i, 2],
                       self.wgts[i]))
        else:
            kp = 'Empty KPOINTS'
        return kp
