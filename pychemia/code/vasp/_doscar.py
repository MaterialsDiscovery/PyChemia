import os
import numpy as np


class VaspDoscar:
    def __init__(self, filename='DOSCAR'):

        if not os.path.isfile(filename):
            raise ValueError('ERROR: DOSCAR file not found')

        self.filename = filename

        self.doscar = self.parse_doscar(self.filename)

        if 'projected' in self.doscar:
            self.has_projected = True
            self.nions = len(self.doscar['projected'])
        else:
            self.has_projected = False
            self.nions = None

        self.ndos, self.total_ncols = self.doscar['total'].shape

        if self.total_ncols == 5:
            self.is_spin_polarized = True
        else:
            self.is_spin_polarized = False

        self.total_dos = self._dos_dict()

        if self.has_projected:
            self.projected_dos = self._proj_dos_dict()

    def _dos_dict(self):

        ret = {'energy': self.doscar['total'][:, 0]}

        if self.is_spin_polarized:
            ret['dos_up'] = self.doscar['total'][:, 1]
            ret['dos_down'] = self.doscar['total'][:, 2]
            ret['integrated_dos_up'] = self.doscar['total'][:, 3]
            ret['integrated_dos_down'] = self.doscar['total'][:, 4]
        else:
            ret['dos'] = self.doscar['total'][:, 1]
            ret['integrated_dos'] = self.doscar['total'][:, 2]
        return ret

    def _proj_dos_dict(self):

        ret = {}
        for i in range(self.nions):
            ret['proj'] = {'energy': self.doscar['projected'][i][:, 0]}
            if self.is_spin_polarized:
                ret['ncols'] = self.doscar['projected'][i].shape[1]
        return ret

    @staticmethod
    def parse_doscar(filename):

        if not os.path.isfile(filename):
            raise ValueError('ERROR: DOSCAR file not found')

        rf = open(filename)
        data = rf.readlines()
        rf.close()

        if len(data) < 5:
            raise ValueError('DOSCAR seems truncated')

        # Skipping the first lines of header
        iline = 5

        header = [float(x) for x in data[iline].strip().split()]
        ndos = int(header[2])
        iline += 1

        total_dos = [[float(x) for x in y.split()] for y in data[iline:iline + ndos]]
        total_dos = np.array(total_dos)

        iline += ndos

        # In case there are more lines of data, they are the projected DOS
        if len(data) > iline:

            projected_dos = []

            while iline < len(data):
                header = [float(x) for x in data[iline].strip().split()]
                ndos = int(header[2])
                iline += 1
                tmp_dos = [[float(x) for x in y.split()] for y in data[iline:iline + ndos]]
                projected_dos.append(tmp_dos)
                iline += ndos

            projected_dos = np.array(projected_dos)

            return {'total': total_dos, 'projected': projected_dos}

        else:

            return {'total': total_dos}
