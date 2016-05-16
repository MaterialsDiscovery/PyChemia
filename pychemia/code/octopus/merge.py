import os as _os
import shutil as _shutil
import numpy as _np


class OctopusStates:
    def __init__(self, basedir):
        rf = open(basedir + '/restart/gs/states')

        self.basedir = basedir
        self.nst = int(rf.readline().split()[1])
        self.dim = int(rf.readline().split()[1])
        self.nik = int(rf.readline().split()[1])

    def __str__(self):
        return "nst= %25d\ndim= %25d\nnik= %25d" % (self.nst, self.dim, self.nik)

    def __repr__(self):
        return "states('%s')" % self.basedir

    def write(self, filename):
        wf = open(filename, 'w')
        wf.write(self.__str__())
        wf.close()


class OctopusOccupancies:
    def __init__(self, basedir):
        rf = open(basedir + '/restart/gs/occs')

        self.basedir = basedir
        self.header = rf.readline() + rf.readline()
        rf.close()

        self.data = _np.loadtxt(basedir + '/restart/gs/occs', ndmin=1, skiprows=2, delimiter=' | ', comments='%',
                                dtype={'names': ['occupations', 'eigenvalue[a.u.]', 'k-points X', 'k-points Y',
                                                 'k-points Z', 'k-weights', 'filename', 'ik', 'ist', 'idim'],
                                       'formats': ['f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'S10', 'i4', 'i4', 'i4']})

    def __repr__(self):
        return "occs('%s')" % self.basedir

    def __str__(self):
        text = self.header
        for line in self.data:
            text += "%21.14e | %21.14e | %21.14e | %21.14e | %21.14e | %21.14e | %10s | %3d | %3d | %3d\n" % tuple(line)
        text += "%"
        return text

    def write(self, filename):
        wf = open(filename, 'w')
        wf.write(self.__str__())
        wf.close()


class OctopusWavefunctions:
    def __init__(self, basedir):

        rf = open(basedir + '/restart/gs/wfns')

        self.basedir = basedir
        self.header = rf.readline() + rf.readline()
        rls = rf.readlines()
        lines = [x.split(' | ') for x in rls]
        nlines = len(lines) - 2
        rf.close()

        self.data = _np.zeros((nlines, 3))
        self.filenames = []
        for i in range(nlines):
            self.data[i, 0] = int(lines[i][0])
            self.data[i, 1] = int(lines[i][1])
            self.data[i, 2] = int(lines[i][2])
            self.filenames.append(lines[i][3].strip())
        self.iter = int(lines[-1][0].strip().split()[-1])

    def __repr__(self):
        return "wfns('%s')" % self.basedir

    def __str__(self):
        text = self.header
        for i in range(len(self.data)):
            text += "%8d | %8d | %8d | %12s\n" % tuple(list(self.data[i]) + [self.filenames[i]])
        text = text + '%\n' + ("Iter = %7d" % self.iter)
        return text

    def write(self, filename):
        wf = open(filename, 'w')
        wf.write(self.__str__())
        wf.close()


class OctopusStatesFile:
    def __init__(self, basedir):
        self.basedir = basedir
        self.obfs = sorted([x for x in _os.listdir(basedir + '/restart/gs')
                            if _os.path.isfile(basedir + '/restart/gs/' + x) and x[-4:] == '.obf' and x[0] == '0'])

    def __repr__(self):
        return "states_obf('%s')" % self.basedir


class OctopusRestart:
    def __init__(self, basedir):
        self.states = OctopusStates(basedir)
        self.occs = OctopusOccupancies(basedir)
        self.wfns = OctopusWavefunctions(basedir)
        self.states_obf = OctopusStatesFile(basedir)


def oct_merge(res1, res2, basedir, mksym1, mksym2):
    for i in [basedir, basedir + '/restart', basedir + '/restart/gs']:
        if not _os.path.isdir(i):
            _os.mkdir(i)

    res1.states.nst += res2.states.nst
    assert res1.states.dim == res2.states.dim
    assert res1.states.nik == res2.states.nik
    res1.states.write(basedir + '/restart/gs/states')

    lastf = int(res1.occs.data[-1][6])  # last filename in res1
    lasts = int(res1.occs.data[-1][8])  # last state number in res1
    lasts = res1.wfns.data[-1][1]  # last state number in res1
    lastf = int(res1.wfns.filenames[-1][1:-1])  # last filename in res1
    lastfilename = lastf
    index = 0
    for i in range(len(res2.occs.data)):
        res2.occs.data[i][6] = str(lastf + 1).zfill(10)
        res2.occs.data[i][8] = lasts + 1
        res2.wfns.data[i][1] = lasts + 1
        res2.wfns.filenames[i] = '"' + str(lastf + 1).zfill(10) + '"'
        lastf += 1
        if index == 0:
            index = 1
        else:
            index = 0
            lasts += 1
    res1.occs.data = _np.concatenate((res1.occs.data, res2.occs.data), axis=0)
    res1.occs.write(basedir + '/restart/gs/occs')
    res1.wfns.data = _np.concatenate((res1.wfns.data, res2.wfns.data), axis=0)
    res1.wfns.filenames += res2.wfns.filenames
    res1.wfns.write(basedir + '/restart/gs/wfns')

    # Copy or make symlinks of the States (*.obf) from first directory
    for i in res1.states_obf.obfs:
        print('FIRST: Copying or Linking ', i)
        if mksym1:
            _os.symlink(_os.getcwd() + '/' + res1.states_obf.basedir + '/restart/gs/' + i, basedir + '/restart/gs/' + i)
        else:
            _shutil.copy2(res1.states_obf.basedir + '/restart/gs/' + i, basedir + '/restart/gs')

    # Copy or make symlinks of the States (*.obf) from second directory
    for i in res2.states_obf.obfs:
        print('SECOND: Copying or Linking ', i)
        if mksym2:
            _os.symlink(_os.getcwd() + '/' + res2.states_obf.basedir + '/restart/gs/' + i,
                        basedir + '/restart/gs/' + str(lastfilename + 1).zfill(10) + '.obf')
        else:
            _shutil.copyfile(res2.states_obf.basedir + '/restart/gs/' + i,
                             basedir + '/restart/gs/' + str(lastfilename + 1).zfill(10) + '.obf')
        lastfilename += 1

    dirs = [x for x in _os.listdir(res1.states.basedir + '/restart/gs/') if
            x[0] != '0' and x not in ['states', 'occs', 'wfns']]

    # Copy the other files
    print(dirs)
    for i in dirs:
        _shutil.copy2(res1.states.basedir + '/restart/gs/' + i, basedir + '/restart/gs')


def merge_restart(dir1, dir2, dest, mksym1=True, mksym2=False):
    """
    Merge the restart directories (dir1 and dir2)
    into the destination (dest)
    """
    res1 = OctopusRestart(dir1)
    res2 = OctopusRestart(dir2)
    oct_merge(res1, res2, dest, mksym1, mksym2)
