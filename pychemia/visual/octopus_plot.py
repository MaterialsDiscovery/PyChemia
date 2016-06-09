import os
import numpy as np


def plot_energy(path):
    import matplotlib.pylab as plt
    if not os.path.isdir(path):
        print('ERROR: No such dir ', path)
        return None
    if not os.path.isfile(path + '/td.general/energy'):
        print('ERROR: No energy file ', path + '/td.general/energy')
        return None

    data = np.loadtxt(path + '/td.general/energy')

    if len(data) == 0:
        print('No data on ', path + '/td.general/energy')
        return None

    plt.clf()
    fig = plt.figure(num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
    labels = ['Total', 'Kinetic (ions)', 'Ion-Ion', 'Electronic',
              'Eigenvalues', 'Hartree', 'Int[n v_xc]', 'Exchange', 'Correlation']

    # For some columns the energy is only computed each
    # certain steps, this will identify the number of
    # steps to get a new non-zero value
    nstep = 0
    for i in range(len(data)):
        if i == 0 and data[i, -3] != 0.0:
            nstep += 1
        if data[i, -3] == 0.0:
            nstep += 1
        if i != 0 and data[i, -3] != 0.0:
            break
    print('Complete data each', nstep, ' steps')

    for i in range(9):
        plt.subplot(330 + i + 1)
        plt.plot(data[::nstep, 1], data[::nstep, i + 2], label=labels[i])
        plt.xlabel('Time [hbar/H]')
        plt.legend(loc=4)

    return fig
