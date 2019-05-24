# cython: language_level=3

import numpy as np


def lj_forces(pos, sigmas, epsilons, cp=0.0):
    #cdef int n, i, j
    #cdef float distance, magnitude
    pos = np.array(pos).reshape((-1, 3))
    n = len(pos)
    ret = np.zeros((n, 3))
    for i in range(n - 1):
        for j in range(i + 1, n):
            vector = pos[j] - pos[i]
            distance = np.linalg.norm(vector)
            if distance < 1E-8:
                print('ERROR: Too small distance between atoms')
                distance = 1E-8
            sr6 = (sigmas[i, j] / distance) ** 6
            magnitude = 24 * epsilons[i, j] / distance * (2.0 * sr6 ** 2 - sr6)
            ret[i] -= magnitude * vector / distance
            ret[j] += magnitude * vector / distance
    if cp > 0.0:
        for i in range(n):
            ret[i] -= cp * pos[i]
    return ret.flatten()


def lj_gradient(pos, sigmas, epsilons, cp=0.0):
    return -lj_forces(pos, sigmas, epsilons, cp)


def lj_energy(pos, sigmas, epsilons, cp=0.0):
    #cdef int n, i, j
    #cdef float distance, sr6, ret
    pos = np.array(pos).reshape((-1, 3))
    ret = 0
    n = len(pos)
    for i in range(n - 1):
        for j in range(i + 1, n):
            vector = pos[j] - pos[i]
            distance = np.linalg.norm(vector)
            if distance > 1E-10:
                sr6 = (sigmas[i, j] / distance) ** 6
                ret += 4 * epsilons[i, j] * (sr6 ** 2 - sr6)
    if cp > 0.0:
        for i in range(n):
            ret += 0.5 * cp * np.linalg.norm(pos[i]) ** 2
    return ret

