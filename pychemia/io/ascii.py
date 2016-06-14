"""
Module for support V_sim ascii fileformat
Contains routines to load .ascii files and
create pychemia Structure objects and to save
back to .ascii files

This code was originally created for ASE
"""

import re as _re
from pychemia.utils.constants import bohr_angstrom
from pychemia.core import Structure


def load(filep):
    """
    Read an V_sim .ascii file and returns a pychemia
    Structure object

    Args:
       filep: (string) Path to a .ascii file or an
             actual file-like object

    Returns:
       struct: (object) A pychemia Structure object
    """

    if isinstance(filep, str):
        f = open(filep)
    else:  # Assume it's a file-like object
        f = filep

    comment = f.readline()

    line = f.readline() + ' ' + f.readline()
    box = line.split()
    for i in range(len(box)):
        box[i] = float(box[i])

    keywords = []
    positions = []
    symbols = []

    re_comment = _re.compile('^\s*[#!]')
    re_node = _re.compile('^\s*\S+\s+\S+\s+\S+\s+\S+')

    while True:
        line = f.readline()

        if line == '':
            break  # EOF -> Exit

        p = re_comment.match(line)
        if p is not None:
            # remove comment character at the beginning of line
            line = line[p.end():].replace(',', ' ').lower()
            if line[:8] == "keyword:":
                keywords.extend(line[8:].split())

        elif re_node.match(line):
            unit = 1.0
            if not ("reduced" in keywords):
                if ("bohr" in keywords) or ("bohrd0" in keywords) or ("atomic" in keywords) or ("atomicd0" in keywords):
                    unit = bohr_angstrom

            fields = line.split()
            positions.append([unit * float(fields[0]),
                              unit * float(fields[1]),
                              unit * float(fields[2])])
            symbols.append(fields[3])

    f.close()

    if ("surface" in keywords) or ("freeBC" in keywords):
        raise NotImplementedError

    # create atoms object based on the information
    if "angdeg" in keywords:
        cell = cell.par_to_cell(box)
    else:
        unit = 1.0
        if ("bohr" in keywords) or ("bohrd0" in keywords) or ("atomic" in keywords) or ("atomicd0" in keywords):
            unit = bohr_angstrom
        cell = [[unit * box[0], 0.0, 0.0],
                [unit * box[1], unit * box[2], 0.0],
                [unit * box[3], unit * box[4], unit * box[5]]]

    if "reduced" in keywords:
        struct = Structure(cell=cell, reduced=positions, symbols=symbols, name=comment)
    else:
        struct = Structure(cell=cell, positions=positions, symbols=symbols, name=comment)

    return struct


def save(struct, filep, cartesian=True, long_format=True, angdeg=False):
    """
    Saves a pychemia Structure object in  V_sim .ascii fileformat
    in the simplest way, i.e. using all
    defaults with no optional keywords. In the first line we add the
    number of atoms, as this is used by certain code
    """

    if isinstance(filep, str):
        f = open(filep, 'w')
    else:  # Assume it's a 'file-like object'
        f = filep

    # write header (treated as a comment by v_sim
    f.write("%s\n" % struct.name)

    # write cell
    cell = struct.cell
    if angdeg:
        ddd = cell_to_cellpar(cell)
    else:
        ddd = cell_to_reduced(cell)

    f.write("%.14f %.14f %.14f\n" % (ddd[0], ddd[1], ddd[2]))
    f.write("%.14f %.14f %.14f\n" % (ddd[3], ddd[4], ddd[5]))

    if angdeg:
        f.write("#keyword: angdeg\n")

    # Write atom positions in scaled or cartesian coordinates
    if cartesian:
        coord = struct.positions
    else:
        f.write("#keyword: reduced\n")
        coord = struct.reduced

    if long_format:
        cform = ' %19.16f'
    else:
        cform = ' %9.6f'

    symbols = struct.symbols
    for iatom, atom in enumerate(coord):
        f.write(' ')
        for dcoord in atom:
            f.write(cform % dcoord)
        f.write('  ' + symbols[iatom] + '\n')


def cell_to_reduced(full):
    """
    Transforms the given matrix full into a reduced array used by
    V_Sim to store box definition.

    translated from src/coreTools/toolMatrix.c
    subroutine tool_matrix_reducePrimitiveVectors
    """

    from numpy import zeros
    from numpy.linalg import norm

    xcoord = full[0].copy()

    # Compute the Y vector of the new basis, orthogonal to xcoord an coplanar with xcoord and old y vector
    u = zeros(3)
    u[0] = full[0][1] * full[1][2] - full[0][2] * full[1][1]
    u[1] = full[0][2] * full[1][0] - full[0][0] * full[1][2]
    u[2] = full[0][0] * full[1][1] - full[0][1] * full[1][0]

    deltaij = xcoord[0] * u[1] - xcoord[1] * u[0]
    if deltaij != 0.0:
        i = 0
        j = 1
        k = 2
    else:
        deltaij = xcoord[0] * u[2] - xcoord[2] * u[0]
        if deltaij != 0.0:
            i = 0
            j = 2
            k = 1
        else:
            deltaij = xcoord[1] * u[2] - xcoord[2] * u[1]
            if deltaij != 0.0:
                i = 1
                j = 2
                k = 0
            else:
                # Error
                return None

    y = zeros(3)
    y[k] = -1.0
    y[i] = (xcoord[k] * u[j] - xcoord[j] * u[k]) / deltaij
    y[j] = (xcoord[i] * u[k] - xcoord[k] * u[i]) / deltaij

    # We need to turn Y if y.Y is negative
    fnorm = full[1][0] * y[0] + full[1][1] * y[1] + full[1][2] * y[2]
    if fnorm < 0.0:
        y *= -1.

    # Compute the new Z vector in order to form a direct orthogonal
    # basis with xcoord and Y
    z = zeros(3)
    z[0] = xcoord[1] * y[2] - xcoord[2] * y[1]
    z[1] = xcoord[2] * y[0] - xcoord[0] * y[2]
    z[2] = xcoord[0] * y[1] - xcoord[1] * y[0]

    # Normalize vectors
    xcoord /= norm(xcoord)
    y /= norm(y)
    z /= norm(z)

    # Compute the reduce value for the basis.
    reduced = zeros(6)
    reduced[0] = xcoord[0] * full[0][0] + xcoord[1] * full[0][1] + xcoord[2] * full[0][2]
    reduced[1] = xcoord[0] * full[1][0] + xcoord[1] * full[1][1] + xcoord[2] * full[1][2]
    reduced[2] = y[0] * full[1][0] + y[1] * full[1][1] + y[2] * full[1][2]
    reduced[3] = xcoord[0] * full[2][0] + xcoord[1] * full[2][1] + xcoord[2] * full[2][2]
    reduced[4] = y[0] * full[2][0] + y[1] * full[2][1] + y[2] * full[2][2]
    reduced[5] = z[0] * full[2][0] + z[1] * full[2][1] + z[2] * full[2][2]

    return reduced
