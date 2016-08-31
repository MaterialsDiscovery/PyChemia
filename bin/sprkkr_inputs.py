#!/usr/bin/env python

import pychemia

st = pychemia.code.vasp.read_poscar()
sk = pychemia.code.sprkkr.CodeSPRKKR(st)
sk.write_sys()
sk.write_pot()
