import pychemia
import os
from .samples import Al2O3


def test_analysis():
    """
    Test (pychemia.analysis.analysis)                           :
    """
    st = Al2O3()
    sa = pychemia.analysis.StructureAnalysis(st, radius=20)
    struc_dist_x, fp_oganov = sa.fp_oganov(delta=0.1, sigma=0.1)
    assert len(fp_oganov) == 3


def test_match():
    """
    Test (pychemia.analysis.match)                              :
    """
    st = Al2O3()
    st2 = st.supercell((2, 3, 4))
    sm = pychemia.analysis.StructureMatch(st, st2)
    sm.match_size()
    sm.match_shape()
    assert sm.structure1.natom == sm.structure2.natom


def test_distances():
    """
    Test (pychemia.analysis.match)                              :
    """
    print(os.getcwd())
    st = pychemia.io.xyz.load('tests/data/xyz/chlorophyll.xyz')
    sa = pychemia.analysis.StructureAnalysis(st)
    distances = sa.all_distances()
    assert len(distances) == int(st.natom*(st.natom+1)/2)

