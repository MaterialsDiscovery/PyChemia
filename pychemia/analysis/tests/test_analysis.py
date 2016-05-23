import pychemia


def notest_analysis():
    """
    Tests (pychemia.analysis.analysis)                           :
    """
    st=pychemia.samples.Al2O3()
    sa=pychemia.analysis.StructureAnalysis(st)
    struc_dist_x, fp_oganov=sa.fp_oganov(delta=0.1, sigma=0.1)
    assert len(fp_oganov)==3

def notest_match():
    """
    Tests (pychemia.analysis.match)                              :
    """
    st = pychemia.samples.Al2O3()
    st2 = st.supercell((2, 3, 4))
    sm = pychemia.analysis.StructureMatch(st, st2)
    sm.match_size()
    sm.match_atoms()
    sm.amtch_shape()
    assert sm.structure1.natom == sm.structure2.natom

