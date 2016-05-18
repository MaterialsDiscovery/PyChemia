import pychemia


def test_analysis():
    """
    Tests (pychemia.analysis.analysis)                           :
    """
    st=pychemia.samples.Al2O3()
    sa=pychemia.analysis.StructureAnalysis(st)
    struc_dist_x, fp_oganov=sa.fp_oganov(delta=0.1, sigma=0.1)
    assert len(fp_oganov)==3
