import pychemia


def test_abinit_parser():
    """
    Tests (pychemia.code.abinit) [parser]                        :
    """
    if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
        from pychemia.code.abinit import parser
        from numpy import array, all, ones
        from math import sqrt
        import tempfile

        wf = tempfile.NamedTemporaryFile(mode='w')
        wf.write('    # Comentario\n')
        wf.write('    ! Comentario\n')
        wf.write('\n')
        wf.write('inputvar1 1 # integer\n')
        wf.write('inputvar2 1.2 # float\n')
        wf.write('inputvar3 3*4 # list of integer\n')
        wf.write('inputvar4 3*4.5 # list of float\n')
        wf.write('inputvar5 3*4.5e6 # list of float\n')
        wf.write('inputvar6 3*4.5d7 # list of float\n')
        wf.write('inputvar7 3*4.5E6 # list of float\n')
        wf.write('inputvar8 3*4.5D7 # list of float\n')
        wf.write('inputvar9 *1\n')
        wf.write('inputvar10 sqrt(2)\n')
        wf.write('inputvar11 6*sqrt(3)\n')
        wf.flush()

        inp = parser(wf.name)
        wf.close()

        assert len(inp.keys()) == 11
        assert inp['inputvar1'] == array([1])
        assert inp['inputvar2'] == array([1.2])
        assert all(inp['inputvar3'] == 4 * ones(3))
        assert all(inp['inputvar4'] == 4.5 * ones(3))
        assert all(inp['inputvar5'] == 4.5e6 * ones(3))
        assert all(inp['inputvar6'] == 4.5e7 * ones(3))
        assert all(inp['inputvar7'] == 4.5e6 * ones(3))
        assert all(inp['inputvar8'] == 4.5e7 * ones(3))
        assert inp['inputvar9'] == '*1'
        assert inp['inputvar10'] == sqrt(2)
        assert all(inp['inputvar11'] == sqrt(3) * ones(6))


def test_abinit_utils():
    """
    Tests (pychemia.code.abinit) [utils]                         :
    """
    if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
        from pychemia.code.abinit import xyz2input, netcdf2dict, psp_name

        filename = "pychemia/test/data/abinit_01/abinit-o_OUT.nc"
        print(len(netcdf2dict(filename)))
        assert len(netcdf2dict(filename)) == 310
        assert psp_name(1, 'LDA', 'FHI') == '01-H.LDA.fhi'
        filename = "pychemia/test/data/abinit_01/abinit_DS11.xyz"
        assert xyz2input(filename).variables['natom'] == 2


def test_abinit_abifiles():
    """
    Tests (pychemia.code.abinit) [abifiles]                      :
    """
    if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
        from pychemia.code.abinit import AbiFiles

        filename = "pychemia/test/data/abinit_01/abinit.files"
        abf = AbiFiles(filename)
        assert abf.filename == "abinit.files"
        assert abf.get_input_filename() == 'pychemia/test/data/abinit_01/abinit.in'


def test_abinit_input():
    """
    Tests (pychemia.code.abinit) [input]                         :
    """
    if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
        from pychemia.code.abinit import AbiFiles, InputVariables

        filename = "pychemia/test/data/abinit_01/abinit.files"
        abf = AbiFiles(filename)
        inp = InputVariables(abf)
        print(inp)
        print(len(inp))
        assert len(inp) == 31
        assert inp.get_value('ecut') == 10
        assert len(inp.get_dtsets_keys()) == 12
        assert inp.get_value('ntime', 41) == 10
        assert inp.get_value('acell', 41)[0] == 14


def test_abinit():
    """
    Tests (pychemia.code.abinit) [general]                       :
    """
    import pychemia
    if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
        import pychemia.code.abinit

        af = pychemia.code.abinit.AbiFiles(basedir='pychemia/test/data/abinit_03')
        iv = pychemia.code.abinit.InputVariables('pychemia/test/data/abinit_03/rnpg.in')
        af.set_input(iv)
        af.set_psps('LDA', 'FHI')
        af.create()
        iv.write(af.get_input_filename())
        assert len(open(af.get_input_filename()).readlines()) == 71
