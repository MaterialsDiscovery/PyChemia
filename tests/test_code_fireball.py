import pychemia
import tempfile
import zipfile
import os


def test_fireball():
    """
    Test (pychemia.code.fireball) [Reading fireball output]     :
    """
    tmpdir = tempfile.mkdtemp()
    path = os.path.abspath('tests/data/1280.zip')
    zf = zipfile.ZipFile(path)
    oldpwd = os.getcwd()
    os.chdir(tmpdir)
    zf.extract('1280/output.log')
    fo = zf.extract('1280/output.log')
    output = pychemia.code.fireball.read_fireball_stdout(fo)

    # Testing results from reading the Grand Total energy
    assert output['grand_total'][-1] == -31202.46383232
    assert len(output['grand_total']) == 313

    # Cleaning files and directory
    os.remove(fo)
    os.rmdir(tmpdir + os.sep + '1280')
    os.rmdir(tmpdir)
    os.chdir(oldpwd)
