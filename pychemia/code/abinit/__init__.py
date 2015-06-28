"""
Set of classes and functions to manipulate

ABINIT '.in' input files
ABINIT '.files' files
ABINIT '_OUT.nc' output files

"""
try:
    import scipy


    from _abifiles import AbiFiles
    from _input import InputVariables, netcdf2dict
    from _utils import get_all_psps, psp_name, xyz2input, plot_simple, abihelp
    from _parser import parser
    from _htmlparser import MyHTMLParser
    from _tasks import RelaxPopulation
    from _abinit import AbinitJob

except ImportError:
    print 'SCIPY not present'

# __all__ = filter(lambda s: not s.startswith('_'), dir())
