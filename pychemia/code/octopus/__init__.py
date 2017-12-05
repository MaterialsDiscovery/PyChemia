"""
Set of classes and functions to manipulate

OCTOPUS 'inp' input files

"""

from .input import OctopusInput
from .output import OctopusOutput
from .run import OctopusRun
from . import analysis

# __all__ = filter(lambda s: not s.startswith('_'), dir())


