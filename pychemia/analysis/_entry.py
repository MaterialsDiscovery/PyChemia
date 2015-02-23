__author__ = 'Guillermo Avendano-Franco'

from pychemia.db import StructureEntry


class EntryAnalysis():
    """
    This class provides a set of methods that uses not only the structural information but also
    tags, philogenetic info and/or computed properties.
    This class is created providing a StructureEntry instead of a simple Structure object
    This class uses lazy evaluation
    """

    def __init__(self, entry):
        assert (isinstance(entry, StructureEntry))
        self.entry = entry

