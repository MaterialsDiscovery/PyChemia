from numpy import array, argsort
from pychemia.utils.periodic import atomic_symbols, electronegativity
from fractions import gcd as _gcd
from pychemia.utils.computing import unicode2string

__author__ = 'Guillermo Avendano-Franco'


def formula_parser(value):
    ret = {}
    jump = False
    for i in range(len(value)):
        if jump > 0:
            jump -= 1
        elif value[i].isupper():
            if i+1 < len(value) and value[i+1].islower():
                if i+2 < len(value) and value[i+2].islower():
                    specie = value[i:i+3]
                    jump = 2
                else:
                    specie = value[i:i+2]
                    jump = 1
            else:
                specie = value[i]
                jump = 0
            j = 1
            number = ''
            while True:
                if i+jump+j < len(value) and value[i+jump+j].isdigit():
                    number += value[i+jump+j]
                    j += 1
                else:
                    break
            if number == '':
                ret[specie] = 1
            else:
                ret[specie] = int(number)
    return ret


class Composition():

    def __init__(self, value):
        """
        Creates a new composition, internally it is a dictionary
        where each specie is the key and the value is an integer
        with the number of atoms of that specie

        :param value: (str, dict)
        """

        if isinstance(value, str):
            self._composition = formula_parser(value)
        elif isinstance(value, dict):
            self.set_composition(value)

    def set_composition(self, value):
        for i in value:
            assert(i in atomic_symbols)
            assert(isinstance(value[i], int))
        self._composition = value.copy

    @property
    def composition(self):
        """
        Return the composition dictionary

        :return: dict
        """
        return self._composition

    @property
    def species(self):
        """
        Return the list of species

        :return: list
        """
        return self._composition.keys()

    @property
    def values(self):
        """
        Return the number of atoms of each specie

        :return: list
        """
        return [self._composition[x] for x in self._composition]

    @property
    def gcd(self):
        """
        Return the Greatest common denominator for the composition

        :return: int
        """
        return reduce(_gcd, self.values)

    @property
    def natom(self):
        """
        Return the number of atoms in the composition

        :return: int
        """
        return sum(self.values)

    @property
    def formula(self):
        """
        Return the formula with atoms sorted alphabetically

        :return: str
        """
        return self.sorted_formula(sortby='alpha', reduced=True)

    def sorted_formula(self, sortby='alpha', reduced=True):
        """
        Return the chemical formula. It could be sorted
        alphabetically sortby='alpha' by electronegativity
        sortby='electroneg' or by Hill System sortby='Hill'

        :param sortby: (str) 'alpha' : Alphabetically
                             'electroneg' : Electronegativity
                             'hill' : Hill System

        :param reduced:
        """
        if reduced and self.gcd > 1:
            comp = Composition(self.composition)
            for i in comp.composition:
                comp._composition[i] /= self.gcd
        else:
            comp = self
        ret = ''
        if sortby == 'electroneg':
            electroneg = electronegativity(comp.species)
            sortedspecies = array(comp.species)[argsort(electroneg)]
        elif sortby == "hill":  # FIXME: Hill system exceptions not implemented
            sortedspecies = []
            presortedspecies = sorted(comp.species)
            if 'C' in presortedspecies:
                sortedspecies.append('C')
                presortedspecies.pop(presortedspecies.index('C'))
            if 'H' in presortedspecies:
                sortedspecies.append('H')
                presortedspecies.pop(presortedspecies.index('H'))
            sortedspecies += presortedspecies
        else:
            sortedspecies = sorted(comp.species)
        for specie in sortedspecies:
            ret += specie
            if comp.composition[specie] > 1:
                ret += str(comp.composition[specie])
        return ret