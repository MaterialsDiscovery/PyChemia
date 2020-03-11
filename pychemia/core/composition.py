"""
Chemical composition is just the description of the amount of atoms of each specie. In the case of clusters or
molecules, ie a finite structure, it represents the complete set of atoms. For periodic structures it represents
the species present on a cell.
"""

import re
from numpy import array, argsort
from math import gcd as _gcd
from math import pi
from pychemia.utils.periodic import atomic_symbols, electronegativity, atomic_number, covalent_radius
from pychemia.utils.computing import deep_unicode
from functools import reduce
from collections.abc import Mapping


class Composition(Mapping):
    """
    A Composition is basically a mapping between a number of species and a integer indicating how many atoms of that
    specie are present in the structure.
    A composition object do not contain geometrical information or bonding.

    The main purpose of this class is to be able to parse formulas into compositions and return string formulas sorted
    in various ways.
    """

    def __init__(self, value=None):
        """
        Creates a new composition,  currently only absolute formulas are supported.

        :param value: (str, dict) The input argument could be a string with a chemical formula or the actual dictionary
        of species and values. The order of species is not guaranteed to be preserved. A iterable of atomic symbols
        is also accepted to build a composition object.

        :rtype: Composition

        >>> cp = Composition({'Ba': 2, 'Cu': 3, 'O': 7, 'Y': 1})
        >>> cp.formula
        'Ba2Cu3O7Y'
        >>> cp = Composition('Ba2Cu3O7Y')
        >>> cp2 = Composition(cp)
        >>> len(cp2)
        4
        >>> cp.nspecies
        4
        >>> cp = Composition(['O', 'H', 'O'])
        >>> len(cp)
        2
        >>> cp['O']
        2

        """
        # The internal dictionary where atom species and numbers of atoms of each specie are stored.
        self._composition = {}
        # Convert strings and dictionaries into unicode
        if value is not None:
            value = deep_unicode(value)
        # Case 1: The input is a formula
        if isinstance(value, str):
            self._set_composition(self.formula_parser(value))
        # Case 2: The input is a dictionary
        elif isinstance(value, dict):
            self._set_composition(value)
        # Case 3: The input is another composition object
        elif isinstance(value, Composition):
            self._set_composition(value.composition)
        # Case 4: The input is an iterable of atomic symbols
        elif hasattr(value, "__len__"):
            dvalue = {}
            for i in value:
                if i in dvalue:
                    dvalue[i] += 1
                else:
                    dvalue[i] = 1
            self._set_composition(dvalue)
        else:
            self._composition = {}

    def __len__(self):
        return len(self._composition)

    def __getitem__(self, specie):
        """
        Returns the number of atoms of a given specie

        :param specie: Atomic Symbol for which the value will be returned
        :return: number of atoms of the given specie
        :rtype: int

        >>> comp = Composition('H2')
        >>> comp['H']
        2
        >>> comp['He']
        0
        """
        if specie in self._composition:
            return self._composition[specie]
        else:
            return 0

    def __repr__(self):
        """
        Evaluable representation of Composition object

        :return: Text representation that can be evaluated
        :rtype: str

        >>> cp1 = Composition('H2O')
        >>> cp2 = eval(repr(cp1))
        >>> cp2 == cp1
        True
        """
        return 'Composition(' + str(self.composition) + ')'

    def __str__(self):
        """

        :return: String representation of the composition

        >>> cp = Composition('YBa2Cu3O7')
        >>> 'Cu' in str(cp)
        True
        """
        ret = ''
        for i in self.species:
            ret += " %3s: %4d  " % (i, self.composition[i])
        return ret

    def __iter__(self):
        return iter(self.composition)

    def __contains__(self, specie):
        """True if 'specie' is present in composition

        :return: True if specie is present
        :param  specie: atomic specie
        :rtype: bool

        >>> cp = Composition('H2O')
        >>> 'He' in cp
        False
        """
        return specie in self._composition

    def _set_composition(self, value):
        """
        Checks the values of a dictionary before setting the actual composition

        :param value: (dict)
        :rtype: None
        """
        for i in value:
            assert (i in atomic_symbols)
            assert (isinstance(value[i], int))
        self._composition = value.copy()

    @property
    def composition(self):
        """Dictionary with composition

        :return: The composition dictionary
        :rtype: dict

        >>> import pprint
        >>> cp = Composition('H2O')
        >>> pprint.pprint(cp.composition)
        {'H': 2, 'O': 1}

        """
        return self._composition

    def covalent_volume(self, packing='cubes'):
        """
        :param packing: The kind of packing could be 'cubes' or 'spheres'
        :type packing: str
        :return: The volume occupied by a given formula assuming a 'cubes' packing or 'spheres' packing
        :rtype: (float)

        >>> cp = Composition('C5H10')
        >>> cp.covalent_volume()
        19.942320000000002
        >>> cp.covalent_volume(packing='spheres')
        10.441774334589468
        """
        if packing == 'cubes':
            factor = 8
        elif packing == 'spheres':
            factor = 4 * pi / 3.0
        else:
            raise ValueError('Non-valid packing: "%s"' % packing)

        # find volume of unit cell by adding cubes
        volume = 0.0
        for specie in self:
            number_atoms_specie = self.composition[specie]
            # Pack each atom in a cube (2*r)^3
            volume += factor * number_atoms_specie * covalent_radius(specie) ** 3
        return volume

    @property
    def formula(self):
        """Chemical formula

        :return: The chemical formula with atoms sorted alphabetically
        :rtype: str

        >>> cp = Composition('NaCl')
        >>> cp.formula
        'ClNa'

        """
        return self.sorted_formula(sortby='alpha', reduced=True)

    @staticmethod
    def formula_parser(value):
        """Return a dictionary from a chemical formula

        :return: Convert an string representing a chemical formula into a dictionary with the species as keys
                 and values as the number of atoms of that specie
        :param value: (str) Chemical formula
        :rtype: dict

        >>> import pprint
        >>> Composition.formula_parser('Au20')
        {'Au': 20}
        >>> ret = Composition.formula_parser('UutUupUusUuo')
        >>> pprint.pprint(ret)
        {'Uuo': 1, 'Uup': 1, 'Uus': 1, 'Uut': 1}
        """
        ret = {}
        jump = False
        for i in range(len(value)):
            if jump > 0:  # This char belongs to the current atom, move on
                jump -= 1
            elif value[i].isupper():  # Atom Name starts with Uppercase
                if i + 1 < len(value) and value[i + 1].islower():  # Atom name has more than 1 char
                    if i + 2 < len(value) and value[i + 2].islower():  # Atom name has more than 2 chars
                        specie = value[i:i + 3]
                        jump = 2
                    else:
                        specie = value[i:i + 2]
                        jump = 1
                else:
                    specie = value[i]
                    jump = 0
                j = 1
                number = ''
                while True:
                    if i + jump + j < len(value) and value[i + jump + j].isdigit():
                        number += value[i + jump + j]
                        j += 1
                    else:
                        break
                if number == '':
                    ret[specie] = 1
                else:
                    ret[specie] = int(number)
        return ret

    @staticmethod
    def formula_to_list(formula, nunits=1):
        """
        Reads a formula and returns a list of atomic symbols consistent with the formula and the number of
        formulas given by nunits

        :param formula: (str) Chemical formula as string
        :param nunits: (int) Number of formulas to apply
        :return: list of atomic symbols
        :rtype: list

        >>> Composition.formula_to_list('NaCl')
        ['Na', 'Cl']
        >>> flist = Composition.formula_to_list(u'Uut2Uup3Uus4Uuo5')
        >>> len(flist)
        14
        >>> flist = Composition.formula_to_list('Uut2Uup3Uus4Uuo5', nunits=2)
        >>> len(flist)
        28
        """
        # decompose composition
        a = re.findall(r"[A-Z][a-z0-9]*", formula)
        composition = []
        for i in a:
            m = re.match(r"([A-Za-z]+)([0-9]*)", i)
            if m.group(2) == "":
                n = int(1)
            else:
                n = int(m.group(2))

            for j in range(n * nunits):
                composition.append(m.group(1))

        return composition

    @property
    def gcd(self):
        """ Number of minimal formulas on a given composition.

        :return: The number of formulas that can be extracted from a composition ie, the greatest common denominator
                 for the composition.
        :rtype: int

        >>> cp = Composition('NaCl')
        >>> cp.gcd
        1
        >>> cp = Composition('Na2Cl2')
        >>> cp.gcd
        2
        >>> cp = Composition()
        >>> cp.gcd is None
        True
        """
        if self.natom > 0:
            return reduce(_gcd, self.values)
        else:
            return None

    @staticmethod
    def get_species_from_hex(arg):
        """List of species encoded for hex string produced by species_hex

        :return: Return a set of species from the encoded species created by the output of "species_hex" method.
        :param arg: str String with hexadecimal representation of list of species.

        >>> Composition.get_species_from_hex('0x38271d08')
        [8, 29, 39, 56]
        """
        num = int(arg, 16)
        ret = []
        while num > 0:
            ret.append(num % 256)
            num = (num-ret[-1])//256
        return ret

    @property
    def natom(self):
        """
        :return: The number of atoms in the composition
        :rtype: int

        >>> cp = Composition('H2O')
        >>> cp.natom
        3
        """
        return sum(self.values)

    @property
    def nspecies(self):
        """
        :return: Number of species in the composition
        :rtype: int

        >>> cp = Composition('H2O')
        >>> cp.nspecies
        2
        """
        return len(self.species)

    @property
    def symbols(self):
        """List of species on the composition

        :return: A list of atomic symbols
        :rtype: list

        >>> cp = Composition('H2O')
        >>> cp.symbols
        ['H', 'H', 'O']
        """
        ret = []
        for specie in self:
            number_atoms_specie = self.composition[specie]
            for i in range(number_atoms_specie):
                ret.append(specie)
        return sorted(deep_unicode(ret))

    @property
    def species(self):
        """List of species on the composition

        :return: The list of species, no particular order but atoms of the same specie are contiguous.
        :rtype: list

        >>> cp = Composition('H2O')
        >>> sorted(cp.species)
        ['H', 'O']
        """
        return [deep_unicode(x) for x in self._composition]

    def sorted_formula(self, sortby='alpha', reduced=True):
        """
        :return: The chemical formula. It could be sorted  alphabetically using sortby='alpha', by electronegativity
                 using sortby='electronegativity' or using Hill System with sortby='Hill'
                 Just the first 3 letters are unambiguous and case is not taken in account so you can use 'alp', 'hil'
                 or 'ele'
        :param sortby: (str) 'alpha' : Alphabetically
                             'electronegativity' : Electronegativity
                             'hill' : Hill System
        :param reduced: (bool) If the formula should be normalized
        :rtype: str

        .. notes: Hill exceptions have not being implemented yet

        >>> cp = Composition('YBa2Cu3O7')
        >>> cp.sorted_formula()
        'Ba2Cu3O7Y'
        >>> cp.sorted_formula(sortby='hill')
        'Ba2Cu3O7Y'
        >>> cp.sorted_formula(sortby='electroneg')
        'Ba2YCu3O7'
        >>> cp = Composition('H10C5')
        >>> cp.sorted_formula(sortby='hill', reduced=True)
        'CH2'
        >>> cp = Composition('IBr')
        >>> cp.sorted_formula(sortby='hill', reduced=False)
        'BrI'
        >>> cp = Composition('Cl4C')
        >>> cp.sorted_formula(sortby='hill', reduced=False)
        'CCl4'
        >>> cp = Composition('IH3C')
        >>> cp.sorted_formula(sortby='hill', reduced=False)
        'CH3I'
        >>> cp = Composition('BrH5C2')
        >>> cp.sorted_formula(sortby='hill', reduced=False)
        'C2H5Br'
        >>> cp = Composition('S04H2')
        >>> cp.sorted_formula(sortby='hill', reduced=False)
        'H2S4'
        >>> cp = Composition('SO4H2')
        >>> cp.sorted_formula(sortby='hill', reduced=False)
        'H2O4S'
        """
        if reduced and self.gcd > 1:
            comp = Composition(self.composition)
            for i in comp.composition:
                comp._composition[i] //= self.gcd
        else:
            comp = self
        if sortby.lower()[:3] == 'ele':
            electroneg = list(electronegativity(comp.species))
            # Not longer needed as electronegativy will return 0 for 'None' values
            # for i in range(len(electroneg)):
            #    if electroneg[i] is None:
            #        electroneg[i] = -1
            sortedspecies = array(comp.species)[argsort(electroneg)]
        elif sortby.lower()[:3] == "hil":  # FIXME: Hill system exceptions not implemented
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
        ret = u''
        for specie in sortedspecies:
            ret += '%s' % specie
            if comp.composition[specie] > 1:
                ret += "%d" % comp.composition[specie]
        return deep_unicode(ret)

    def species_encoded(self, base):
        """Encode the list of species with a number

        :return: Encodes the species as a number.
        :param base: Integer used as base for encoding.
        :rtype: int

        >>> cp = Composition('H2O')
        >>> cp.species_encoded(100)
        801
        """
        ret = 0
        i = 0
        for atom_number in sorted(atomic_number(self.species)):
            ret += atom_number * (base ** i)
            i += 1
        return ret

    def species_hex(self):
        """Encoding in hexadecimal with 2 bytes per specie (base 256)

        :return: Encodes the species into a hexadecimal representation where each specie is stored on a 2-Byte slot
                 ordered by atomic number.
                 The output produces a unique encoding where each 2 character from the hexadecimal will encode a single
                 species and the species are ordered by atomic number making the codification unique.
        :rtype: str

        >>> cp = Composition('YBa2Cu3O7')
        >>> cp.species_hex()
        '0x38271d08'
        """
        enc = self.species_encoded(256)
        return hex(enc)

    @property
    def values(self):
        """
        :return: The number of atoms of each specie
        :rtype: list

        >>> cp = Composition('YBa2Cu3O7')
        >>> sorted(cp.values)
        [1, 2, 3, 7]
        """
        return [self._composition[x] for x in self._composition]


