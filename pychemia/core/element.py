from collections.abc import Mapping
from pychemia.utils.periodic import *

madelung_exceptions = {'Cr': ['Ar', '4s1', '3d5'],
                       'Cu': ['Ar', '4s1', '3d10'],
                       'Nb': ['Kr', '5s1', '4d4'],
                       'Mo': ['Kr', '5s1', '4d5'],
                       'Ru': ['Kr', '5s1', '4d7'],
                       'Rh': ['Kr', '5s1', '4d8'],
                       'Pd': ['Kr', '4d10'],
                       'Ag': ['Kr', '5s1', '4d10'],
                       'La': ['Xe', '6s2', '5d1'],
                       'Ce': ['Xe', '6s2', '4f1', '5d1'],
                       'Gd': ['Xe', '6s2', '4f7', '5d1'],
                       'Pt': ['Xe', '6s1', '4f14', '5d9'],
                       'Au': ['Xe', '6s1', '4f14', '5d10'],
                       'Ac': ['Rn', '7s2', '6d1'],
                       'Th': ['Rn', '7s2', '6d2'],
                       'Pa': ['Rn', '7s2', '5f2', '6d1'],
                       'U': ['Rn', '7s2', '5f3', '6d1'],
                       'Np': ['Rn', '7s2', '5f4', '6d1'],
                       'Cm': ['Rn', '7s2', '5f7', '6d1'],
                       'Lr': ['Rn', '7s2', '5f14', '7p1']}


class Element:

    def __init__(self, value=None):

        if value in atomic_symbols:
            self.symbol = value
        elif value.capitalize() in atomic_symbols:
            self.symbol = value.capitalize()
        else:
            raise ValueError('Symbol %s does not appear on the periodic table' % value)

    @property
    def name(self):
        return atomic_name(self.symbol)

    @property
    def atomic_number(self):
        return atomic_number(self.symbol)

    @property
    def group(self):
        return group(self.symbol)

    @property
    def period(self):
        return period(self.symbol)

    @property
    def block(self):
        return block(self.symbol)

    @property
    def valence(self):
        return valence(self.symbol)

    @property
    def valence_nominal(self):
        return valence_nominal(self.symbol)

    @property
    def mass(self):
        return mass(self.symbol)

    @property
    def covalent_radius(self):
        return covalent_radius(self.symbol)

    @property
    def electronegativity(self):
        return electronegativity(self.symbol)

    @property
    def crystal_structure(self):
        return crystal_structure(self.symbol)

    @property
    def phase(self):
        return phase(self.symbol)

    @property
    def boiling_point(self):
        return boiling_point(self.symbol)

    @property
    def melting_point(self):
        return melting_point(self.symbol)

    @property
    def oxidation_states(self):
        return oxidation_state(self.symbol)

    @property
    def oxidation_states_common(self):
        return oxidation_state(self.symbol, common=True)

    def __str__(self):

        ret = """
Symbol: %s
Name  : %s       

Z     : %d
Group : %d
Period: %d
Block : %s

Valence          : %f
Valence (Nominal): %f

Mass              : %f
Covalent Radius   : %f
Electronegativity : %f

Crystal Structure : %s
Phase             : %s
Boiling Point     : %f
Melting Point     : %f
        """ % (self.symbol, self.name, self.atomic_number, self.group, self.period, self.block,
               self.valence, self.valence_nominal,
               self.mass, self.covalent_radius, self.electronegativity, self.crystal_structure, self.phase,
               self.boiling_point, self.melting_point)
        return ret

    def previous_inert(self):

        inerts = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'Og']

        if self.period == 1:
            return None
        else:
            # In case the element is already a noble gas the previous one look one period above in Periodic Table
            return inerts[self.period - 2]

    @property
    def madelung_filling(self):

        order = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d',
                 '5p', '6s', '4f', '5d', '6p', '7s', '5f', '6d', '7p', '8s',
                 '5g', '6f', '7d', '8p', '9s']

        # We start with the total number of electron and get those associated to orbitals following the order
        capacities = {}
        # l quantum number
        for lqn in range(4):
            label = Element.orbital_label_from_number(lqn)
            nele = Element.max_electrons_subshell(label)
            capacities[label] = nele

        max_electrons = {}
        for ishell in order:
            label = ishell[-1]
            maxele = Element.max_electrons_subshell(label)
            max_electrons[ishell] = maxele

        ret = []
        if self.previous_inert() is not None:
            inert = self.__class__(self.previous_inert())
            ret.append(inert.symbol)
            inert_remain = inert.atomic_number
            # Consume the shells up to the previous inert atom
            numele = self.atomic_number - inert.atomic_number
        else:
            numele = self.atomic_number
            inert_remain = 0

        for i in order:
            if inert_remain >= max_electrons[i]:
                inert_remain -= max_electrons[i]
            elif inert_remain == 0:
                if numele >= max_electrons[i]:
                    numele -= max_electrons[i]
                    ret.append(i + str(max_electrons[i]))
                elif numele == 0:
                    break
                elif numele < max_electrons[i]:
                    ret.append(i + str(numele))
                    break
        return ret

    @staticmethod
    def azimuthal_quantum_number(label):
        aqn = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}
        if label not in ['s', 'p', 'd', 'f', 'g']:
            raise ValueError('Not such label for an orbital: %s' % label)
        return aqn[label]

    # lqn = angular momemtum
    @staticmethod
    def orbital_label_from_number(lqn):
        orbitals = ['s', 'p', 'd', 'f', 'g']
        if lqn not in range(4):
            raise ValueError('Not such azimuthal quantum number: %s' % lqn)
        return orbitals[lqn]

    @staticmethod
    def max_electrons_subshell(subshell):

        if subshell in ['s', 'p', 'd', 'f', 'g']:
            ll = Element.azimuthal_quantum_number(subshell)
        elif subshell in [0, 1, 2, 3, 4]:
            ll = subshell
        else:
            raise ValueError('Not a valid subshell: %s' % subshell)
        return 2 * (2 * ll + 1)

    @property
    def electronic_configuration(self):
        """
        Return the known electronic configuration including exceptions to Madelung's rule
        Based on:

        https://en.wikipedia.org/wiki/Electron_configuration#Atoms:_Aufbau_principle_and_Madelung_rule

        :return:
        """
        if self.symbol in madelung_exceptions:
            return madelung_exceptions[self.symbol]
        else:
            return self.madelung_filling

    @property
    def is_madelung_exception(self):

        return self.symbol in madelung_exceptions
