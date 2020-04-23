
import numpy as _np

# 118 Elements to date
atomic_names = ['', 'Hydrogen', 'Helium',
                'Lithium', 'Beryllium', 'Boron', 'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon',
                'Sodium', 'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur', 'Chlorine', 'Argon',
                'Potassium', 'Calcium', 'Scandium', 'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron', 'Cobalt',
                'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium', 'Arsenic', 'Selenium', 'Bromine', 'Krypton',
                'Rubidium', 'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum', 'Technetium', 'Ruthenium',
                'Rhodium', 'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium', 'Iodine',
                'Xenon', 'Caesium', 'Barium', 'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
                'Samarium', 'Europium', 'Gadolinium', 'Terbium', 'Dysprosium', 'Holmium', 'Erbium', 'Thulium',
                'Ytterbium', 'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium', 'Osmium', 'Iridium', 'Platinum',
                'Gold', 'Mercury', 'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine', 'Radon', 'Francium',
                'Radium', 'Actinium', 'Thorium', 'Protactinium', 'Uranium', 'Neptunium', 'Plutonium', 'Americium',
                'Curium', 'Berkelium', 'Californium', 'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
                'Lawrencium', 'Rutherfordium', 'Dubnium', 'Seaborgium', 'Bohrium', 'Hassium', 'Meitnerium',
                'Darmstadtium', 'Roentgenium', 'Copernicium', 'Nihonium', 'Flerovium', 'moscovium', 'Livermorium',
                'Tennessine', 'Oganesson']

electronegativities = [None, 2.2, None, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98,
                       None, 0.93, 1.31, 1.61, 1.9, 2.19, 2.58, 3.16, None, 0.82,
                       1.0, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.9,
                       1.65, 1.81, 2.01, 2.18, 2.55, 2.96, 3.0, 0.82, 0.95, 1.22,
                       1.33, 1.6, 2.16, 1.9, 2.2, 2.28, 2.2, 1.93, 1.69, 1.78,
                       1.96, 2.05, 2.1, 2.66, 2.6, 0.79, 0.89, 1.1, 1.12, 1.13,
                       1.14, 1.13, 1.17, 1.2, 1.2, 1.2, 1.22, 1.23, 1.24, 1.25,
                       1.1, 1.27, 1.3, 1.5, 2.36, 1.9, 2.2, 2.2, 2.28, 2.54,
                       2.0, 1.62, 1.87, 2.02, 2.0, 2.2, 2.2, 0.7, 0.9, 1.1,
                       1.3, 1.5, 1.38, 1.36, 1.28, 1.3, 1.3, 1.3, 1.3, 1.3,
                       1.3, 1.3, 1.3, None, None, None, None, None, None, None,
                       None, None, None, None, None, None, None, None, None]

# Melting points in Kelvin (Extracted from Wikipedia)
melting_points = [None, 13.99, 0.95, 453.65, 1560, 2349, None, 63.15, 54.36, 53.48,
                  24.56, 370.944, 923, 933.47, 1687, None, 388.36, 171.6, 83.81, 336.7,
                  1115, 1814, 1941, 2183, 2180, 1519, 1811, 1768, 1728, 1357.77,
                  692.68, 302.9146, 1211.4, None, 494, 265.8, 115.78, 312.45, 1050, 1799,
                  2128, 2750, 2896, 2430, 2607, 2237, 1828.05, 1234.93, 594.22, 429.7485,
                  505.08, 903.78, 722.66, 386.85, 161.4, 301.7, 1000, 1193, 1068, 1208,
                  1297, 1315, 1345, 1099, 1585, 1629, 1680, 1734, 1802, 1818,
                  1097, 1925, 2506, 3290, 3695, 3459, 3306, 2739, 2041.4, 1337.33,
                  234.321, 577, 600.61, 544.7, 527, 575, 202.0, 300, 973, 1323,
                  2115, 1841, 1405.3, 912, 912.5, 1449, 1613, 1259, 1173, 1133,
                  1125, 1100, 1100, None, 2400, None, None, None, None, None,
                  None, None, None, 700, 340, 670, None, None, None]

# Boiling points in Kelvin (Extracted from Wikipedia)
boiling_points = [None, 20.271, 4.222, 1603, 3243, 4200, None, 77.355, 90.188, 85.03,
                  27.104, 1156.09, 1363, 2743, 3538, 553.7, 717.8, 239.11, 87.302, 1032,
                  1757, 3109, 3560, 3680, 2944, 2334, 3134, 3200, 3003, 2835,
                  1180, 2673, 3106, None, 958, 332.0, 119.93, 961, 1650, 3203,
                  4650, 5017, 4912, 4538, 4423, 3968, 3236, 2435, 1040, 2345,
                  2875, 1908, 1261, 457.4, 165.051, 944, 2118, 3737, 3716, 3403,
                  3347, 3273, 2173, 1802, 3273, 3396, 2840, 2873, 3141, 2003,
                  1469, 3675, 4876, 5731, 6203, 5869, 5285, 4403, 4098, 3243,
                  629.88, 1746, 2022, 1837, 1235, 610, 211.5, 950, 2010, 3471,
                  5061, 4300, 4404, 4447, 3505, 2880, 3383, 2900, 1743, 1269,
                  None, None, None, None, 5800, None, None, None, None, None,
                  None, None, None, 1430, 420, 1400, 1035, 823, 350]

crystal_structures = [None, 'hexagonal', 'hcp',
                      'bcc', 'hcp', 'rhombohedral', 'diamond', 'hexagonal', 'cubic', 'monoclinic base-centered',
                      'fcc', 'bcc', 'hcp', 'fcc', 'diamond cubic', 'simple triclinic', 'orthorhombic', 'orthorhombic',
                      'fcc', 'bcc', 'fcc', 'hcp', 'hcp', 'bcc', 'bcc', 'bcc', 'bcc', 'hcp', 'fcc', 'fcc',
                      'hcp', 'orthorhombic', 'diamond cubic', 'simple trigonal', 'hexagonal', 'orthorhombic',
                      'cubic face-centered', 'bcc', 'fcc', 'hcp', 'hcp', 'cubic body-centered', 'bcc', 'hcp', 'hcp',
                      'fcc', 'fcc', 'fcc', 'hcp', 'tetragonal', 'tetragonal', 'simple trigonal', 'hexagonal',
                      'orthorhombic', 'fcc', 'bcc', 'bcc', 'hexagonal', 'fcc', 'hexagonal', 'hexagonal', 'hexagonal',
                      'rhombohedral', 'bcc', 'hcp', 'hcp', 'hcp', 'hcp', 'hcp', 'hcp', 'fcc', 'hcp', None, 'bcc',
                      'bcc', 'hcp', 'hcp', 'fcc', 'fcc', 'fcc', 'rhombohedral', 'hcp', 'fcc', 'rhombohedral', 'cubic',
                      None, None, 'bcc', 'bcc', 'fcc', 'fcc', 'tetragonal', 'orthorhombic', 'orthorhombic',
                      'monoclinic', 'hexagonal', 'hcp', 'hcp', 'simple hexagonal', 'fcc', None, None, None, 'hcp',
                      'hcp', 'bcc', 'bcc', 'hcp', 'hcp', 'fcc', 'bcc', 'bcc', 'hcp', None, None, None, None, None, None]

phases = [None, 'gas', 'gas', 'solid', 'solid', 'solid', 'solid', 'gas', 'gas', 'gas',
          'gas', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'gas', 'gas', 'solid',
          'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid',
          'solid', 'solid', 'solid', 'solid', 'solid', 'liquid', 'gas', 'solid', 'solid', 'solid',
          'solid', 'solid', 'solid', 'solid', None, 'solid', 'solid', 'solid', 'solid', 'solid',
          'solid', 'solid', 'solid', 'solid', 'gas', 'solid', 'solid', 'solid', 'solid', 'solid',
          'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid',
          'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid',
          'liquid', 'solid', 'solid', 'solid', 'solid', 'solid', 'gas', 'solid', 'solid', 'solid',
          'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid',
          'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid',
          'solid', 'solid', 'unknown', 'solid', 'solid', 'solid', 'solid', 'solid', 'solid']

# Van der Waals radii in Angstroms
van_der_waals_radii = [None, 1.20, 1.40, 1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.35,
                       1.54, 2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75,
                       2.31, 2.11, None, None, None, None, None, None, 1.63, 1.40,
                       1.39, 1.87, 2.11, 1.85, 1.90, 1.85, 2.02, 3.03, 2.49, None,
                       None, None, None, None, None, None, 1.63, 1.72, 1.58, 1.93,
                       2.17, 2.06, 2.06, 1.98, 2.16, 3.43, 2.68, None, None, None,
                       None, None, None, None, None, None, None, None, None, None,
                       None, None, None, None, None, None, None, None, 1.75, 1.66,
                       1.55, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, None,
                       None, None, 1.86, None, None, None, None, None, None, None,
                       None, None, None, None, None, None, None, None, None, None,
                       None, None, None, None, None, None, None, None, None]

# 118 Elements to date
atomic_symbols = ['',
                  'H', 'He',
                  'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                  'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                  'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
                  'Kr',
                  'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
                  'Xe',
                  'Cs', 'Ba',
                  'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                  'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
                  'Fr', 'Ra',
                  'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
                  'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

valences_nominal = [0,
                    1, 0,
                    1, 2, 3, 4, 5, 2, 1, 0,
                    1, 2, 3, 4, 5, 6, 7, 0,
                    1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 4, 2, 3, 4, 5, 6, 7, 0,
                    1, 2, 3, 4, 5, 6, 7, 8, 6, 4, 4, 2, 3, 4, 5, 6, 7, 0,
                    1, 2,
                    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                    4, 5, 6, 7, 8, 8, 6, 5, 4, 3, 4, 5, 6, 7, 0,
                    1, 2,
                    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                    4, 5, 6, 7, 8,
                    None, None, None, None, None, None, None, None, None, None]

valences = [0,
            1, 0,
            1, 2, 3, 4, 5, 2, 1, 0,
            1, 2, 3, 4, 5, 6, 7, 0,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 0,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 0,
            1, 2,
            3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 3,
            8, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 0,
            1, 2,
            3, 4, 5, 6, 7, 8, 9, 9, 11, 12, 13, 14, 15, 16, 17,
            4, 5, 6, 7, 8,
            None, None, None, None, None, None, None, None, None, None]

periods = [None, 1, 1,
           2, 2, 2, 2, 2, 2, 2, 2,
           3, 3, 3, 3, 3, 3, 3, 3,
           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
           5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
           6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
           7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]

groups = [None,
          1, 18,
          1, 2, 13, 14, 15, 16, 17, 18,
          1, 2, 13, 14, 15, 16, 17, 18,
          1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
          1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
          1, 2,
          -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
          4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
          1, 2,
          -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
          4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]

blocks = [None,
          's', 'p',
          's', 's', 'p', 'p', 'p', 'p', 'p', 'p',
          's', 's', 'p', 'p', 'p', 'p', 'p', 'p',
          's', 's', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p',
          's', 's', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p',
          's', 's',
          'd', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f',
          'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p',
          's', 's',
          'd', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f',
          'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p']

# covalent_radii=[0.00,
#                 0.31,0.28,1.28,0.96,0.76,0.73,0.71,0.66,0.57,0.58,
#                 1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,2.03,1.76,
#                 1.70,1.60,1.53,1.39,1.39,1.32,1.26,1.24,1.32,1.22,
#                 1.22,1.20,1.19,1.20,1.20,1.16,2.20,1.95,1.90,1.75,
#                 1.64,1.54,1.47,.146,1.42,1.39,1.45,1.44,1.42,1.39,
#                 1.39,1.38,1.39,1.40,2.44,2.15,2.07,2.04,2.03,2.01,
#                 1.99,1.98,1.98,1.96,1.94,1.92,1.92,1.89,1.90,1.87,
#                 1.87,1.75,1.70,1.62,1.51,1.44,1.41,1.36,1.36,1.32,
#                 1.45,1.46,1.48,1.40,1.50,1.50,2.60,2.21,2.15,2.06,
#                 2.00,1.96,1.90,1.87,1.80,1.69,1.66,1.68,1.65,1.67,
#                 1.73,1.76,1.61,1.57,1.49,1.43,1.41,1.34,1.29,1.28,
#                 1.21,1.22,1.36,1.43,1.62,1.75,1.65,1.57]

# Values up to 96:
#
# Beatriz Cordero; Ver\'onica G\'omez; Ana E. Platero-Prats; Marc Rev\'es; Jorge Echeverr\'\ia; Eduard Cremades;
# Flavia Barrag\'an; Santiago Alvarez (2008).
# "Covalent radii revisited". Dalton Trans. (21): 2832-2838. doi:10.1039/b801115j
# 
# Higher values (Single bond values were used):
# P. Pyykk\:o; M. Atsumi (2009). "Molecular Single-Bond Covalent Radii for Elements 1-118".
# Chemistry: A European Journal. 15: 186-197. doi:10.1002/chem.200800987.
#
# Also see:
# https://en.wikipedia.org/wiki/Covalent_radius

covalent_radii = [0.20, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66,  # H .. O
                  0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02,  # F .. Cl
                  1.06, 2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32,  # Ar .. Fe
                  1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20,  # Co .. Br
                  1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46,  # Kr .. Ru
                  1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39,  # Rh .. I
                  1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98,  # Xe .. Sm
                  1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,  # Eu .. Lu
                  1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,  # Hf .. Hg
                  1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15,  # Tl .. Ac
                  2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 1.68, 1.68,  # Th .. Cf
                  1.65, 1.67, 1.73, 1.76, 1.61, 1.57, 1.49, 1.43, 1.41,  # Es .. Bh
                  1.34, 1.29, 1.28, 1.21, 1.22, 1.36, 1.43, 1.62, 1.75,  # Hs .. Lv
                  1.65, 1.57]  # Ts .. Og

masses = [0.0, 1.00794, 4.002602, 6.941, 9.012182,   # None + H - Be
          10.811, 12.011, 14.00674, 15.9994, 18.9984032,    # B - F
          20.1797, 22.989768, 24.3050, 26.981539, 28.0855,  # Ne -Si
          30.973762, 32.066, 35.4527, 39.948, 39.0983,      # P - K
          40.078, 44.955910, 47.88, 50.9415, 51.9961,       # Ca - Cr
          54.93805, 55.847, 58.93320, 58.69, 63.546,        # Mn - Cu
          65.39, 69.723, 72.61, 74.92159, 78.96,            # Zn - Se
          79.904, 83.80, 85.4678, 87.62, 88.90585,          # Br - Y
          91.224, 92.90638, 95.94, 98.9062, 101.07,         # Zr - Ru
          102.9055, 106.42, 107.8682, 112.411, 114.82,      # Rh - In
          118.710, 121.753, 127.60, 126.90447, 131.29,      # Sn - Xe
          132.90543, 137.327, 138.9055, 140.115, 140.90765,  # Cs - Pr
          144.24, 147.91, 150.36, 151.965, 157.25,          # Nd - Gd
          158.92534, 162.50, 164.93032, 167.26, 168.93421,  # Tb - Tm
          173.04, 174.967, 178.49, 180.9479, 183.85,        # Yb - W
          186.207, 190.2, 192.22, 195.08, 196.96654,        # Re - Au
          200.59, 204.3833, 207.2, 208.98037, 209.0,        # Hg - Po
          210.0, 222.0, 223.0, 226.0254, 230.0,             # At - Ac
          232.0381, 231.0359, 238.0289, 237.0482, 242.0,    # Th - Pu
          243.0, 247.0, 247.0, 249.0, 254.0,                # Am - Es
          253.0, 256.0, 254.0, 257.0, 260.0,                # Fm - Rf
          268.0,  # Db
          269.0,  # Sg
          270.0,  # Bh
          270.0,  # Hs
          278.0,  # Mt
          281.0,  # Ds
          282.0,  # Rg
          285.0,  # Cn
          286.0,  # Nh
          290.0,  # Fl
          290.0,  # Mc
          293.0,  # Lv
          294.0,  # Ts
          294.0]  # Og

cpk_colors = [(0.0, 0.0, 0.0),   #   0 None
              (1.0, 1.0, 1.0),     #   1 H   Hydrogen
              (1.0, 0.75, 0.8),    #   2 He  Helium
              (0.7, 0.13, 0.13),   #   3 Li  Lithium
              (1.0, 0.08, 0.58),   #   4 Be  Beryllium
              (0.0, 1.0, 0.0),     #   5 B   Boron
              (0.78, 0.78, 0.78),  #   6 C   Carbon
              (0.56, 0.56, 1.0),   #   7 N   Nitrogen
              (0.94, 0.0, 0.0),    #   8 O   Oxygen
              (0.85, 0.65, 0.13),  #   9 F   Fluorine
              (1.0, 0.08, 0.58),   #  10 Ne  Neon
              (0.0, 0.0, 1.0),     #  11 Na  Sodium
              (0.13, 0.55, 0.13),  #  12 Mg  Magnesium
              (0.5, 0.5, 0.56),    #  13 Al  Aluminium
              (0.85, 0.65, 0.13),  #  14 Si  Silicon
              (1.0, 0.65, 0.0),    #  15 P   Phosphorus
              (1.0, 0.78, 0.2),    #  16 S   Sulfur
              (0.0, 1.0, 0.0),     #  17 Cl  Chlorine
              (1.0, 0.08, 0.58),   #  18 Ar  Argon
              (1.0, 0.08, 0.58),   #  19 K   Potassium
              (0.5, 0.5, 0.56),    #  20 Ca  Calcium
              (1.0, 0.08, 0.58),   #  21 Sc  Scandium
              (0.5, 0.5, 0.56),    #  22 Ti  Titanium
              (1.0, 0.08, 0.58),   #  23 V   Vanadium
              (0.5, 0.5, 0.56),    #  24 Cr  Chromium
              (0.5, 0.5, 0.56),    #  25 Mn  Manganese
              (1.0, 0.65, 0.0),    #  26 Fe  Iron
              (1.0, 0.08, 0.58),   #  27 Co  Cobalt
              (0.65, 0.16, 0.16),  #  28 Ni  Nickel
              (0.65, 0.16, 0.16),  #  29 Cu  Copper
              (0.65, 0.16, 0.16),  #  30 Zn  Zinc
              (1.0, 0.08, 0.58),   #  31 Ga  Gallium
              (1.0, 0.08, 0.58),   #  32 Ge  Germanium
              (1.0, 0.08, 0.58),   #  33 As  Arsenic
              (1.0, 0.08, 0.58),   #  34 Se  Selenium
              (0.65, 0.16, 0.16),  #  35 Br  Bromine
              (1.0, 0.08, 0.58),   #  36 Kr  Krypton
              (1.0, 0.08, 0.58),   #  37 Rb  Rubidium
              (1.0, 0.08, 0.58),   #  38 Sr  Strontium
              (1.0, 0.08, 0.58),   #  39 Y   Yttrium
              (1.0, 0.08, 0.58),   #  40 Zr  Zirconium
              (1.0, 0.08, 0.58),   #  41 Nb  Niobium
              (1.0, 0.08, 0.58),   #  42 Mo  Molybdenum
              (1.0, 0.08, 0.58),   #  43 Tc  Technetium
              (1.0, 0.08, 0.58),   #  44 Ru  Ruthenium
              (1.0, 0.08, 0.58),   #  45 Rh  Rhodium
              (1.0, 0.08, 0.58),   #  46 Pd  Palladium
              (0.5, 0.5, 0.56),    #  47 Ag  Silver
              (1.0, 0.08, 0.58),   #  48 Cd  Cadmium
              (1.0, 0.08, 0.58),   #  49 In  Indium
              (1.0, 0.08, 0.58),   #  50 Sn  Tin
              (1.0, 0.08, 0.58),   #  51 Sb  Antimony
              (1.0, 0.08, 0.58),   #  52 Te  Tellurium
              (0.63, 0.13, 0.94),  #  53 I   Iodine
              (1.0, 0.08, 0.58),   #  54 Xe  Xenon
              (1.0, 0.08, 0.58),   #  55 Cs  Caesium
              (1.0, 0.65, 0.0),    #  56 Ba  Barium
              (1.0, 0.08, 0.58),   #  57 La  Lanthanum
              (1.0, 0.08, 0.58),   #  58 Ce  Cerium
              (1.0, 0.08, 0.58),   #  59 Pr  Praseodymium
              (1.0, 0.08, 0.58),   #  60 Nd  Neodymium
              (1.0, 0.08, 0.58),   #  61 Pm  Promethium
              (1.0, 0.08, 0.58),   #  62 Sm  Samarium
              (1.0, 0.08, 0.58),   #  63 Eu  Europium
              (1.0, 0.08, 0.58),   #  64 Gd  Gadolinium
              (1.0, 0.08, 0.58),   #  65 Tb  Terbium
              (1.0, 0.08, 0.58),   #  66 Dy  Dysprosium
              (1.0, 0.08, 0.58),   #  67 Ho  Holmium
              (1.0, 0.08, 0.58),   #  68 Er  Erbium
              (1.0, 0.08, 0.58),   #  69 Tm  Thulium
              (1.0, 0.08, 0.58),   #  70 Yb  Ytterbium
              (1.0, 0.08, 0.58),   #  71 Lu  Lutetium
              (1.0, 0.08, 0.58),   #  72 Hf  Hafnium
              (1.0, 0.08, 0.58),   #  73 Ta  Tantalum
              (1.0, 0.08, 0.58),   #  74 W   Tungsten
              (1.0, 0.08, 0.58),   #  75 Re  Rhenium
              (1.0, 0.08, 0.58),   #  76 Os  Osmium
              (1.0, 0.08, 0.58),   #  77 Ir  Iridium
              (1.0, 0.08, 0.58),   #  78 Pt  Platinum
              (0.85, 0.65, 0.13),  #  79 Au  Gold
              (0.72, 0.72, 0.82),  #  80 Hg  Mercury (elemen
              (1.0, 0.08, 0.58),   #  81 Tl  Thallium
              (1.0, 0.08, 0.58),   #  82 Pb  Lead
              (1.0, 0.08, 0.58),   #  83 Bi  Bismuth
              (1.0, 0.08, 0.58),   #  84 Po  Polonium
              (1.0, 0.08, 0.58),   #  85 At  Astatine
              (1.0, 0.08, 0.58),   #  86 Rn  Radon
              (1.0, 0.08, 0.58),   #  87 Fr  Francium
              (1.0, 0.08, 0.58),   #  88 Ra  Radium
              (1.0, 0.08, 0.58),   #  89 Ac  Actinium
              (1.0, 0.08, 0.58),   #  90 Th  Thorium
              (1.0, 0.08, 0.58),   #  91 Pa  Protactinium
              (1.0, 0.08, 0.58),   #  92 U   Uranium
              (1.0, 0.08, 0.58),   #  93 Np  Neptunium
              (1.0, 0.08, 0.58),   #  94 Pu  Plutonium
              (1.0, 0.08, 0.58),   #  95 Am  Americium
              (1.0, 0.08, 0.58),   #  96 Cm  Curium
              (1.0, 0.08, 0.58),   #  97 Bk  Berkelium
              (1.0, 0.08, 0.58),   #  98 Cf  Californium
              (1.0, 0.08, 0.58),   #  99 Es  Einsteinium
              (1.0, 0.08, 0.58),   # 100 Fm  Fermium
              (1.0, 0.08, 0.58),   # 101 Md  Mendelevium
              (1.0, 0.08, 0.58),   # 102 No  Nobelium
              (1.0, 0.08, 0.58),   # 103 Lr  Lawrencium
              (1.0, 0.08, 0.58),   # 104 Rf  Rutherfordium
              (1.0, 0.08, 0.58),   # 105 Db  Dubnium
              (1.0, 0.08, 0.58),   # 106 Sg  Seaborgium
              (1.0, 0.08, 0.58),   # 107 Bh  Bohrium
              (1.0, 0.08, 0.58),   # 108 Hs  Hassium
              (1.0, 0.08, 0.58),   # 109 Mt  Meitnerium
              (1.0, 0.08, 0.58),   # 110 Ds  Darmstadtium
              (1.0, 0.08, 0.58),   # 111 Rg  Roentgenium
              (1.0, 0.08, 0.58),   # 112 Cn  Copernicium
              (1.0, 0.08, 0.58),   # 113 Nh  Nihonium
              (1.0, 0.08, 0.58),   # 114 Fl  Flerovium
              (1.0, 0.08, 0.58),   # 115 Mc  Moscovium
              (1.0, 0.08, 0.58),   # 116 Lv  Livermorium
              (1.0, 0.08, 0.58),   # 117 Ts  Tennessine
              (1.0, 0.08, 0.58)]   # 118 Og  Oganesson

# Data from:
# https://en.wikipedia.org/wiki/List_of_oxidation_states_of_the_elements

oxidation_states = [(),
                    (-1, 1),  # H
                    (),  # He
                    (1,),  # Li
                    (1, 2),  # Be
                    (-5, -1, 1, 2, 3),  # B
                    (-4, -3, -2, -1, 1, 2, 3, 4),  # C
                    (-3, -2, -1, 1, 2, 3, 4, 5),  # N
                    (-2, -1, 1, 2,),  # O
                    (-1,),  # F
                    (),  # Ne
                    (-1, 1),  # Na
                    (1, 2),  # Mg
                    (-2, -1, 1, 2, 3),  # Al
                    (-4, -3, -2, -1, 1, 2, 3, 4),  # Si
                    (-3, -2, -1, 1, 2, 3, 4, 5),  # P
                    (-2, -1, 1, 2, 3, 4, 5, 6),  # S
                    (-1, 1, 2, 3, 4, 5, 6, 7),  # Cl
                    (),  # Ar
                    (-1, 1),  # K
                    (1, 2),  # Ca
                    (1, 2, 3),  # Sc
                    (-2, -1, 1, 2, 3, 4),  # Ti
                    (-3, -1, 1, 2, 3, 4, 5),  # V
                    (-4, -2, -1, 1, 2, 3, 4, 5, 6),  # Cr
                    (-3, -2, -1, 1, 2, 3, 4, 5, 6, 7),  # Mn
                    (-4, -2, -1, 1, 2, 3, 4, 5, 6, 7),  # Fe
                    (-3, -1, 1, 2, 3, 4, 5),  # Co
                    (-2, -1, 1, 2, 3, 4),  # Ni
                    (-2, 1, 2, 3, 4),  # Cu
                    (-2, 1, 2),  # Zn
                    (-5, -4, -2, -1, 1, 2, 3),  # Ga
                    (-4, -3, -2, -1, 1, 2, 3, 4),  # Ge
                    (-3, -2, -1, 1, 2, 3, 4, 5),  # As
                    (-2, -1, 1, 2, 3, 4, 5, 6),  # Se
                    (-1, 1, 3, 4, 5, 7),  # Br
                    (2,),  # Kr
                    (-1, 1),  # Rb
                    (1, 2),  # Sr
                    (1, 2, 3),  # Y
                    (-2, 1, 2, 3, 4),  # Zr
                    (-3, -1, 1, 2, 3, 4, 5),  # Nb
                    (-4, -2, -1, 1, 2, 3, 4, 5, 6),  # Mo
                    (-3, -1, 1, 2, 3, 4, 5, 6, 7),  # Tc
                    (-4, -2, 1, 2, 3, 4, 5, 6, 7, 8),  # Ru
                    (-3, -1, 1, 2, 3, 4, 5, 6),  # Rh
                    (1, 2, 3, 4),  # Pd
                    (-2, -1, 1, 2, 3),  # Ag
                    (-2, 1, 2),  # Cd
                    (-5, -2, -1, 1, 2, 3),  # In
                    (-4, -3, -2, -1, 1, 2, 3, 4),  # Sn
                    (-3, -2, -1, 1, 2, 3, 4, 5),  # Sb
                    (-2, -1, 1, 2, 3, 4, 5, 6),  # Te
                    (-1, 1, 3, 4, 5, 6, 7),  # I
                    (2, 4, 6, 8),  # Xe
                    (-1, 1),  # Cs
                    (1, 2),  # Ba
                    (1, 2, 3),  # La
                    (2, 3, 4),  # Ce
                    (2, 3, 4, 5),  # Pr
                    (2, 3, 4),  # Nd
                    (2, 3),  # Pm
                    (2, 3),  # Sm
                    (2, 3),  # Eu
                    (1, 2, 3),  # Gd
                    (1, 2, 3, 4),  # Tb
                    (2, 3, 4),  # Dy
                    (2, 3),  # Ho
                    (2, 3),  # Er
                    (2, 3),  # Tm
                    (2, 3),  # Yb
                    (2, 3),  # Lu
                    (-2, 1, 2, 3, 4),  # Hf
                    (-3, -1, 1, 2, 3, 4, 5),  # Ta
                    (-4, -2, -1, 1, 2, 3, 4, 5, 6),  # W
                    (-3, -1, 1, 2, 3, 4, 5, 6, 7),  # Re
                    (-4, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8),  # Os
                    (-3, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9),  # Ir
                    (-3, -2, -1, 1, 2, 3, 4, 5, 6),  # Pt
                    (-3, -2, -1, 1, 2, 3, 5),  # Au
                    (-2, 1, 2),  # Hg
                    (-5, -2, -1, 1, 2, 3),  # Tl
                    (-4, -2, -1, 1, 2, 3, 4),  # Pb
                    (-3, -2, -1, 1, 2, 3, 4, 5),  # Bi
                    (-2, 2, 4, 5, 6),  # Po
                    (-1, 1, 3, 5, 7),  # At
                    (2, 6),  # Rn
                    (1,),  # Fr
                    (2,),  # Ra
                    (3,),  # Ac
                    (1, 2, 3, 4),  # Th
                    (3, 4, 5),  # Pa
                    (1, 2, 3, 4, 5, 6),  # U
                    (2, 3, 4, 5, 6, 7),  # Np
                    (2, 3, 4, 5, 6, 7),  # Pu
                    (2, 3, 4, 5, 6, 7),  # Am
                    (3, 4, 6),  # Cm
                    (3, 4),  # Bk
                    (2, 3, 4),  # Cf
                    (2, 3, 4),  # Es
                    (2, 3),  # Fm
                    (2, 3),  # Md
                    (2, 3),  # No
                    (3,),  # Lr
                    (4,),  # Rf
                    (5,),  # Db
                    (6,),  # Sg
                    (7,),  # Bh
                    (8,),  # Hs
                    (),  # Mt
                    (),  # Ds
                    (),  # Rg
                    (2,),  # Cn
                    (),  # Nh
                    (),  # Fl
                    (),  # Mc
                    (),  # Lv
                    (),  # Ts
                    ()]  # Og

oxidation_states_common = [(),
                           (-1, 1),  # H
                           (),  # He
                           (1,),  # Li
                           (2,),  # Be
                           (3,),  # B
                           (-4, -3, -2, -1, 1, 2, 3, 4),  # C
                           (-3, 3, 5),  # N
                           (-2,),  # O
                           (-1,),  # F
                           (),  # Ne
                           (1,),  # Na
                           (2,),  # Mg
                           (3,),  # Al
                           (-4, 4),  # Si
                           (-3, 3, 5),  # P
                           (-2, 2, 4, 6),  # S
                           (-1, 1, 3, 5, 7),  # Cl
                           (),  # Ar
                           (1,),  # K
                           (2,),  # Ca
                           (3,),  # Sc
                           (4,),  # Ti
                           (5,),  # V
                           (3, 6),  # Cr
                           (2, 4, 7),  # Mn
                           (2, 3, 6),  # Fe
                           (2, 3),  # Co
                           (2,),  # Ni
                           (2,),  # Cu
                           (2,),  # Zn
                           (3,),  # Ga
                           (-4, 2, 4),  # Ge
                           (-3, 3, 5),  # As
                           (-2, 2, 4, 6),  # Se
                           (-1, 1, 3, 5),  # Br
                           (2,),  # Kr
                           (1,),  # Rb
                           (2,),  # Sr
                           (3,),  # Y
                           (4,),  # Zr
                           (5,),  # Nb
                           (4, 6),  # Mo
                           (4, 7),  # Tc
                           (3, 4),  # Ru
                           (3,),  # Rh
                           (2, 4),  # Pd
                           (1,),  # Ag
                           (2,),  # Cd
                           (3,),  # In
                           (-4, 2, 4),  # Sn
                           (-3, 3, 5),  # Sb
                           (-2, 2, 4, 6),  # Te
                           (-1, 1, 3, 5, 7),  # I
                           (2, 4, 6),  # Xe
                           (1,),  # Cs
                           (2,),  # Ba
                           (3,),  # La
                           (3, 4),  # Ce
                           (3,),  # Pr
                           (3,),  # Nd
                           (3,),  # Pm
                           (3,),  # Sm
                           (2, 3),  # Eu
                           (3,),  # Gd
                           (3,),  # Tb
                           (3,),  # Dy
                           (3,),  # Ho
                           (3,),  # Er
                           (3,),  # Tm
                           (3,),  # Yb
                           (3,),  # Lu
                           (4,),  # Hf
                           (5,),  # Ta
                           (4, 6),  # W
                           (4,),  # Re
                           (4,),  # Os
                           (3, 4),  # Ir
                           (2, 4),  # Pt
                           (3,),  # Au
                           (1, 2),  # Hg
                           (1, 3),  # Tl
                           (2, 4),  # Pb
                           (3,),  # Bi
                           (-2, 2, 4),  # Po
                           (-1, 1),  # At
                           (2,),  # Rn
                           (1,),  # Fr
                           (2,),  # Ra
                           (3,),  # Ac
                           (4,),  # Th
                           (5,),  # Pa
                           (6,),  # U
                           (5,),  # Np
                           (4,),  # Pu
                           (3,),  # Am
                           (3,),  # Cm
                           (3,),  # Bk
                           (3,),  # Cf
                           (3,),  # Es
                           (3,),  # Fm
                           (3,),  # Md
                           (2,),  # No
                           (3,),  # Lr
                           (4,),  # Rf
                           (5,),  # Db
                           (6,),  # Sg
                           (7,),  # Bh
                           (8,),  # Hs
                           (),  # Mt
                           (),  # Ds
                           (),  # Rg
                           (2,),  # Cn
                           (),  # Nh
                           (),  # Fl
                           (),  # Mc
                           (),  # Lv
                           (),  # Ts
                           ()]  # Og


def cpk_color(arg):
    return _get_property(cpk_colors, arg)


def _get_property(table, value=None, scale_factor=1):
    """
    Return a property from a given table of values.

    Args:
        value: (int,float,list) Value or list of values
           for the atomic numbers for which the
           property will be searched
        table: (list) List of properties for all atoms
           in the periodic table
           The index of the table must be equal to the atomic
           number of the element
        scale_factor: (float) Scaling factor for the
    """
    if hasattr(value, 'decode'):
        value = value.decode()

    if value is None:
        ret = {}
        for i in range(1, len(table)):
            ret[atomic_symbols[i]] = scale_factor * table[i]
    elif isinstance(value, int):
        ret = (scale_factor * table[value]) if value is not None else None
    elif isinstance(value, float):
        ret = (scale_factor * table[int(value)]) if value is not None else None
    elif isinstance(value, str) and value in atomic_symbols:
        if table[atomic_number(value)] is None:
            ret = float('nan')
        else:
            ret = scale_factor * table[atomic_number(value)]
    else:
        try:
            ret = [scale_factor * table[int(x)] for x in value]
        except ValueError:
            if not all([(x in atomic_symbols) for x in value]):
                raise ValueError('Not all the values are valid:', value)
            else:
                ret = [scale_factor * table[int(x)] for x in
                       atomic_number(value)]
    return ret


def valence(value=None):
    """
    Return the 'valences' as taken from
    http://www.ptable.com/#Property/Valence

    Args:
        value: (int,float,str,list) The value/s
               for which the valence will be evaluated
               Float values will be cast into integer

    Return:
        The valence of the element/s

    Examples:
    >>> valence(1)
    1
    >>> valence([1, 2])
    [1, 0]
    >>> valence(['H', 'He'])
    [1, 0]
    """
    return _get_property(valences, value)


def valence_nominal(value=None):
    """
    Return the 'valences' as taken from
    http://www.ptable.com/#Property/Valence

    Args:
        value: (int,float,str,list) The value/s
               for which the valence will be evaluated
               Float values will be cast into integer

    Return:
        The valence of the element/s

    Examples:
    >>> valence(1)
    1
    >>> valence([1, 2])
    [1, 0]
    >>> valence(['H', 'He'])
    [1, 0]
    """
    return _get_property(valences_nominal, value)


def period(value=None):
    # type: (str, int, list) -> (int, list)
    """
    Return the period of the element(s)

    :param value: atom symbol, atom number or list of atoms
    :return: The period of the atom or atoms
    :rtype: (int, list)

    Examples:
    >>> period('C')
    2
    >>> period(6)
    2
    >>> period(['Au', 'Sb'])
    [6, 5]
    >>> period(['Na', 'O', 'Ag', 'La'])
    [3, 2, 5, 6]
    >>> period(range(1, 12))
    [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3]

    """
    return _get_property(periods, value)


def group(value=None):
    # type: (str, int, list) -> (int, list)
    """
    Return the group of the element(s)
    The group for lanctanides and actinides is identified with -3

    :param value: Atom symbol, atom number or list of atoms
    :return: The group of the atom or atoms
    :rtype: (int, list)

    Examples:
    >>> group('C')
    14
    >>> group(6)
    14
    >>> group(['Au', 'Sb'])
    [11, 15]
    >>> group(['Na', 'O', 'Ag', 'La'])
    [1, 16, 11, -3]
    >>> group(range(1, 12))
    [1, 18, 1, 2, 13, 14, 15, 16, 17, 18, 1]

    """
    table = []
    for i in groups:
        if i is None:
            table.append(0)
        else:
            table.append(i)

    return _get_property(table, value)


def block(value=None):
    # type: (str, int, list) -> (str, list)
    """
    Return the orbital block of the element(s)
    For example, elements on group 1 and 2 are
    s-block. Transition metals are on the d-block,

    :param value: Atom symbol, atom number or list of atoms
    :return: The block of the atom or atoms
    :rtype: (str, list)

    Examples:
    >>> block('C')
    'p'
    >>> block(6)
    'p'
    >>> block(['Au', 'Sb'])
    ['d', 'p']
    >>> block(['Na', 'O', 'Ag', 'La'])
    ['s', 'p', 'd', 'd']
    >>> block(range(1, 12))
    ['s', 'p', 's', 's', 'p', 'p', 'p', 'p', 'p', 'p', 's']
    """
    table = []
    for i in blocks:
        if i is None:
            table.append(0)
        else:
            table.append(i)

    return _get_property(table, value)


def electronegativity(value=None):
    # type: (str, int, list) -> (float, list)
    """
    Electronegativity of the atom or atoms

    :param value: Atom symbol, atom number or list of atoms
    :return: The value of electronegativity for the atom or atoms
    :rtype: (float, list)

    Examples:
    >>> electronegativity('C')
    2.55
    >>> electronegativity(6)
    2.55
    >>> electronegativity(['Au', 'Sb'])
    [2.54, 2.05]
    >>> electronegativity(['Na', 'O', 'Ag', 'La'])
    [0.93, 3.44, 1.93, 1.1]
    >>> electronegativity(range(1, 12))
    [2.2, 0, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0, 0.93]

    """
    table = []
    for i in electronegativities:
        if i is None:
            table.append(0)
        else:
            table.append(i)

    return _get_property(table, value)


def covalent_radius(value=None):
    # type: (str, int, list) -> (float, list)
    """
    Covalent radius in Angstrom of an atom or list of atoms.
    The atoms could be given as symbols or by atomic number

    :param value: Atom symbol, atom number or list of atoms
    :return: Covalent radius in angstrom for the atom or atoms
    :rtype: float

    Examples:
    >>> covalent_radius('C')
    0.76
    >>> covalent_radius(6)
    0.76
    >>> covalent_radius(['Au', 'Sb'])
    [1.36, 1.39]
    >>> covalent_radius(['Na', 'O', 'Ag', 'La'])
    [1.66, 0.66, 1.45, 2.07]
    >>> covalent_radius(range(1, 12))
    [0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66]

    """

    return _get_property(covalent_radii, value)


def atomic_symbol(value=None):
    # type: (int, float, list) -> (str, list)
    """
    Atomic symbol for a atom or atoms given by atomic number
    Float values will be rounded to the lowest integer

    :param value: Atomic number or list of atomic number
    :return: Atomic symbol(s) for the atom(s)
    :rtype: (str, list)

    Examples:
    >>> atomic_symbol(6)
    'C'
    >>> atomic_symbol(6.9)
    'C'
    >>> atomic_symbol([79, 51])
    ['Au', 'Sb']
    >>> atomic_symbol([11, 8, 47, 57])
    ['Na', 'O', 'Ag', 'La']
    >>> atomic_symbol(range(1, 12))
    ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na']

    """
    ret = None
    if value is None:
        ret = {}
        for i in range(1, len(atomic_symbols)):
            ret[atomic_symbols[i]] = i
    elif isinstance(value, int):
        ret = atomic_symbols[value]
    elif isinstance(value, float):
        ret = atomic_symbols[int(value)]
    elif _np.iterable(value):
        try:
            ret = [atomic_symbols[int(x)] for x in value]
        except ValueError:
            ret = [atomic_symbols[int(x)] for x in atomic_number(value)]
    return ret


def mass(value=None):
    # type: (int, str, list) -> (float, list)
    """
    Atomic mass for a atom or atoms given by atomic number
    Float values will be rounded to the lowest integer
    Atomic masses are returned in Standard Atomic weight

    :param value: Atomic mass or list of atomic masses
    :return: Atomic mass(es) for the atom(s)
    :rtype: float

    Examples:
    >>> mass(6)
    12.011
    >>> mass(6.9)
    12.011
    >>> mass(['Au', 'Sb'])
    [196.96654, 121.753]
    >>> mass(['Na', 'O', 'Ag', 'La'])
    [22.989768, 15.9994, 107.8682, 138.9055]
    >>> mass(range(1, 12))
    [1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.011, 14.00674, 15.9994, 18.9984032, 20.1797, 22.989768]

    """
    return _get_property(masses, value)


def atomic_number(arg):
    # type: (str, list) -> (int, list)
    """
    Atomic number(s) of a symbol or list of atomic symbols

    :param arg: Atomic number or list of atomic numbers
    :return: Atomic number(s) for the atom(s)
    :rtype: int, list

    Examples:
    >>> atomic_number('C')
    6
    >>> atomic_number(['Au', 'Sb'])
    [79, 51]
    >>> atomic_number(['Na', 'O', 'Ag', 'La'])
    [11, 8, 47, 57]

    """
    symbol_dict = dict(atomic_symbol())
    if hasattr(arg, 'decode'):
        arg = arg.decode()
    if isinstance(arg, str):
        if arg not in atomic_symbols:
            raise ValueError('Atomic symbol not found')
        return symbol_dict[str(arg)]
    try:
        return [atomic_number(x) for x in arg]
    except TypeError:  # catch when for loop fails
        raise ValueError('Argument not recognized as atomic symbol or list of atomic symbols')


def atomic_name(value):
    return _get_property(atomic_names, value)


def melting_point(value):
    return _get_property(melting_points, value)


def boiling_point(value):
    return _get_property(boiling_points, value)


def crystal_structure(value):
    return _get_property(crystal_structures, value)


def phase(value):
    return _get_property(phases, value)


def oxidation_state(value, common=False):
    if common:
        return _get_property(oxidation_states_common, value)
    else:
        return _get_property(oxidation_states, value)
