from __future__ import print_function
from builtins import str
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
                'Darmstadtium', 'Roentgenium', 'Copernicium', 'Ununtrium', 'Flerovium', 'Ununpentium', 'Livermorium',
                'Ununseptium', 'Ununoctium']

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
                  'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']

valences = [0,
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
            4, 5, 6, 7, 8]

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
          'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f',
          'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p',
          's', 's',
          'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f',
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

covalent_radii = [0.2, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66,
                  0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02,
                  1.06, 2.03, 1.76, 1.7, 1.6, 1.53, 1.39, 1.39, 1.32,
                  1.26, 1.24, 1.32, 1.22, 1.22, 1.2, 1.19, 1.2, 1.2,
                  1.16, 2.2, 1.95, 1.9, 1.75, 1.64, 1.54, 1.47, 1.46,
                  1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39,
                  1.4, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98,
                  1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.9, 1.87, 1.87,
                  1.75, 1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,
                  1.45, 1.46, 1.48, 1.4, 1.5, 1.5, 2.6, 2.21, 2.15,
                  2.06, 2., 1.96, 1.9, 1.87, 1.8, 1.69, 0.2, 0.2,
                  0.2, 0.2, 0.2, 0.2, 0.2]

masses = [0.0, 1.00794, 4.002602, 6.941, 9.012182,
          10.811, 12.011, 14.00674, 15.9994, 18.9984032,
          20.1797, 22.989768, 24.3050, 26.981539, 28.0855,
          30.973762, 32.066, 35.4527, 39.948, 39.0983,
          40.078, 44.955910, 47.88, 50.9415, 51.9961,
          54.93805, 55.847, 58.93320, 58.69, 63.546,
          65.39, 69.723, 72.61, 74.92159, 78.96,
          79.904, 83.80, 85.4678, 87.62, 88.90585,
          91.224, 92.90638, 95.94, 98.9062, 101.07,
          102.9055, 106.42, 107.8682, 112.411, 114.82,
          118.710, 121.753, 127.60, 126.90447, 131.29,
          132.90543, 137.327, 138.9055, 140.115, 140.90765,
          144.24, 147.91, 150.36, 151.965, 157.25,
          158.92534, 162.50, 164.93032, 167.26, 168.93421,
          173.04, 174.967, 178.49, 180.9479, 183.85,
          186.207, 190.2, 192.22, 195.08, 196.96654,
          200.59, 204.3833, 207.2, 208.98037, 209.0,
          210.0, 222.0, 223.0, 226.0254, 230.0,
          232.0381, 231.0359, 238.0289, 237.0482, 242.0,
          243.0, 247.0, 247.0, 249.0, 254.0,
          253.0, 256.0, 254.0, 257.0, 260.0]

cpk_colors = [None, (1, 1, 1),  # Hydrogen
              (0.0, 1.0, 1.0),  # Helium
              (0.4666666666666667, 0.0, 1.0),  # Lithium
              (0.0, 0.4666666666666667, 0.0),  # Berilium
              (1.0, 0.6666666666666666, 0.4666666666666667),  # Boron
              (0.1333333333333333, 0.133333333333333, 0.1333333333333333),  # Carbon
              (0.5294117647058824, 0.807843137254902, 0.9215686274509803),  # Nitrogen
              (1.0, 0.13333333333333333, 0.0),  # Oxygen
              (0.12156862745098039, 0.9411764705882353, 0.12156862745098039),  # Fluorine ]
              (0.0, 1.0, 1.0),  # Neon
              (0.4666666666666667, 0.0, 1.0),  # Sodium
              (0.0, 0.4666666666666667, 0.0),  # Magnesium
              (0.8666666666666667, 0.4666666666666667, 1.0),
              (0.8666666666666667, 0.4666666666666667, 1.0),
              (1.0, 0.6, 0.0),  # Phosphorus
              (0.8666666666666667, 0.4666666666666667, 1.0),
              (0.12156862745098039, 0.9411764705882353, 0.12156862745098039),  # Clorine ]
              (0.0, 1.0, 1.0),  # Argon
              (0.4666666666666667, 0.0, 1.0),  # Potasium
              (0.0, 0.4666666666666667, 0.0),  # Calsium
              (0.8666666666666667, 0.4666666666666667, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0),
              (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]


def cpk_color(arg):
    return _get_property(cpk_colors, arg)


def _get_property(table, value=None, scale_factor=1):
    """
    Return a property from a given table of values.

    Args:
        Z: (int,float,list) Value or list of values
           for the atomic numbers for which the
           property will be searched
        table: (list) List of properties for all atoms
           in the periodic table
           The index of the table must be equal to the atomic
           number of the element
        scaling factor: (float) Scaling factor for the
    """
    if hasattr(value, 'decode'):
        value = value.decode()

    ret = None
    if value is None:
        ret = {}
        for i in range(1, len(table)):
            ret[atomic_symbols[i]] = scale_factor * table[i]
    elif isinstance(value, int):
        ret = (scale_factor * table[value]) if value is not None else None
    elif isinstance(value, float):
        ret = (scale_factor * table[int(value)]) if value is not None else None
    elif isinstance(value, str) and value in atomic_symbols:
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
    ['s', 'p', 'd', 'f']
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
