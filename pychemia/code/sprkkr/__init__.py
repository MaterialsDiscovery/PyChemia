import pychemia
import datetime


class CodeSPRKKR:

    def __init__(self, structure):
        self.structure = structure

    def write_sys(self):
        name = self.structure.formula
        wf = open(name + '.sys', 'w')
        wf.write('system data-file created by pychemia on %s\n' % str(datetime.datetime.now()))
        wf.write(name + '.sys\n')
        wf.write('xband-version\n')
        wf.write('6.3\n')
        wf.write('dimension\n')
        wf.write('3D\n')
        sym = pychemia.crystal.CrystalSymmetry(self.structure)
        wf.write("Bravais lattice \n")
        wf.write('12 %s primitive %s %s\n' % (sym.crystal_system().lower(), sym.symbol(),
                                              sym.get_symmetry_dataset()['pointgroup']))
        wf.write('space group number (ITXC and AP)\n')
        wf.write(' %4d %4d\n' % (0, 0))
        wf.write('structure type\n')
        wf.write('UNKNOWN\n')
        wf.write('lattice parameter A  [a.u.] \n')
        a, b, c = self.structure.lattice.lengths
        au = pychemia.utils.constants.angstrom_bohr
        wf.write(' %17.12f\n' % (a * au))
        wf.write('ratio of lattice parameters  b/a  c/a\n')
        wf.write(" %17.12f %17.12f\n" % (b/a, c/a))
        wf.write("lattice parameters  a b c  [a.u.] \n")
        wf.write(" %17.12f %17.12f %17.12f\n" % (a*au, b*au, c*au))
        alpha, beta, gamma = self.structure.lattice.angles
        wf.write("lattice angles  alpha beta gamma  [deg] \n")
        wf.write(' %17.12f %17.12f %17.12f\n' % (alpha, beta, gamma))
        wf.write("primitive vectors     (cart. coord.) [A]\n")
        wf.write(" %17.12f %17.12f %17.12f\n" % tuple(self.structure.lattice.cell[0]/a))
        wf.write(" %17.12f %17.12f %17.12f\n" % tuple(self.structure.lattice.cell[1]/a))
        wf.write(" %17.12f %17.12f %17.12f\n" % tuple(self.structure.lattice.cell[2]/a))
        wf.write("number of sites NQ\n")
        wf.write(" %2d\n" % self.structure.natom)
        wf.write(" IQ ICL     basis vectors     (cart. coord.) [A] %s RWS [a.u.]  NLQ  NOQ ITOQ\n" % (20*' '))

        dm = self.structure.distance_matrix().flatten()
        dmin = min(dm[dm > 0.0]) * pychemia.utils.constants.angstrom_bohr
        rws = dmin / 2.0
        reduced = self.structure.reduced
        for i in range(self.structure.natom):
            wf.write(" %2d %3d %17.12f %17.12f %17.12f %17.12f %2d %2d %d\n" % (i+1, i+1,
                                                                                reduced[i, 0],
                                                                                reduced[i, 1],
                                                                                reduced[i, 2],
                                                                                rws, 3, 1, i+1))

        wf.write("number of sites classes NCL\n")
        wf.write(" %2d\n" % self.structure.natom)
        wf.write("ICL WYCK NQCL IQECL (equivalent sites)\n")
        for i in range(self.structure.natom):
            wf.write(" %2d   -   %2d %2d\n" % (i+1, 1, i+1))
        wf.write("number of atom types NT\n")
        wf.write(" %2d\n" % self.structure.natom)
        wf.write(" IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n")
        for i in range(self.structure.natom):
            atomic_z = pychemia.utils.periodic.atomic_number(self.structure.symbols[i])
            wf.write(" %2d %3d %6s %7d %5.3f %2d\n" % (i+1, atomic_z, self.structure.symbols[i], 1, 1.0, i+1))
        wf.close()

    def write_pot(self):
        strdate = str(datetime.datetime.now())
        natom = self.structure.natom
        name = self.structure.formula
        ret = """*******************************************************************************
HEADER    'SCF-start data created by pychemia %s'
*******************************************************************************
TITLE     'SPR-KKR calculation for %s'
SYSTEM    %s
PACKAGE   SPRKKR
FORMAT    6  (21.05.2007)
*******************************************************************************
GLOBAL SYSTEM PARAMETER
NQ                 %d
NT                 %d
NM                 %d
IREL               3
*******************************************************************************
SCF-INFO
INFO      NONE
SCFSTATUS START
FULLPOT   F
BREITINT  F
NONMAG    F
ORBPOL    NONE
EXTFIELD  F
BLCOUPL   F
BEXT          0.0000000000
SEMICORE  F
LLOYD     F
NE                32
IBZINT             2
NKTAB              0
XC-POT    VWN
SCF-ALG   BROYDEN2
SCF-ITER           0
SCF-MIX       0.2000000000
SCF-TOL       0.0000100000
RMSAVV    999999.0000000000
RMSAVB    999999.0000000000
EF            0.0000000000
VMTZ          0.0000000000
""" % (strdate, name, name, natom, natom, natom)

        sym = pychemia.crystal.CrystalSymmetry(self.structure)
        a, b, c = self.structure.lattice.lengths
        au = pychemia.utils.constants.angstrom_bohr

        ret += """*******************************************************************************
LATTICE
SYSDIM       3D
SYSTYPE      BULK
BRAVAIS           12        %s       primitive      %s    %s
""" % (sym.crystal_system().lower(), sym.symbol(), sym.get_symmetry_dataset()['pointgroup'])

        ret += "ALAT          %15.10f\n" % (a*au)

        ret += "A(1)          %15.10f %15.10f %15.10f\n" % tuple(self.structure.lattice.cell[0] / a)
        ret += "A(2)          %15.10f %15.10f %15.10f\n" % tuple(self.structure.lattice.cell[1] / a)
        ret += "A(3)          %15.10f %15.10f %15.10f\n" % tuple(self.structure.lattice.cell[2] / a)

        ret += """*******************************************************************************
SITES
CARTESIAN T
BASSCALE      1.0000000000    1.0000000000    1.0000000000
        IQ      QX              QY              QZ
"""
        reduced = self.structure.reduced
        for i in range(natom):
            ret += "%d  %15.10f %15.10f %15.10f\n" % (i+1, reduced[i, 0], reduced[i, 1], reduced[i, 2])

        ret += """*******************************************************************************
OCCUPATION
        IQ     IREFQ       IMQ       NOQ  ITOQ  CONC
"""

        for i in range(natom):
            ret += "%10d %10d %10d %10d %5d %5f\n" % (i+1, i+1, i+1, 1, i+1, 1.0)

        ret += """*******************************************************************************
REFERENCE SYSTEM
NREF              %d
      IREF      VREF            RMTREF
""" % natom

        for i in range(natom):
            ret += " %10d %15.10f %15.10f\n" % (i+1, 4.0, 0.0)

        ret += """*******************************************************************************
MAGNETISATION DIRECTION
KMROT              0
QMVEC         0.0000000000    0.0000000000    0.0000000000
        IQ      QMTET           QMPHI
"""
        for i in range(natom):
            ret += " %10d %15.10f %15.10f\n" % (i+1, 0.0, 0.0)

        ret += """*******************************************************************************
MESH INFORMATION
MESH-TYPE EXPONENTIAL
   IM      R(1)            DX         JRMT      RMT        JRWS      RWS
        """

        for i in range(natom):
            ret += " %5d %15.10f %15.10f %d %15.10f %d %15.10f\n" % (i+1, 0.0000010000, 0.0205501520,    0,
                                                                     2.2661445689, 721, 2.6660524340)

        ret += """*******************************************************************************
TYPES
   IT     TXTT        ZT     NCORT     NVALT    NSEMCORSHLT
"""
        for i in range(natom):
            atomic_z = pychemia.utils.periodic.atomic_number(self.structure.symbols[i])
            ret += " %2d %6s %5d %7d %5d %2d\n" % (i+1, self.structure.symbols[i], atomic_z, 18, 8-i, 0)

        wf = open(name + '.pot', 'w')
        wf.write(ret)
        wf.close()
