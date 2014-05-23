# Real physical constants
# Revised fundamental constants from http://physics.nist.gov/cuu/Constants/index.html
# (from 2006 least squares adjustment)

#: 1 Bohr, in Angstrom
bohr_angstrom = 0.52917720859
angstrom_bohr = 1.0 / bohr_angstrom

#: 1 Hartree, in cm^-1
Ha_cmm1 = 219474.6313705

#: 1 Hartree, in eV
Ha_eV = 27.21138386
eV_Ha = 1 / Ha_eV

#: 1 Hartree, in meV
Ha_meV = Ha_eV * 1000

#: 1 Hartree, in Kelvin
Ha_K = 315774.65

# 1 Hartree, in THz
Ha_THz = 6579.683920722

#: 1 Hartree, in J
Ha_J = 4.35974394E-18

#:  Minus the electron charge, in Coulomb
e_Cb = 1.602176487E-19

#: Boltzmann constant in Ha/K
kb_HaK = 8.617343E-5 / Ha_eV

#: 1 atomic mass unit, in electronic mass
amu_emass = 1.660538782E-27 / 9.10938215E-31

#: 1 Ha/Bohr^3, in GPa
HaBohr3_GPa = Ha_eV / bohr_angstrom ** 3 * e_Cb * 1.0E+21

#: Avogadro number
Avogadro = 6.02214179E+23

#: Inverse of fine structure constant
InvFineStruct = 137.035999679

#: Speed of light in atomic units
Sp_Lt = 2.99792458E+8 / 2.1876912633E+6

#: Atomic unit of time, in seconds
Time_Sec = 2.418884326505E-17

#: Atomic unit of induction field (in Tesla) * mu_B (in atomic units).
BField_Tesla = 4.254383E-6 * 0.5
