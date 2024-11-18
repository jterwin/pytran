
from __future__ import print_function


import pytran
import numpy as np


if __name__ == '__main__':

    print("\ntest get_molecule_id:")
    for mol in ('CO2', 'CH4', 'CO', 'HCN', 'H2O'):
        print("\t%s: %d" % (mol, pytran.get_molecule_id(mol)))

    print("\ntest get_iso_id:")
    for name in ('CO2_626', 'CO2_838', 'CO_38', 'C2H6_1221', 'HCN_124'):
        nmol, niso = pytran.get_iso_id(name)
        print("\t%s: %d,%d" % (name, nmol, niso))

    print("\ntest get_iso_name:")
    for nmol, niso in zip((1,2,3,4,5,6),(1,2,1,2,1,2)):
        print("\t%d, %d: %s" % (nmol, niso, pytran.get_iso_name(nmol, niso)))
 
    print("\ntest get_molecule_nisops:")
    for nmol in (1,2,32,44):
        print("\t%d: %d" % (nmol,pytran.get_molecule_nisops(nmol)))

    print("\ntest get_iso_mass:")
    for nmol, niso in zip((1,2,3,4,5,6),(1,2,1,2,1,2)):
        print("\t%d, %d: %.4f" % (nmol, niso, pytran.get_iso_mass(nmol, niso)))
    
    print("\ntest get_molecule_mass:")
    for nmol in (1,2,32,44):
        print("\t%d: %.4f" % (nmol,pytran.get_molecule_mass(nmol)))

    print("\ntest qtips")
    for tmp in np.logspace(0,3,13):
        print("\t%.1f: %.1f" % (tmp, pytran.qtips(tmp, 2, 1)))

        

