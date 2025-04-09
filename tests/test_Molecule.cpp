#include "Molecule.h"
#include <cassert>
#include <iostream>

int main() {
    Molecule mol;
    
    Atom atom1 = {0.0f, 0.0f, 0.0f, "C"};
    Atom atom2 = {1.0f, 1.0f, 1.0f, "O"};

    mol.addAtom(atom1);
    mol.addAtom(atom2);

    assert(mol.getAtoms().size() == 2);
    std::cout << "test_Molecule passed." << std::endl;
    return 0;
}