#include "Docking.h"
#include "Molecule.h"
#include <cassert>
#include <iostream>

int main() {
    Molecule protein;
    Molecule ligand;

    // Creamos proteína con 3 átomos y ligando con 2 átomos.
    protein.addAtom({0, 0, 0, "N"});
    protein.addAtom({1, 0, 0, "C"});
    protein.addAtom({0, 1, 0, "O"});
    
    ligand.addAtom({0, 0, 1, "H"});
    ligand.addAtom({1, 1, 1, "C"});

    float score = performDocking(protein, ligand);
    //Según la función dummy: score = 100 - |3 - 2| = 99
    assert(score == 99.0f);
    std::cout << "test_Docking passed." << std::endl;
    return 0;
}