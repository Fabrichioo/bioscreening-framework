#include "DataManager.h"
#include "Molecule.h"
#include <vector>
#include <cassert>
#include <iostream>

int main() {
    DataManager dm;
    std::vector<Molecule> proteins;
    std::vector<Molecule> ligands;

    // Aunque la función es dummy, esperamos que cargue al menos una molécula.
    bool resProteins = dm.loadProteins("data/proteins/", proteins);
    bool resLigands = dm.loadLigands("data/ligands/", ligands);

    assert(resProteins);
    assert(resLigands);
    assert(!proteins.empty());
    assert(!ligands.empty());
    std::cout << "test_DataManager passed." << std::endl;
    return 0;
}