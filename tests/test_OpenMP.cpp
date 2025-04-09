#include "Parallel.h"
#include "Molecule.h"
#include "Docking.h"
#include <vector>
#include <cassert>
#include <iostream>

int main() {
    // Creamos datos de prueba: 2 proteínas y 2 ligandos
    std::vector<Molecule> proteins(2);
    std::vector<Molecule> ligands(2);

    // Cada molécula tendrá 1 átomo, de modo que la diferencia es 0 y score = 100.
    proteins[0].addAtom({0, 0, 0, "C"});
    proteins[1].addAtom({1, 0, 0, "O"});
    ligands[0].addAtom({0, 1, 0, "N"});
    ligands[1].addAtom({0, 0, 1, "H"});

    std::vector<float> scores;
    openmpDocking(proteins, ligands, scores);

    // Se evalúan 4 pares (2 x 2) y se espera score = 100 para cada uno.
    assert(scores.size() == 4);
    for (auto s : scores) {
        assert(s == 100.0f);
    }
    std::cout << "test_OpenMP passed." << std::endl;
    return 0;
}