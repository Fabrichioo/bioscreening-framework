#include "Parallel.h"
#include "Molecule.h"
#include <vector>
#include <cassert>
#include <iostream>

int main() {
    std::vector<Molecule> proteins(1);
    std::vector<Molecule> ligands(1);
    
    proteins[0].addAtom({0, 0, 0, "C"});
    ligands[0].addAtom({1, 0, 0, "H"});
    
    std::vector<float> scores;
    cudaDocking(proteins, ligands, scores);

    // El kernel CUDA dummy asigna 42.0f a cada score.
    assert(scores.size() == 1);
    assert(scores[0] == 42.0f);
    std::cout << "test_CUDA passed." << std::endl;
    return 0;
}