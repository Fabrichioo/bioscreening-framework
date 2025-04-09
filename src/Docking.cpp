#include "Docking.h"
#include <cmath>

// Implementación dummy: el score se calcula a partir de la diferencia de número de átomos
float performDocking(const Molecule& protein, const Molecule& ligand) {
    int diff = std::abs(static_cast<int>(protein.getAtoms().size()) -
                        static_cast<int>(ligand.getAtoms().size()));
    return static_cast<float>(100 - diff);
}