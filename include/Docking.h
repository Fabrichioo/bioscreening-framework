#ifndef DOCKING_H
#define DOCKING_H

#include "Molecule.h"

// Función básica de docking que compara una proteína y un ligando.
// Devuelve un score numérico (se usa el dummy implementado en performDocking)
float performDocking(const Molecule& protein, const Molecule& ligand);

// Estructura para almacenar el resultado del docking: 
// índice de proteína, índice de ligando y score de docking.
struct DockingResult {
    int proteinIndex;
    int ligandIndex;
    float score;
};

#endif // DOCKING_H