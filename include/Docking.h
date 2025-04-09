#ifndef DOCKING_H
#define DOCKING_H

#include "Molecule.h"

// Función básica de docking que compara una proteína y un ligando.
// Devuelve un score numérico (a modo de ejemplo, un scoring dummy)
float performDocking(const Molecule& protein, const Molecule& ligand);

// Estructura para almacenar el resultado del docking (opcional)
struct DockingResult {
    float score;
    // Se podrían agregar otros campos, por ejemplo: orientación, posición, etc.
};

#endif // DOCKING_H
