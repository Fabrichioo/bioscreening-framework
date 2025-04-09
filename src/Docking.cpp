#include "Docking.h"
#include <cmath>

// Re-implementación real de performDocking basada en un potencial de Lennard-Jones
float performDocking(const Molecule& protein, const Molecule& ligand) {
    // Verificar que ambas moléculas contienen átomos.
    if (protein.getAtoms().empty() || ligand.getAtoms().empty()) {
        return 0.0f;
    }
    
    float energy = 0.0f;
    // Parámetros del potencial Lennard-Jones (para simplificar, sigma = 1 y epsilon = 1)
    const float epsilon = 1.0f;
    // sigma^6 y sigma^12 son 1 cuando sigma = 1, por lo que se omiten.
    
    // Se recorre cada par de átomos (ligando-proteína)
    for (const auto& atomL : ligand.getAtoms()) {
        for (const auto& atomP : protein.getAtoms()) {
            float dx = atomL.x - atomP.x;
            float dy = atomL.y - atomP.y;
            float dz = atomL.z - atomP.z;
            float r2 = dx * dx + dy * dy + dz * dz;
            // Evitar división por cero o distancias extremadamente cortas
            if (r2 < 1e-6f)
                continue;
            // Calcular r^6 y r^12
            float r6 = r2 * r2 * r2;
            float r12 = r6 * r6;
            // Cálculo del potencial: 4 * epsilon * ((1 / r^12) - (1 / r^6))
            float ljEnergy = 4.0f * epsilon * ((1.0f / r12) - (1.0f / r6));
            energy += ljEnergy;
        }
    }
    return energy;
}