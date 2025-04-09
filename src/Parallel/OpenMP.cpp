#include "Parallel.h"
#include "Docking.h"
#include <omp.h>
#include <vector>
#include <iostream>

void openmpDocking(const std::vector<Molecule>& proteins,
                   const std::vector<Molecule>& ligands,
                   std::vector<float>& scores) {
    std::cout << "Ejecutando docking con OpenMP..." << std::endl;
    size_t total = proteins.size() * ligands.size();
    scores.resize(total);

    // Paralelización doble: cada iteración se ejecuta en paralelo
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < proteins.size(); ++i) {
        for (size_t j = 0; j < ligands.size(); ++j) {
            size_t idx = i * ligands.size() + j;
            scores[idx] = performDocking(proteins[i], ligands[j]);
        }
    }
}