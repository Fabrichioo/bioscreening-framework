#ifndef PARALLEL_H
#define PARALLEL_H

#include <vector>
#include "Molecule.h"

// Cada función paralela realiza el docking para todos los pares proteína-ligando 
// y guarda los resultados (scores) en el vector "scores".

#ifdef _OPENMP
#include <omp.h>
#endif

void openmpDocking(const std::vector<Molecule>& proteins,
                   const std::vector<Molecule>& ligands,
                   std::vector<float>& scores);

void mpiDocking(const std::vector<Molecule>& proteins,
                const std::vector<Molecule>& ligands,
                std::vector<float>& scores);

#ifdef __CUDACC__
void cudaDocking(const std::vector<Molecule>& proteins,
                 const std::vector<Molecule>& ligands,
                 std::vector<float>& scores);
#endif

#endif // PARALLEL_H