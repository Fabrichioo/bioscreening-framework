#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "DataManager.h"
#include "Molecule.h"
#include "Docking.h"
#include "Utils.h"

/* src/Parallel/MPI_OpenMP.cpp */
// #include "Parallel.h"
// #include "Docking.h"    // Para performDocking
// #include "Molecule.h"
// #include <mpi.h>
// #ifdef _OPENMP
//   #include <omp.h>
// #endif
// #include <vector>
// #include <iostream>
// #include <cstdlib>

void hybridDocking(const std::vector<Molecule>& proteins,
                   const std::vector<Molecule>& ligands,
                   std::vector<float>& scores) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        std::cout << "Ejecutando docking híbrido con MPI + OpenMP..." << std::endl;

    // Calcular el total de evaluaciones: (proteínas × ligandos).
    size_t total = proteins.size() * ligands.size();
    // Distribuir el total entre procesos.
    size_t chunk = total / size;
    size_t remainder = total % size;
    size_t start = rank * chunk + (rank < remainder ? rank : remainder);
    size_t end = start + chunk + (rank < remainder ? 1 : 0);

    // Cada proceso calculará sus docking en el rango [start, end)
    std::vector<float> localScores(end - start);

    // Utilizar OpenMP para paralelizar el loop en el proceso.
    #pragma omp parallel for schedule(dynamic)
    for (size_t idx = start; idx < end; ++idx) {
        int i = idx / ligands.size();   // índice de la proteína
        int j = idx % ligands.size();     // índice del ligando
        localScores[idx - start] = performDocking(proteins[i], ligands[j]);
    }

    // Recolectar los tamaños locales de cada proceso.
    int localSize = localScores.size();
    std::vector<int> recvCounts;
    std::vector<int> displs;
    if (rank == 0) {
        recvCounts.resize(size);
    }
    MPI_Gather(&localSize, 1, MPI_INT,
               recvCounts.data(), 1, MPI_INT,
               0, MPI_COMM_WORLD);

    // El proceso maestro prepara el vector final de resultados.
    if (rank == 0) {
        displs.resize(size);
        displs[0] = 0;
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + recvCounts[i - 1];
        }
        scores.resize(total);
    }

    // Se recopilan los resultados locales en el vector scores en el proceso 0.
    MPI_Gatherv(localScores.data(), localSize, MPI_FLOAT,
                scores.data(), recvCounts.data(), displs.data(), MPI_FLOAT,
                0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[]) {

    DataManager dataManager;
    std::vector<Molecule> proteins;
    std::vector<Molecule> ligands;

    if (!dataManager.loadProteins("data/proteins/", proteins)) {
        std::cerr << "Error cargando proteínas." << std::endl;
        return 1;
    }
    if (!dataManager.loadLigands("data/ligands/", ligands)) {
        std::cerr << "Error cargando ligandos." << std::endl;
        return 1;
    }

    Timer timer;
    timer.start();

    std::cout << "Sequential Mode" << std::endl;
    std::vector<float> scores;
    for (const auto &protein : proteins) {
        for (const auto &ligand : ligands) {
            float score = performDocking(protein, ligand);
            scores.push_back(score);
        }
    }

    timer.stop();
    std::cout << "Tiempo de ejecución: " << timer.elapsedMilliseconds() << " ms" << std::endl;

    analyzeDockingResults(scores, proteins.size(), ligands.size());
    
    exit(EXIT_SUCCESS);
}