#include "Parallel.h"
#include "Docking.h"
#include <mpi.h>
#include <vector>
#include <iostream>

void mpiDocking(const std::vector<Molecule>& proteins,
                const std::vector<Molecule>& ligands,
                std::vector<float>& scores) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0)
        std::cout << "Ejecutando docking con MPI..." << std::endl;

    size_t total = proteins.size() * ligands.size();
    // Distribución del trabajo entre procesos
    size_t chunk = total / size;
    size_t remainder = total % size;

    size_t start = rank * chunk + (rank < remainder ? rank : remainder);
    size_t end = start + chunk + (rank < remainder ? 1 : 0);

    std::vector<float> localScores(end - start);

    for (size_t idx = start; idx < end; ++idx) {
        size_t i = idx / ligands.size();
        size_t j = idx % ligands.size();
        localScores[idx - start] = performDocking(proteins[i], ligands[j]);
    }

    if(rank == 0) {
        scores.resize(total);
    }

    // Preparar recuento y desplazamientos para la recolección de datos
    std::vector<int> recvCounts(size);
    std::vector<int> displs(size);
    int localSize = localScores.size();
    MPI_Gather(&localSize, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + recvCounts[i - 1];
        }
    }

    MPI_Gatherv(localScores.data(), localSize, MPI_FLOAT,
                scores.data(), recvCounts.data(), displs.data(), MPI_FLOAT, 0, MPI_COMM_WORLD);
}