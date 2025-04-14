#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "DataManager.h"
#include "Molecule.h"
#include "Docking.h"
#include "Utils.h"
#include <mpi.h>

void mpi_docking(const std::vector<Molecule>& proteins,
                 const std::vector<Molecule>& ligands,
                 std::vector<float>& scores) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "Ejecutando docking con MPI..." << std::endl;
        std::cout << "Número total de procesos: " << size << std::endl;
    }

    size_t total = proteins.size() * ligands.size();
    size_t chunk = total / size;
    size_t remainder = total % size;

    // Calcular el rango de índices a procesar por cada proceso
    size_t start = rank * chunk + (rank < remainder ? rank : remainder);
    size_t end = start + chunk + (rank < remainder ? 1 : 0);

    std::vector<float> localScores(end - start);
    for (size_t idx = start; idx < end; ++idx) {
        size_t i = idx / ligands.size();
        size_t j = idx % ligands.size();
        localScores[idx - start] = performDocking(proteins[i], ligands[j]);
    }

    // El proceso 0 reserva el espacio para almacenar todos los resultados
    if (rank == 0) {
        scores.resize(total);
    }

    // Preparar arreglos para el envío de datos
    std::vector<int> recvCounts(size);
    std::vector<int> displs(size);
    int localSize = static_cast<int>(localScores.size());
    MPI_Gather(&localSize, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + recvCounts[i - 1];
        }
    }

    MPI_Gatherv(localScores.data(), localSize, MPI_FLOAT,
                scores.data(), recvCounts.data(), displs.data(), MPI_FLOAT, 0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    std::string proteinsDir;
    std::string ligandsDir;
    bool verbose;
    
    parseArguments(argc, argv, proteinsDir, ligandsDir, verbose);

    DataManager dataManager;
    std::vector<Molecule> proteins;
    std::vector<Molecule> ligands;

    if (!dataManager.loadProteins(proteinsDir, proteins)) {
        std::cerr << "Error loading proteins." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!dataManager.loadLigands(ligandsDir, ligands)) {
        std::cerr << "Error loading ligands." << std::endl;
        exit(EXIT_FAILURE);
    }

    double t1, t2;
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    
    std::vector<float> scores;
    mpi_docking(proteins, ligands, scores);
    
    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        std::cout << "Tiempo de ejecución: " << t2 - t1 << " s" << std::endl;
    if (rank == 0 && verbose)
        analyzeDockingResults(scores, proteins.size(), ligands.size());
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}
