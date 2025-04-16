#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "DataManager.h"
#include "Molecule.h"
#include "Docking.h"
#include "Utils.h"
#include <mpi.h>
#include <omp.h>

std::vector<float> hybrid_docking(const std::vector<Molecule>& proteins, 
                                  const std::vector<Molecule>& ligands) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t total = proteins.size() * ligands.size();

    // Distribuir trabajo entre procesos MPI
    size_t chunk     = total / size;
    size_t remainder = total % size;
    size_t start = rank * chunk + (rank < remainder ? rank : remainder);
    size_t end   = start + chunk + (rank < remainder ? 1 : 0);
    
    std::vector<float> localScores(end - start);
    
    #pragma omp parallel
    {
        #pragma omp master
        {
            std::cout << "Process " << rank << " running hybrid docking with " 
                      << omp_get_num_threads() << " OpenMP threads." << std::endl;
        }
    }
    
    double t1 = omp_get_wtime();

    // Docking paralelo
    #pragma omp parallel for schedule(static)
    for (size_t idx = start; idx < end; ++idx) {
        size_t i = idx / ligands.size();
        size_t j = idx % ligands.size();
        localScores[idx - start] = performDocking(proteins[i], ligands[j]);
    }
    
    double t2 = omp_get_wtime();
    std::cout << "Process " << rank << " local execution time: " 
              << (t2 - t1)*1000 << " ms" << std::endl;

    // Recopilar resultados
    std::vector<float> scores;
    if (rank == 0) {
        scores.resize(total);
    }
    
    int localSize = static_cast<int>(localScores.size());
    std::vector<int> recvCounts(size);
    MPI_Gather(&localSize, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    std::vector<int> displs(size, 0);
    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + recvCounts[i - 1];
        }
    }
    
    MPI_Gatherv(localScores.data(), localSize, MPI_FLOAT,
                scores.data(), recvCounts.data(), displs.data(), MPI_FLOAT,
                0, MPI_COMM_WORLD);

    return scores;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    std::string proteinsDir;
    std::string ligandsDir;
    bool verbose;

    parseArguments(argc, argv, proteinsDir, ligandsDir, verbose);

    DataManager dataManager;
    std::vector<Molecule> proteins, ligands;

    if (!dataManager.loadProteins(proteinsDir, proteins)) {
        std::cerr << "Error loading proteins." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (!dataManager.loadLigands(ligandsDir, ligands)) {
        std::cerr << "Error loading ligands." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    double t1, t2;
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();

    std::vector<float> scores = hybrid_docking(proteins, ligands);

    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        std::cout << "Execution time: " << (t2 - t1) * 1000 << " ms" << std::endl;
    if (rank == 0 && verbose)
        analyzeDockingResults(scores, proteins.size(), ligands.size());

    MPI_Finalize();
    return EXIT_SUCCESS;
}