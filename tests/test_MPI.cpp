#include "Parallel.h"
#include "Molecule.h"
#include "Docking.h"
#include <mpi.h>
#include <vector>
#include <cassert>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    std::vector<Molecule> proteins(1);
    std::vector<Molecule> ligands(1);
    
    proteins[0].addAtom({0, 0, 0, "C"});
    ligands[0].addAtom({1, 0, 0, "H"});

    std::vector<float> scores;
    mpiDocking(proteins, ligands, scores);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        // Con 1 proteína y 1 ligando, se espera 1 score = 100
        assert(scores.size() == 1);
        assert(scores[0] == 100.0f);
        std::cout << "test_MPI passed." << std::endl;
    }

    MPI_Finalize();
    return 0;
}