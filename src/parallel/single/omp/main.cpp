#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "DataManager.h"
#include "Molecule.h"
#include "Docking.h"
#include "Utils.h"
#include <omp.h>

// void omp_docking(std::vector<Molecule> proteins, std::vector<Molecule> ligands, std::vector<float> scores){ 
// }

int main(int argc, char* argv[]) {

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

    Timer timer;
    timer.start();

    std::vector<float> scores;
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

    timer.stop();
    std::cout << "Tiempo de ejecución: " << timer.elapsedMilliseconds() << " ms" << std::endl;

    if(verbose)
        analyzeDockingResults(scores, proteins.size(), ligands.size());
    
    exit(EXIT_SUCCESS);
}