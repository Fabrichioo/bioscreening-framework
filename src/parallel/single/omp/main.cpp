#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "DataManager.h"
#include "Molecule.h"
#include "Docking.h"
#include "Utils.h"
#include <omp.h>

std::vector<float> omp_docking(const std::vector<Molecule>& proteins, 
                               const std::vector<Molecule>& ligands) {
    size_t total = proteins.size() * ligands.size();
    std::vector<float> scores(total);
    
    double t1 = omp_get_wtime();

    #pragma omp parallel 
    {
        #pragma omp master
        {
            std::cout << "Running docking with OpenMP with " 
                    << omp_get_num_threads() << " active threads..." << std::endl;
        }
        
        #pragma omp for collapse(2)
        for (size_t i = 0; i < proteins.size(); ++i) {
            for (size_t j = 0; j < ligands.size(); ++j) {
                size_t idx = i * ligands.size() + j;
                scores[idx] = performDocking(proteins[i], ligands[j]);
            }
        }
    }
    double t2 = omp_get_wtime();
    std::cout << "Execution time: " << t2 - t1 << " s" << std::endl;

    return scores;
}


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

    std::vector<float> scores = omp_docking(proteins, ligands);

    if(verbose)
        analyzeDockingResults(scores, proteins.size(), ligands.size());
    
    exit(EXIT_SUCCESS);
}