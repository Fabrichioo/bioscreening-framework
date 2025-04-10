#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "DataManager.h"
#include "Molecule.h"
#include "Docking.h"
#include "Utils.h"

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

    std::cout << "Sequential Mode" << std::endl;
    std::vector<float> scores;
    for (const auto &protein : proteins) {
        for (const auto &ligand : ligands) {
            float score = performDocking(protein, ligand);
            scores.push_back(score);
        }
    }

    timer.stop();
    std::cout << "Execution time: " << timer.elapsedMilliseconds() << " ms" << std::endl;

    if(verbose)
        analyzeDockingResults(scores, proteins.size(), ligands.size());
    
    exit(EXIT_SUCCESS);
}