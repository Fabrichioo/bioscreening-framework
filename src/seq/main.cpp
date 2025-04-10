#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "DataManager.h"
#include "Molecule.h"
#include "Docking.h"
#include "Utils.h"

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