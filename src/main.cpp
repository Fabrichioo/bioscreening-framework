#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "DataManager.h"
#include "Molecule.h"
#include "Docking.h"
#include "Parallel.h"
#include "Utils.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Uso: ./bioscreening [sequential|openmp|mpi|cuda]" << std::endl;
        return 1;
    }
    
    std::string mode = argv[1];
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

    std::vector<float> scores;
    if (mode == "sequential") {
        for (const auto &protein : proteins) {
            for (const auto &ligand : ligands) {
                float score = performDocking(protein, ligand);
                scores.push_back(score);
            }
        }
    } else if (mode == "openmp") {
        openmpDocking(proteins, ligands, scores);
    } else if (mode == "mpi") {
        // mpiDocking(proteins, ligands, scores);
    } else if (mode == "cuda") {
        // cudaDocking(proteins, ligands, scores);
    } else {
        std::cerr << "Modo '" << mode << "' no reconocido." << std::endl;
        return 1;
    }

    timer.stop();
    std::cout << "Tiempo de ejecución: " << timer.elapsedMilliseconds() << " ms" << std::endl;

    // Llamada a la función de Utils que analiza y muestra los resultados del docking.
    analyzeDockingResults(scores, proteins.size(), ligands.size());
    
    return 0;
}