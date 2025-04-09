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

    // Cargar datos (rutas de ejemplo)
    dataManager.loadProteins("data/proteins/", proteins);
    dataManager.loadLigands("data/ligands/", ligands);

    Timer timer;
    timer.start();

    std::vector<float> scores;

    if (mode == "sequential") {
        // Versión secuencial: bucle anidado sobre todos los pares proteína-ligando.
        for (const auto& protein : proteins) {
            for (const auto& ligand : ligands) {
                float score = performDocking(protein, ligand);
                scores.push_back(score);
            }
        }
    } else if (mode == "openmp") {
        openmpDocking(proteins, ligands, scores);
    } else if (mode == "mpi") {
        mpiDocking(proteins, ligands, scores);
    } else if (mode == "cuda") {
        // cudaDocking(proteins, ligands, scores);
    } else {
        std::cerr << "Modo '" << mode << "' no reconocido." << std::endl;
        return 1;
    }

    timer.stop();
    std::cout << "Tiempo de ejecución: " << timer.elapsedMilliseconds() << " ms" << std::endl;

    // Ejemplo: se imprimen los primeros 10 scores
    for (size_t i = 0; i < std::min(scores.size(), size_t(10)); ++i) {
        std::cout << "Score[" << i << "] = " << scores[i] << std::endl;
    }

    return 0;
}