#include "Utils.h"
#include "Docking.h"   // Para usar DockingResult
#include <algorithm>
#include <vector>
#include <iostream>

Timer::Timer() {}

void Timer::start() {
    m_start = std::chrono::high_resolution_clock::now();
}

void Timer::stop() {
    m_end = std::chrono::high_resolution_clock::now();
}

double Timer::elapsedMilliseconds() const {
    return std::chrono::duration<double, std::milli>(m_end - m_start).count();
}

void logMessage(const std::string& msg) {
    std::cout << msg << std::endl;
}

void analyzeDockingResults(const std::vector<float>& scores, int numProteins, int numLigands) {
    if (scores.empty() || numProteins <= 0 || numLigands <= 0) {
        std::cerr << "No hay resultados de docking para analizar." << std::endl;
        return;
    }
    
    // Determinar el par con el score óptimo (score menor = mejor afinidad).
    int bestIndex = 0;
    float bestScore = scores[0];
    for (size_t idx = 1; idx < scores.size(); ++idx) {
        if (scores[idx] < bestScore) {
            bestScore = scores[idx];
            bestIndex = idx;
        }
    }
    int bestProteinIndex = bestIndex / numLigands;
    int bestLigandIndex = bestIndex % numLigands;
    
    std::cout << "\nMejor docking:" << std::endl;
    std::cout << "  Proteína " << bestProteinIndex 
              << " con Ligando " << bestLigandIndex 
              << " (Score: " << bestScore << ")" << std::endl;
    
    // Crear un ranking de los resultados.
    std::vector<DockingResult> results;
    results.reserve(scores.size());
    for (int i = 0; i < numProteins; ++i) {
        for (int j = 0; j < numLigands; ++j) {
            int idx = i * numLigands + j;
            // Usar constructor explícito mediante lista de inicializadores
            results.push_back(DockingResult{i, j, scores[idx]});
        }
    }
    
    // Ordenar los resultados: menor score es mejor.
    std::sort(results.begin(), results.end(), [](const DockingResult &a, const DockingResult &b) {
        return a.score < b.score;
    });
    
    // Mostrar un ranking de los 10 mejores resultados.
    std::cout << "\nTop 10 resultados de docking:" << std::endl;
    int top = std::min(10, static_cast<int>(results.size()));
    for (int i = 0; i < top; i++) {
        std::cout << (i + 1) << ") Proteína " << results[i].proteinIndex
                  << ", Ligando " << results[i].ligandIndex
                  << ", Score: " << results[i].score << std::endl;
    }
}