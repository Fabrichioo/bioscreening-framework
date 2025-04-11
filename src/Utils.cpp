#include "Utils.h"
#include "Docking.h"   // To use DockingResult
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
        std::cerr << "No docking results to analyze." << std::endl;
        return;
    }
    
    // Determine the pair with the optimal score (lower score = better affinity).
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
    
    std::cout << "\nBest docking:" << std::endl;
    std::cout << "  Protein " << bestProteinIndex 
              << " with Ligand " << bestLigandIndex 
              << " (Score: " << bestScore << ")" << std::endl;
    
    // Create a ranking of the results.
    std::vector<DockingResult> results;
    results.reserve(scores.size());
    for (int i = 0; i < numProteins; ++i) {
        for (int j = 0; j < numLigands; ++j) {
            int idx = i * numLigands + j;
            // Use explicit constructor via initializer list
            results.push_back(DockingResult{i, j, scores[idx]});
        }
    }
    
    // Sort the results: lower score is better.
    std::sort(results.begin(), results.end(), [](const DockingResult &a, const DockingResult &b) {
        return a.score < b.score;
    });
    
    // Show a ranking of the top 10 results.
    std::cout << "\nTop 10 docking results:" << std::endl;
    int top = std::min(10, static_cast<int>(results.size()));
    for (int i = 0; i < top; i++) {
        std::cout << (i + 1) << ") Protein " << results[i].proteinIndex
                  << ", Ligand " << results[i].ligandIndex
                  << ", Score: " << results[i].score << std::endl;
    }
}

void parseArguments(int argc, char* argv[], std::string &proteinsDir, std::string &ligandsDir, bool &verbose) {
    proteinsDir = DEFAULT_PROTEINS_DIR;
    ligandsDir  = DEFAULT_LIGANDS_DIR;
    verbose = false;
    
    int dirCount = 0;
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printHelp();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-v") {
            verbose = true;
        } else {
            if (dirCount == 0) {
                proteinsDir = arg;
                dirCount++;
            } else if (dirCount == 1) {
                ligandsDir = arg;
                dirCount++;
            }
        }
    }

    if(verbose){
        std::cout << "Current configuration:" << std::endl;
        std::cout << " Proteins path: " << proteinsDir << std::endl;
        std::cout << " Ligands path: " << ligandsDir << std::endl;
        std::cout << " Verbose mode: " << (verbose ? "enabled" : "disabled") << std::endl;
        std::cout << std::endl;
    }
}

void printHelp() {
    std::cout << "Usage: program [options] [proteins_path] [ligands_path]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -h, --help Displays this help and exits." << std::endl;
    std::cout << " -v Enables verbose mode and prints docking analysis." << std::endl;
    std::cout << "If no paths are specified, the following defaults will be used:" << std::endl;
    std::cout << " Proteins: " << DEFAULT_PROTEINS_DIR << std::endl;
    std::cout << " Ligands: " << DEFAULT_LIGANDS_DIR << std::endl;
}