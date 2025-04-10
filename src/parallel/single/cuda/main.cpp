#ifdef __CUDACC__
#include "Parallel.h"
#include "Docking.h"
#include "Molecule.h"
#include <cuda_runtime.h>
#include <cmath>
#include <iostream>
#include <vector>

// Función dispositivo para calcular el potencial de Lennard-Jones.
// Usamos ε = 1 y σ = 1 para simplificar.
__device__ inline float lennardJones(float r) {
    // Evitar división por cero
    if (r < 1e-6f) r = 1e-6f;
    float inv_r = 1.0f / r;
    float inv_r6 = inv_r * inv_r * inv_r * inv_r * inv_r * inv_r;
    float inv_r12 = inv_r6 * inv_r6;
    return 4.0f * (inv_r12 - inv_r6);
}

// Kernel: Cada hilo calcula el score de docking para una combinación proteína-ligando.
// Se reciben arrays con información “aplanada” de los átomos de las proteínas y ligandos.
__global__ void dockingKernel(
    int numProteins,
    int numLigands,
    const int* proteinAtomCounts,    // Número de átomos por proteína
    const int* proteinAtomOffsets,   // Offset (índice inicial) de cada proteína en el array de átomos
    const float* proteinAtoms,       // Arreglo aplanado: [x0,y0,z0, x1,y1,z1, ...] de todas las proteínas
    const int* ligandAtomCounts,     // Número de átomos por ligando
    const int* ligandAtomOffsets,    // Offsets en el arreglo de átomos de ligandos
    const float* ligandAtoms,        // Arreglo aplanado de coordenadas para ligandos
    float* scores                    // Array de scores: una posición por cada docking (proteína X ligando)
) {
    int dockingIdx = blockIdx.x * blockDim.x + threadIdx.x;
    int totalDockings = numProteins * numLigands;
    if (dockingIdx < totalDockings) {
        // Determinar el índice de proteína y ligando para este hilo
        int i = dockingIdx / numLigands;   // índice de proteína
        int j = dockingIdx % numLigands;     // índice de ligando

        int pCount = proteinAtomCounts[i];
        int pOffset = proteinAtomOffsets[i];
        int lCount = ligandAtomCounts[j];
        int lOffset = ligandAtomOffsets[j];

        float score = 0.0f;
        // Se suman las contribuciones de Lennard-Jones para cada par de átomos
        for (int a = 0; a < pCount; a++) {
            // Índice del átomo 'a' de la proteína (cada átomo tiene 3 coordenadas)
            int pAtomIdx = (pOffset + a) * 3;
            float px = proteinAtoms[pAtomIdx];
            float py = proteinAtoms[pAtomIdx + 1];
            float pz = proteinAtoms[pAtomIdx + 2];
            for (int b = 0; b < lCount; b++) {
                int lAtomIdx = (lOffset + b) * 3;
                float lx = ligandAtoms[lAtomIdx];
                float ly = ligandAtoms[lAtomIdx + 1];
                float lz = ligandAtoms[lAtomIdx + 2];
                float dx = px - lx;
                float dy = py - ly;
                float dz = pz - lz;
                float r = sqrtf(dx*dx + dy*dy + dz*dz);
                score += lennardJones(r);
            }
        }
        scores[dockingIdx] = score;
    }
}

// Función host que prepara los datos, lanza el kernel y recopila resultados.
// Se recorre el vector de Molecule y se generan arreglos “aplanados” para su transferencia.
void cudaDocking(const std::vector<Molecule>& proteins,
                 const std::vector<Molecule>& ligands,
                 std::vector<float>& scores) {
    std::cout << "Ejecutando docking con CUDA (potencial de Lennard-Jones)..." << std::endl;
    int numProteins = proteins.size();
    int numLigands = ligands.size();
    int totalDockings = numProteins * numLigands;
    scores.resize(totalDockings);

    // 1. Procesar proteínas: generar arrays con número de átomos, offsets y coordenadas
    std::vector<int> h_proteinAtomCounts(numProteins);
    std::vector<int> h_proteinAtomOffsets(numProteins);
    int totalProteinAtoms = 0;
    for (int i = 0; i < numProteins; i++) {
        h_proteinAtomCounts[i] = proteins[i].getAtoms().size();
        h_proteinAtomOffsets[i] = totalProteinAtoms;
        totalProteinAtoms += proteins[i].getAtoms().size();
    }
    std::vector<float> h_proteinAtoms(totalProteinAtoms * 3);
    int idx = 0;
    for (int i = 0; i < numProteins; i++) {
        const std::vector<Atom>& atoms = proteins[i].getAtoms();
        for (size_t j = 0; j < atoms.size(); j++) {
            h_proteinAtoms[idx++] = atoms[j].x;
            h_proteinAtoms[idx++] = atoms[j].y;
            h_proteinAtoms[idx++] = atoms[j].z;
        }
    }

    // 2. Procesar ligandos: extraer número de átomos, offsets y coordenadas
    std::vector<int> h_ligandAtomCounts(numLigands);
    std::vector<int> h_ligandAtomOffsets(numLigands);
    int totalLigandAtoms = 0;
    for (int i = 0; i < numLigands; i++) {
        h_ligandAtomCounts[i] = ligands[i].getAtoms().size();
        h_ligandAtomOffsets[i] = totalLigandAtoms;
        totalLigandAtoms += ligands[i].getAtoms().size();
    }
    std::vector<float> h_ligandAtoms(totalLigandAtoms * 3);
    idx = 0;
    for (int i = 0; i < numLigands; i++) {
        const std::vector<Atom>& atoms = ligands[i].getAtoms();
        for (size_t j = 0; j < atoms.size(); j++) {
            h_ligandAtoms[idx++] = atoms[j].x;
            h_ligandAtoms[idx++] = atoms[j].y;
            h_ligandAtoms[idx++] = atoms[j].z;
        }
    }

    // 3. Reservar memoria en el dispositivo para cada uno de los arreglos.
    int *d_proteinAtomCounts, *d_proteinAtomOffsets;
    float *d_proteinAtoms;
    int *d_ligandAtomCounts, *d_ligandAtomOffsets;
    float *d_ligandAtoms;
    float *d_scores;

    cudaMalloc((void**)&d_proteinAtomCounts, numProteins * sizeof(int));
    cudaMalloc((void**)&d_proteinAtomOffsets, numProteins * sizeof(int));
    cudaMalloc((void**)&d_proteinAtoms, totalProteinAtoms * 3 * sizeof(float));

    cudaMalloc((void**)&d_ligandAtomCounts, numLigands * sizeof(int));
    cudaMalloc((void**)&d_ligandAtomOffsets, numLigands * sizeof(int));
    cudaMalloc((void**)&d_ligandAtoms, totalLigandAtoms * 3 * sizeof(float));

    cudaMalloc((void**)&d_scores, totalDockings * sizeof(float));

    // 4. Transferir datos del host al dispositivo.
    cudaMemcpy(d_proteinAtomCounts, h_proteinAtomCounts.data(), numProteins * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_proteinAtomOffsets, h_proteinAtomOffsets.data(), numProteins * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_proteinAtoms, h_proteinAtoms.data(), totalProteinAtoms * 3 * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(d_ligandAtomCounts, h_ligandAtomCounts.data(), numLigands * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ligandAtomOffsets, h_ligandAtomOffsets.data(), numLigands * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ligandAtoms, h_ligandAtoms.data(), totalLigandAtoms * 3 * sizeof(float), cudaMemcpyHostToDevice);

    // 5. Lanzar el kernel:
    int blockSize = 256;
    int gridSize = (totalDockings + blockSize - 1) / blockSize;
    dockingKernel<<<gridSize, blockSize>>>(
        numProteins, numLigands,
        d_proteinAtomCounts, d_proteinAtomOffsets, d_proteinAtoms,
        d_ligandAtomCounts, d_ligandAtomOffsets, d_ligandAtoms,
        d_scores
    );
    cudaDeviceSynchronize();

    // 6. Copiar los resultados (scores) de vuelta al host.
    cudaMemcpy(scores.data(), d_scores, totalDockings * sizeof(float), cudaMemcpyDeviceToHost);

    // 7. Liberar la memoria reservada en el dispositivo.
    cudaFree(d_proteinAtomCounts);
    cudaFree(d_proteinAtomOffsets);
    cudaFree(d_proteinAtoms);
    cudaFree(d_ligandAtomCounts);
    cudaFree(d_ligandAtomOffsets);
    cudaFree(d_ligandAtoms);
    cudaFree(d_scores);
}

#endif // __CUDACC__


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
    std::cout << "Tiempo de ejecución: " << timer.elapsedMilliseconds() << " ms" << std::endl;

    if(verbose)
        analyzeDockingResults(scores, proteins.size(), ligands.size());
    
    exit(EXIT_SUCCESS);
}