// src/Parallel/Hybrid_OpenMP_CUDA.cpp

#ifdef __CUDACC__

#include "Parallel.h"
#include "Docking.h"      // Para performDocking (versión CPU)
#include "Molecule.h"
#include <cuda_runtime.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>

// ---------------------------------------------------------------------
// Kernel CUDA que calcula el docking para un subrango global.
// Se reciben las mismas entradas que el kernel completo y dos parámetros
// 'start' y 'end' que definen el rango (en el índice global) a procesar.
__global__ void dockingKernelSubset(
    int numProteins,
    int numLigands,
    const int* proteinAtomCounts,    // Número de átomos por proteína
    const int* proteinAtomOffsets,   // Offset (índice inicial) de cada proteína en el array de átomos
    const float* proteinAtoms,       // Arreglo aplanado: [x, y, z, …] de todas las proteínas
    const int* ligandAtomCounts,     // Número de átomos por ligando
    const int* ligandAtomOffsets,    // Offsets en el arreglo de átomos de ligandos
    const float* ligandAtoms,        // Arreglo aplanado de coordenadas para ligandos
    float* scores,                   // Array de resultados locales (tamaño = gpu_count)
    int start,                       // Índice global de inicio (dentro del total de evaluaciones)
    int end                          // Índice global final (no inclusivo)
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;  // índice local en el subrango
    int globalIdx = idx + start;
    if (globalIdx >= end) return;

    // Mapear el índice global en (i, j)
    int i = globalIdx / numLigands;   // índice de proteína
    int j = globalIdx % numLigands;     // índice de ligando

    float score = 0.0f;
    int pCount = proteinAtomCounts[i];
    int pOffset = proteinAtomOffsets[i];
    int lCount = ligandAtomCounts[j];
    int lOffset = ligandAtomOffsets[j];

    for (int a = 0; a < pCount; a++) {
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
            float r = sqrtf(dx * dx + dy * dy + dz * dz);
            if (r < 1e-6f) r = 1e-6f;
            float inv_r = 1.0f / r;
            float inv_r6 = inv_r * inv_r * inv_r * inv_r * inv_r * inv_r;
            float inv_r12 = inv_r6 * inv_r6;
            score += 4.0f * (inv_r12 - inv_r6);
        }
    }
    // Guardar el resultado en la posición local (0 <= idx < gpu_count)
    scores[idx] = score;
}

// ---------------------------------------------------------------------
// Función híbrida que reparte las evaluaciones de docking entre
// la CPU (con OpenMP) y la GPU (con CUDA).
//
// Se divide el total de evaluaciones (numProteins * numLigands) en dos partes:
//   • La primera parte se procesa en la CPU usando OpenMP.
//   • La segunda parte se procesa en la GPU.
// Los resultados se combinan en el vector 'scores' (de tamaño total).
void hybridDockingOpenmpCuda(const std::vector<Molecule>& proteins,
                               const std::vector<Molecule>& ligands,
                               std::vector<float>& scores) {
    int numProteins = proteins.size();
    int numLigands = ligands.size();
    int total = numProteins * numLigands;
    if(total == 0) {
        scores.clear();
        return;
    }

    // Dividir el trabajo: por ejemplo, la mitad para CPU y la otra mitad para GPU.
    int cpu_count = total / 2;         // evaluaciones en CPU
    int gpu_count = total - cpu_count;   // evaluaciones en GPU

    // --- PARTE CPU: OpenMP ---
    std::vector<float> cpu_scores(cpu_count);
    #pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < cpu_count; ++idx) {
        int i = idx / numLigands;
        int j = idx % numLigands;
        // Se utiliza la función CPU de docking (por ejemplo, performDocking)
        cpu_scores[idx] = performDocking(proteins[i], ligands[j]);
    }

    // --- PARTE GPU: Preparar datos ---
    // Se requiere aplanar la información de las moléculas para transferirla a la GPU.
    std::vector<int> h_proteinAtomCounts(numProteins);
    std::vector<int> h_proteinAtomOffsets(numProteins);
    int totalProteinAtoms = 0;
    for (int i = 0; i < numProteins; i++) {
        h_proteinAtomCounts[i] = proteins[i].getAtoms().size();
        h_proteinAtomOffsets[i] = totalProteinAtoms;
        totalProteinAtoms += proteins[i].getAtoms().size();
    }
    std::vector<float> h_proteinAtoms(totalProteinAtoms * 3);
    int idxAtom = 0;
    for (int i = 0; i < numProteins; i++) {
        const std::vector<Atom>& atoms = proteins[i].getAtoms();
        for (size_t j = 0; j < atoms.size(); j++) {
            h_proteinAtoms[idxAtom++] = atoms[j].x;
            h_proteinAtoms[idxAtom++] = atoms[j].y;
            h_proteinAtoms[idxAtom++] = atoms[j].z;
        }
    }

    std::vector<int> h_ligandAtomCounts(numLigands);
    std::vector<int> h_ligandAtomOffsets(numLigands);
    int totalLigandAtoms = 0;
    for (int i = 0; i < numLigands; i++) {
        h_ligandAtomCounts[i] = ligands[i].getAtoms().size();
        h_ligandAtomOffsets[i] = totalLigandAtoms;
        totalLigandAtoms += ligands[i].getAtoms().size();
    }
    std::vector<float> h_ligandAtoms(totalLigandAtoms * 3);
    idxAtom = 0;
    for (int i = 0; i < numLigands; i++) {
        const std::vector<Atom>& atoms = ligands[i].getAtoms();
        for (size_t j = 0; j < atoms.size(); j++) {
            h_ligandAtoms[idxAtom++] = atoms[j].x;
            h_ligandAtoms[idxAtom++] = atoms[j].y;
            h_ligandAtoms[idxAtom++] = atoms[j].z;
        }
    }

    // --- PARTE GPU: Reservar y transferir memoria en el dispositivo ---
    int *d_proteinAtomCounts, *d_proteinAtomOffsets;
    float *d_proteinAtoms;
    int *d_ligandAtomCounts, *d_ligandAtomOffsets;
    float *d_ligandAtoms;
    float *d_gpu_scores;  // Array para almacenar los resultados GPU (tamaño = gpu_count)

    cudaMalloc((void**)&d_proteinAtomCounts, numProteins * sizeof(int));
    cudaMalloc((void**)&d_proteinAtomOffsets, numProteins * sizeof(int));
    cudaMalloc((void**)&d_proteinAtoms, totalProteinAtoms * 3 * sizeof(float));
    cudaMalloc((void**)&d_ligandAtomCounts, numLigands * sizeof(int));
    cudaMalloc((void**)&d_ligandAtomOffsets, numLigands * sizeof(int));
    cudaMalloc((void**)&d_ligandAtoms, totalLigandAtoms * 3 * sizeof(float));
    cudaMalloc((void**)&d_gpu_scores, gpu_count * sizeof(float));

    cudaMemcpy(d_proteinAtomCounts, h_proteinAtomCounts.data(), numProteins * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_proteinAtomOffsets, h_proteinAtomOffsets.data(), numProteins * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_proteinAtoms, h_proteinAtoms.data(), totalProteinAtoms * 3 * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(d_ligandAtomCounts, h_ligandAtomCounts.data(), numLigands * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ligandAtomOffsets, h_ligandAtomOffsets.data(), numLigands * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ligandAtoms, h_ligandAtoms.data(), totalLigandAtoms * 3 * sizeof(float), cudaMemcpyHostToDevice);

    // --- PARTE GPU: Lanzar el kernel para el subrango [cpu_count, total) ---
    int blockSize = 256;
    int gridSize = (gpu_count + blockSize - 1) / blockSize;
    int start = cpu_count;  // primer índice global a procesar en GPU
    int end = total;        // último índice (no inclusivo)
    dockingKernelSubset<<<gridSize, blockSize>>>(numProteins, numLigands,
         d_proteinAtomCounts, d_proteinAtomOffsets, d_proteinAtoms,
         d_ligandAtomCounts, d_ligandAtomOffsets, d_ligandAtoms,
         d_gpu_scores, start, end);
    cudaDeviceSynchronize();

    // --- PARTE GPU: Copiar resultados de la GPU al host ---
    std::vector<float> gpu_scores(gpu_count);
    cudaMemcpy(gpu_scores.data(), d_gpu_scores, gpu_count * sizeof(float), cudaMemcpyDeviceToHost);

    // Liberar memoria en el dispositivo.
    cudaFree(d_proteinAtomCounts);
    cudaFree(d_proteinAtomOffsets);
    cudaFree(d_proteinAtoms);
    cudaFree(d_ligandAtomCounts);
    cudaFree(d_ligandAtomOffsets);
    cudaFree(d_ligandAtoms);
    cudaFree(d_gpu_scores);

    // --- COMBINAR RESULTADOS ---
    scores.resize(total);
    // Inserta los resultados de la CPU en las posiciones [0, cpu_count)
    for (int i = 0; i < cpu_count; i++) {
        scores[i] = cpu_scores[i];
    }
    // Inserta los resultados de la GPU en las posiciones [cpu_count, total)
    for (int i = 0; i < gpu_count; i++) {
        scores[cpu_count + i] = gpu_scores[i];
    }
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