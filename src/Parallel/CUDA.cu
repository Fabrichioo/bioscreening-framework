#ifdef __CUDA_ARCH__
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif

#include "Parallel.h"
#include "Docking.h"
#include <iostream>
#include <vector>
#include <cuda_runtime.h>

// Kernel de ejemplo: cada hilo asigna un score dummy
__global__ void dockingKernel(const int numProteins, const int numLigands, float* scores) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = numProteins * numLigands;
    if (idx < total) {
        // Cálculo dummy: en una implementación real se usarían los datos de las moléculas.
        scores[idx] = 42.0f;
    }
}

void cudaDocking(const std::vector<Molecule>& proteins,
                 const std::vector<Molecule>& ligands,
                 std::vector<float>& scores) {
    std::cout << "Ejecutando docking con CUDA..." << std::endl;
    size_t total = proteins.size() * ligands.size();
    scores.resize(total);

    float* d_scores = nullptr;
    cudaMalloc((void**)&d_scores, total * sizeof(float));

    int blockSize = 256;
    int gridSize = (total + blockSize - 1) / blockSize;
    dockingKernel<<<gridSize, blockSize>>>(proteins.size(), ligands.size(), d_scores);
    cudaDeviceSynchronize();

    cudaMemcpy(scores.data(), d_scores, total * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(d_scores);
}