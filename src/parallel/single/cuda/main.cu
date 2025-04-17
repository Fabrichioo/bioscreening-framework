#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cassert>
#include <cuda_runtime.h>

// Se incluyen las implementaciones existentes para manejo de datos y moléculas.
// Se asume que DataManager, Molecule, Atom, Docking y Utils se han implementado en CPU.
#include "DataManager.h"   // Funciones de carga de proteínas y ligandos
#include "Molecule.h"      // Definición de Molecule y Atom (struct con x, y, z y element)
#include "Docking.h"       // Versión secuencial (opcional para comparar)
#include "Utils.h"         // parseArguments, Timer, analyzeDockingResults, etc.

using namespace std;

//-----------------------------------------------------------------------------
// Estructura que usaremos en device para almacenar las coordenadas de un átomo.
// Se omite información adicional (por ejemplo, el elemento) que no se usa en el cálculo.
struct AtomGPU {
    float x, y, z;
};

//-----------------------------------------------------------------------------
// Kernel CUDA que ejecuta el docking entre cada proteína y ligando, sabiendo que
// cada molécula tiene un número constante de átomos.
// Los parámetros atomsPerProtein y atomsPerLigand son valores conocidos.
__global__
void dockingKernelConstant(const AtomGPU* proteinAtoms,
                           const AtomGPU* ligandAtoms,
                           float* scores,
                           int numProteins, int numLigands,
                           int atomsPerProtein, int atomsPerLigand)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int totalDockings = numProteins * numLigands;
    if (tid >= totalDockings)
        return;
    
    // A partir de tid se determina la proteína y el ligando correspondientes
    int proteinIdx = tid / numLigands;
    int ligandIdx  = tid % numLigands;
    
    float energy = 0.0f;
    const float epsilon = 1.0f;
    
    // Calcular los offsets en los arreglos aplanados, pues cada molécula es contigua
    int proteinOffset = proteinIdx * atomsPerProtein;
    int ligandOffset  = ligandIdx  * atomsPerLigand;
    
    // Cálculo del potencial de Lennard-Jones para cada par de átomos de la pareja (ligando, proteína)
    for (int i = 0; i < atomsPerLigand; i++) {
        AtomGPU aL = ligandAtoms[ligandOffset + i];
        for (int j = 0; j < atomsPerProtein; j++) {
            AtomGPU aP = proteinAtoms[proteinOffset + j];
            float dx = aL.x - aP.x;
            float dy = aL.y - aP.y;
            float dz = aL.z - aP.z;
            float r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < 1e-6f)
                continue;
            float r6  = r2 * r2 * r2;
            float r12 = r6 * r6;
            float ljEnergy = 4.0f * epsilon * ((1.0f / r12) - (1.0f / r6));
            energy += ljEnergy;
        }
    }
    
    scores[tid] = energy;
}

//-----------------------------------------------------------------------------
// Función principal donde se asume que el número de átomos por molécula es constante.
// Se define, por ejemplo, ATOMS_PER_PROTEIN = 2000 y ATOMS_PER_LIGAND = 30.
int main(int argc, char* argv[])
{
    // Definición de constantes para el número de átomos
    const int ATOMS_PER_PROTEIN = 1000; 
    const int ATOMS_PER_LIGAND  = 100;
    
    // Variables de entrada (directorios, modo verbose, etc.)
    string proteinsDir, ligandsDir;
    bool verbose = false;
    
    parseArguments(argc, argv, proteinsDir, ligandsDir, verbose);
    
    // Carga de moléculas usando DataManager (implementado en CPU)
    DataManager dataManager;
    vector<Molecule> proteins;
    vector<Molecule> ligands;
    
    if (!dataManager.loadProteins(proteinsDir, proteins)) {
        cerr << "Error loading proteins." << endl;
        exit(EXIT_FAILURE);
    }
    if (!dataManager.loadLigands(ligandsDir, ligands)) {
        cerr << "Error loading ligands." << endl;
        exit(EXIT_FAILURE);
    }
    
    int numProteins = proteins.size();
    int numLigands  = ligands.size();
    int totalDockings = numProteins * numLigands;
    
    // Se valida que cada proteína y ligando tenga el número de átomos esperado:
    for (const auto &protein : proteins) {
         assert(protein.getAtoms().size() == ATOMS_PER_PROTEIN);
    }
    for (const auto &ligand : ligands) {
         assert(ligand.getAtoms().size() == ATOMS_PER_LIGAND);
    }
    
    // Se construyen arreglos aplanados de átomos para proteínas y ligandos.
    // Cada proteína ocupará un bloque de ATOMS_PER_PROTEIN y similarmente para ligandos.
    vector<AtomGPU> flatProteinAtoms(numProteins * ATOMS_PER_PROTEIN);
    vector<AtomGPU> flatLigandAtoms(numLigands * ATOMS_PER_LIGAND);
    
    // Transferir datos de proteínas a un arreglo lineal
    for (int i = 0; i < numProteins; i++) {
         const vector<Atom>& atoms = proteins[i].getAtoms();
         for (int j = 0; j < ATOMS_PER_PROTEIN; j++) {
             flatProteinAtoms[i * ATOMS_PER_PROTEIN + j].x = atoms[j].x;
             flatProteinAtoms[i * ATOMS_PER_PROTEIN + j].y = atoms[j].y;
             flatProteinAtoms[i * ATOMS_PER_PROTEIN + j].z = atoms[j].z;
         }
    }
    
    // Transferir datos de ligandos a un arreglo lineal
    for (int i = 0; i < numLigands; i++) {
         const vector<Atom>& atoms = ligands[i].getAtoms();
         for (int j = 0; j < ATOMS_PER_LIGAND; j++) {
             flatLigandAtoms[i * ATOMS_PER_LIGAND + j].x = atoms[j].x;
             flatLigandAtoms[i * ATOMS_PER_LIGAND + j].y = atoms[j].y;
             flatLigandAtoms[i * ATOMS_PER_LIGAND + j].z = atoms[j].z;
         }
    }
    
    // Reservar memoria en la GPU para los átomos aplanados y para el vector de scores
    AtomGPU *d_proteinAtoms, *d_ligandAtoms;
    float* d_scores;
    size_t sizeProteinAtoms = flatProteinAtoms.size() * sizeof(AtomGPU);
    size_t sizeLigandAtoms  = flatLigandAtoms.size() * sizeof(AtomGPU);
    size_t sizeScores       = totalDockings * sizeof(float);
    
    cudaError_t err;
    err = cudaMalloc((void**)&d_proteinAtoms, sizeProteinAtoms);
    if (err != cudaSuccess) {
         cerr << "Error allocando memoria para d_proteinAtoms: " << cudaGetErrorString(err) << endl;
         exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void**)&d_ligandAtoms, sizeLigandAtoms);
    if (err != cudaSuccess) {
         cerr << "Error allocando memoria para d_ligandAtoms: " << cudaGetErrorString(err) << endl;
         exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void**)&d_scores, sizeScores);
    if (err != cudaSuccess) {
         cerr << "Error allocando memoria para d_scores: " << cudaGetErrorString(err) << endl;
         exit(EXIT_FAILURE);
    }
    
    // Transferir datos desde host a device
    cudaMemcpy(d_proteinAtoms, flatProteinAtoms.data(), sizeProteinAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ligandAtoms, flatLigandAtoms.data(), sizeLigandAtoms, cudaMemcpyHostToDevice);
    
    // Configurar los parámetros para el lanzamiento del kernel: cada hilo calcula una combinación proteína-ligando
    int threadsPerBlock = 256;
    int blocksPerGrid = (totalDockings + threadsPerBlock - 1) / threadsPerBlock;
    
    cout << "Modo CUDA con átomos constantes por molécula" << endl;
    
    // Medición del tiempo de ejecución en la GPU usando eventos
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    cudaEventRecord(start);
    dockingKernelConstant<<<blocksPerGrid, threadsPerBlock>>>(d_proteinAtoms,
                                                              d_ligandAtoms,
                                                              d_scores,
                                                              numProteins,
                                                              numLigands,
                                                              ATOMS_PER_PROTEIN,
                                                              ATOMS_PER_LIGAND);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    cout << "Tiempo de ejecución (kernel GPU): " << milliseconds << " ms" << endl;
    
    // Copiar resultados (scores) desde la GPU al host
    vector<float> scores(totalDockings);
    cudaMemcpy(scores.data(), d_scores, sizeScores, cudaMemcpyDeviceToHost);
    
    // (Opcional) Análisis de resultados
    if (verbose)
         analyzeDockingResults(scores, numProteins, numLigands);
    
    // Liberar memoria en la GPU y destruir eventos
    cudaFree(d_proteinAtoms);
    cudaFree(d_ligandAtoms);
    cudaFree(d_scores);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    return 0;
}
