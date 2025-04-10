#!/bin/bash

echo "Seleccione la versión que desea compilar:"
echo "1. Secuencial"
echo "2. OpenMP"
echo "3. MPI"

read -p "Ingrese el número de su elección: " choice

case $choice in
    1)
        echo "Compilando versión Secuencial..."
        g++ -std=c++14 -Iinclude -o bioscreening src/seq/*.cpp src/*.cpp
        ;;
    2)
        echo "Compilando versión OpenMP..."
        g++ -std=c++14 -Iinclude -o bioscreening src/parallel/single/omp/*.cpp src/*.cpp -lm -fopenmp -O3
        ;;
    3)
        echo "Compilando versión MPI..."
        mpic++ -std=c++14 -Iinclude -o bioscreening src/parallel/single/mpi/*.cpp src/*.cpp
        ;;
    *)
        echo "Opción no válida. Por favor, seleccione 1, 2 o 3."
        ;;
esac