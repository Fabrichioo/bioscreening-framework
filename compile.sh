#!/bin/bash

# Función para mostrar el menú y leer la opción desde la entrada estándar
mostrar_menu() {
    echo "Seleccione la versión que desea compilar:"
    echo "1. Secuencial"
    echo "2. OpenMP"
    echo "3. MPI"
    echo "4. CUDA"
    echo "5. OpenMP + MPI"
    read -p "Ingrese el número de su elección: " choice
}

# Verifica si se pasó un argumento y si es válido (1, 2, 3 o 4)
if [ $# -eq 0 ] || [[ "$1" != "1" && "$1" != "2" && "$1" != "3" && "$1" != "4" && "$1" != "5" ]]; then
    # Si no hay parámetro o es inválido, se muestra el menú interactivo
    mostrar_menu
else
    # Si se pasó un parámetro válido, se toma ese valor y se omite el menú
    choice="$1"
fi

# Ejecutar la acción correspondiente según la opción elegida
case $choice in
    1)
        echo "Compilando versión Secuencial..."
        g++ -std=c++14 -Iinclude -o bioscreening src/seq/*.cpp src/*.cpp -O3
        ;;
    2)
        echo "Compilando versión OpenMP..."
        g++ -std=c++14 -Iinclude -o bioscreening src/parallel/single/omp/*.cpp src/*.cpp -lm -fopenmp -O3
        ;;
    3)
        echo "Compilando versión MPI..."
        mpic++ -std=c++14 -Iinclude -o bioscreening src/parallel/single/mpi/*.cpp src/*.cpp -O3
        ;;
    4)
        echo "Compilando versión CUDA..."
        /usr/local/cuda-12/bin/nvcc -std=c++14 -Iinclude -o bioscreening src/*.cpp src/parallel/single/cuda/*
        ;;
    5)
        echo "Compilando versión OpenMP + MPI..."
        mpic++ -std=c++14 -Iinclude -o bioscreening src/parallel/hibrid/omp_mpi/*.cpp src/*.cpp -O3 -lm -fopenmp
        ;;
    *)
        echo "Opción no válida. Por favor, seleccione 1, 2, 3 o 4."
        ;;
esac