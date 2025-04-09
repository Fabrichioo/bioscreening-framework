#!/bin/bash
# Script para ejecutar la versión OpenMP
# Ajusta OMP_NUM_THREADS según el número de hilos deseado
export OMP_NUM_THREADS=4
./bioscreening openmp