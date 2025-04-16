#!/bin/bash
# Script para compilar la versión híbrida y ejecutar un barrido parametrizable de pruebas.
# Uso: ./run_test_hybrid.sh [iteraciones]
#   iteraciones -> (opcional) número de veces que se ejecutará cada prueba (default: 10)

# Definición de los arreglos a probar
omp_nums=(8 4 8)
mpi_nums=(1 4 8)

# Número de iteraciones por prueba (default 10)
ITERATIONS=10
if [ $# -ge 1 ]; then
    ITERATIONS="$1"
fi

echo "Se ejecutará cada prueba $ITERATIONS veces."
echo "Configuraciones de MPI: ${mpi_nums[*]}"
echo "Configuraciones de OMP: ${omp_nums[*]}"

# Compilamos la versión MPI (asumiendo que incluye OpenMP dentro)
./compile.sh 5

# Archivo de salida para resultados en formato CSV
csv_out="resultados_hibrido.csv"
echo "mpi_procs,omp_threads,average_ms,minimum_ms,maximum_ms" > "$csv_out"

# Bucle anidado para recorrer TODAS las combinaciones de MPI y OMP
for mpi_procs in "${mpi_nums[@]}"; do
    for omp_threads in "${omp_nums[@]}"; do
        echo "Ejecutando pruebas con MPI_PROCS=$mpi_procs, OMP_NUM_THREADS=$omp_threads ..."
        times=""
        for i in $(seq 1 $ITERATIONS); do
            output=$(OMP_NUM_THREADS=$omp_threads mpirun -np $mpi_procs -host localhost:$mpi_procs ./bioscreening)

            echo $output
            # Se espera que la salida contenga una línea como:
            # "Execution time: 1939.12 ms"
            time_ms=$(echo "$output" | grep -oE "Execution time: [0-9]+\.*[0-9]*" | awk '{print $3}')
            if [ -z "$time_ms" ]; then
                echo "No se pudo extraer el tiempo en la iteración $i (MPI=$mpi_procs OMP=$omp_threads)"
                continue
            fi
            times="$times $time_ms"
        done
        # Procesado estadístico de los resultados 
        if [ -n "$times" ]; then
            read avg min max <<< $(echo "$times" | LC_NUMERIC=C awk '{
                s=0; for(i=1; i<=NF; i++){
                    s+=$i;
                    if(i==1){min=$i; max=$i;} 
                    else {if($i<min) min=$i; if($i>max) max=$i;}
                } 
                avg=s/NF; printf "%.2f %.2f %.2f", avg,min,max;
            }')
            echo "$mpi_procs,$omp_threads,$avg,$min,$max" >> "$csv_out"
            echo "MPI $mpi_procs / OMP $omp_threads -> Promedio: $avg ms, Mínimo: $min ms, Máximo: $max ms"
        fi
    done
done

echo "Se han guardado los resultados en $csv_out"