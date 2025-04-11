#!/bin/bash
# Script para compilar la versión OpenMP y ejecutar pruebas parametrizables.
# Uso: ./run_test.sh [iteraciones]
#   iteraciones -> (opcional) número de veces que se ejecutará cada prueba (default: 10)

# Número de iteraciones por prueba (default 10)
ITERATIONS=10
if [ $# -ge 1 ]; then
    # Se asume que el primer parámetro es el número de iteraciones
    ITERATIONS="$1"
fi

echo "Se ejecutará cada prueba $ITERATIONS veces."

# Compilamos la versión OpenMP
./compile.sh 2

# Definición del arreglo de hilos a probar
omp_nums=(1 4 8)

# Archivo de salida para resultados en formato CSV
csv_out="resultados.csv"

# Escribimos la cabecera en el CSV
echo "threads,average_ms,minimum_ms,maximum_ms" > "$csv_out"

# Iteramos para cada cantidad de hilos definida en omp_nums
for threads in "${omp_nums[@]}"; do
    echo "Ejecutando pruebas con OMP_NUM_THREADS=$threads ..."
    # Variable para almacenar los tiempos (en milisegundos)
    times=""

    # Ejecutamos la prueba ITERATIONS veces
    for i in $(seq 1 $ITERATIONS); do
        # Exportamos OMP_NUM_THREADS y ejecutamos la aplicación
        output=$(OMP_NUM_THREADS=$threads ./bioscreening)
        
        # Se espera que la salida contenga una línea similar a:
        # "Execution time: 1939.12 ms"
        # Extraemos el valor numérico usando grep y awk
        time_ms=$(echo "$output" | grep -oE "Execution time: [0-9]+\.*[0-9]*" | awk '{print $3}')
        
        # Verificamos que se haya extraído el tiempo, de lo contrario lanzamos una advertencia
        if [ -z "$time_ms" ]; then
            echo "No se pudo extraer el tiempo en la iteración $i para $threads hilos."
            continue
        fi

        # Acumulamos el tiempo en la variable times
        times="$times $time_ms"
    done

    # Usamos awk para calcular el promedio, el mínimo y el máximo de los tiempos obtenidos
    read avg min max <<< $(echo "$times" | awk '{
        s=0;
        for(i=1; i<=NF; i++){
            s+=$i;
            if(i==1) {
                min=$i;
                max=$i;
            } else {
                if($i<min) min=$i;
                if($i>max) max=$i;
            }
        }
        avg=s/NF;
        # Imprimimos con dos decimales
        printf "%.2f %.2f %.2f", avg, min, max;
    }')

    # Guardamos el resultado en el archivo CSV
    echo "$threads,$avg,$min,$max" >> "$csv_out"
    echo "Con $threads hilos --> Promedio: $avg ms, Mínimo: $min ms, Máximo: $max ms"
done

echo "Se han guardado los resultados en $csv_out"