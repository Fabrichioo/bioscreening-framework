# Guía de Uso

## Ejecución del Framework

El ejecutable `bioscreening` admite distintos modos de ejecución:
- **sequential:** Ejecución secuencial.
- **openmp:** Paralelización en multi-hilo con OpenMP.
- **mpi:** Paralelización distribuida con MPI.
- **cuda:** Aceleración con CUDA en GPU.

### Ejemplo de Ejecución

- Versión Secuencial:
  ```
  ./bioscreening sequential
  ```

- Versión OpenMP (establece el número de hilos con la variable OMP_NUM_THREADS):
  ```
  export OMP_NUM_THREADS=4
  ./bioscreening openmp
  ```

- Versión MPI (ejecuta con 4 procesos):
  ```
  mpirun -np 4 ./bioscreening mpi
  ```

- Versión CUDA:
  ```
  ./bioscreening cuda
  ```

## Uso de los Scripts de Examples

Dentro del directorio `examples/` encontrarás scripts de shell para ejecutar cada versión:
- `run_sequential.sh`
- `run_openmp.sh`
- `run_mpi.sh`
- `run_cuda.sh`

Asegúrate de otorgar permisos de ejecución:
```
chmod +x examples/*.sh
```
Y ejecuta el script deseado:
```
./examples/run_openmp.sh
```

## Ejecución de Tests

Para verificar la funcionalidad de cada módulo, revisa la carpeta `tests/` y ejecuta los tests de manera individual. Se recomienda compilar con un comando similar a:
```
g++ -std=c++14 tests/test_Molecule.cpp src/Molecule.cpp -Iinclude -o test_Molecule && ./test_Molecule
```
