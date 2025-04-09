# BioScreeningFramework

BioScreeningFramework es un framework didáctico y modular para el cribado virtual de proteínas, orientado a la aplicación de técnicas de Computación de Alto Rendimiento (HPC) en biotecnología.  
Este proyecto forma parte de la asignatura de Programación de Arquitecturas Multinúcleos y permite comparar una solución secuencial con implementaciones utilizando OpenMP, MPI y CUDA.

---

## Características

- **Multimodalidad:** Soporta ejecución secuencial, OpenMP (paralelización en multi-hilo), MPI (paralelización distribuida) y CUDA (aceleración en GPU).
- **Modularidad:** Diseño orientado a objetos, con separación en módulos (DataManager, Molecule, Docking, Parallel, Utils) que facilitan su extensión o modificación.
- **Enfoque didáctico:** Incluye documentación, tests unitarios y ejemplos de ejecución para comprender el funcionamiento y la integración de técnicas HPC.
- **Generación de dataset intensivo:** Se incluye un script en Python (generate_dataset.py) que permite generar un conjunto de datos grande (proteínas y ligandos) configurable mediante parámetros y constantes, lo que facilita pruebas de rendimiento en escenarios intensivos.
- **Herramientas de Build:** Uso de CMake para la generación de makefiles y la integración con MPI, CUDA y OpenMP.

---

## Estructura del Proyecto

- **docs/**  
  Documentación que abarca:
  - Overview: Visión general y objetivos.
  - Architecture: Descripción de la arquitectura y organización de módulos.
  - Installation: Guía de instalación y compilación.
  - Usage: Instrucciones y ejemplos de ejecución.
  - DeveloperGuide: Guía para desarrolladores y colaboradores.

- **data/**  
  Directorio que contiene archivos de prueba y datos reales o generados para el cribado (formatos PDB para proteínas y SDF para ligandos).

- **examples/**  
  Scripts de ejecución para cada versión disponible:
  - run_sequential.sh
  - run_openmp.sh
  - run_mpi.sh
  - run_cuda.sh

- **include/**  
  Archivos de cabecera (headers) en C/C++ con la API del framework.

- **src/**  
  Código fuente del framework, con implementaciones de módulos y subdirectorios para estrategias paralelas (OpenMP, MPI, CUDA).

- **tests/**  
  Tests unitarios e integración para validar cada módulo y la correcta paralelización (se pueden compilar individualmente).

- **generate_dataset.py**  
  Script en Python que permite generar un conjunto de datos (dataset) intensivo, configurable mediante parámetros (número de proteínas, ligandos, átomos por molécula y directorios) y utilizando constantes para definir los valores por defecto.

- **CMakeLists.txt**  
  Archivo de configuración para compilar el proyecto usando CMake.

- **README.md**  
  Archivo actual que ofrece una visión global del proyecto y las instrucciones generales.

---

## Requisitos

- C++ (compatible con C++14 o superior; algunos módulos usan exclusivamente C++14 en este ejemplo)
- [CMake](https://cmake.org/) 3.10 o superior
- [MPI](https://www.mpi-forum.org/) (por ejemplo, OpenMPI o MPICH)
- [CUDA](https://developer.nvidia.com/cuda-zone) (si se desea compilar la versión CUDA)
- Soporte de OpenMP (generalmente incluido en compiladores modernos como GCC o Clang)

---

## Instrucciones de Instalación y Compilación

1. Clona el repositorio:
   ```
   git clone https://tu-repositorio-url.git
   cd BioScreeningFramework
   ```

2. Crea y accede a un directorio para la compilación:
   ```
   mkdir build && cd build
   ```

3. Configura el proyecto usando CMake:
   ```
   cmake ..
   ```
   *Para un build en modo Release:*
   ```
   cmake -DCMAKE_BUILD_TYPE=Release ..
   ```

4. Compila el proyecto:
   ```
   make
   ```
   El ejecutable se llamará `bioscreening`.

---

## Generación del Dataset Intensivo

Para simular entornos con cálculos intensivos, se incluye el script `generate_dataset.py` que permite generar un gran número de archivos de proteínas (en formato PDB) y ligandos (en formato SDF), de forma configurable por parámetros.

### Características del Script
- **Parámetros Opcionales:** Se pueden especificar el número de proteínas, ligandos, y la cantidad de átomos por cada uno, así como los directorios de salida.
- **Constantes por Defecto:** Todas las opciones tienen un valor por defecto definido mediante constantes al inicio del script, lo que evita literales dispersos y facilita modificaciones futuras.
- **Función de Ayuda:** Al ejecutar el script con la opción `--help` o `-h`, se mostrará un mensaje que explica el uso de cada parámetro.

### Ejemplos de Uso
- Mostrar la ayuda:
  ```
  python3 generate_dataset.py --help
  ```

- Generar el dataset usando los valores por defecto:
  ```
  python3 generate_dataset.py
  ```

- Generar el dataset especificando parámetros personalizados (por ejemplo, 200 proteínas, 150 ligandos, 1000 átomos por proteína y 75 átomos por ligando):
  ```
  python3 generate_dataset.py --num_proteins 200 --num_ligands 150 --num_atoms_protein 1000 --num_atoms_ligand 75
  ```

El script generará los archivos en los directorios configurados (por defecto, `data/proteins` y `data/ligands`).

---

## Ejecución

Se pueden ejecutar las diferentes versiones del framework usando los scripts en la carpeta `examples/` o directamente mediante línea de comandos. Por ejemplo:

- **Versión Secuencial:**  
  ```
  ./bioscreening sequential
  ```
  O:
  ```
  cd examples
  ./run_sequential.sh
  ```

- **Versión OpenMP:**  
  ```
  export OMP_NUM_THREADS=4
  ./bioscreening openmp
  ```
  O:
  ```
  cd examples
  ./run_openmp.sh
  ```

- **Versión MPI:**  
  Ejecución con 4 procesos:
  ```
  mpirun -np 4 ./bioscreening mpi
  ```
  O:
  ```
  cd examples
  ./run_mpi.sh
  ```

- **Versión CUDA:**  
  ```
  ./bioscreening cuda
  ```
  O:
  ```
  cd examples
  ./run_cuda.sh
  ```

---

## Ejecución de Tests

Dentro de la carpeta `tests/` se encuentran varios archivos fuente con tests unitarios para cada módulo y modalidad paralela. Puedes compilar y ejecutar estos tests individualmente mediante CMake o compilándolos manualmente.

Por ejemplo, para testear el DataManager:
```
g++ -std=c++14 tests/test_DataManager.cpp src/DataManager.cpp -Iinclude -o test_DataManager
./test_DataManager
```

---

## Documentación Adicional

La carpeta `docs/` contiene documentación adicional sobre:
- Visión general y objetivos del proyecto.
- Arquitectura y diseño del framework.
- Guía de instalación, uso y ejecución.
- Guía para desarrolladores y colaboradores.

Se recomienda revisar dichos documentos para comprender a fondo el funcionamiento del framework y poder extenderlo o modificarlo según las necesidades.

---

## Contribuciones

Las contribuciones y sugerencias son bienvenidas. Para colaborar:
- Realiza un fork del repositorio.
- Crea una rama con las modificaciones necesarias.
- Envía un pull request detallando las mejoras implementadas.

---

## Licencia

Este proyecto se distribuye bajo la Licencia MIT.