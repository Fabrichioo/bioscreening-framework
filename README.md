# BioScreeningFramework

BioScreeningFramework es un framework didáctico y modular para el cribado virtual de proteínas, orientado a la aplicación de técnicas de Computación de Alto Rendimiento (HPC) en biotecnología.  
Este proyecto forma parte de la asignatura de Programación de Arquitecturas Multinúcleos y permite comparar una solución secuencial con implementaciones utilizando OpenMP, MPI y CUDA.

---

## Características

- **Multimodalidad:** Soporta ejecución secuencial, OpenMP (paralelización en multi-hilo), MPI (paralelización distribuida) y CUDA (aceleración en GPU).
- **Modularidad:** Diseño orientado a objetos, dividiendo la funcionalidad en módulos (DataManager, Molecule, Docking, Parallel, Utils).
- **Extensibilidad:** La estructura del framework permite ampliar o reemplazar partes específicas (por ejemplo, añadir nuevos formatos de entrada o algoritmos de scoring).
- **Enfoque didáctico:** Documentación detallada y tests unitarios e integración para facilitar la comprensión de cada componente.
- **Herramientas de Build:** Uso de CMake para la generación de makefiles y la integración de dependencias como MPI, CUDA y OpenMP.

---

## Estructura del Proyecto

El repositorio se organiza de la siguiente manera:

- **docs/**  
  Documentación exhaustiva del framework. Puedes navegar por los documentos:
  - [Overview](docs/Overview.md): Visión general y objetivos del proyecto.
  - [Architecture](docs/Architecture.md): Descripción de la arquitectura, módulos y organización.
  - [Installation](docs/Installation.md): Guía de instalación, requisitos y compilación.
  - [Usage](docs/Usage.md): Instrucciones de uso y ejemplos de ejecución.
  - [DeveloperGuide](docs/DeveloperGuide.md): Guía para desarrolladores y colaboradores.

- **data/**  
  Archivos de datos mockup para proteínas y ligandos (formatos simplificados para la fase de prueba).

- **examples/**  
  Scripts de ejecución para cada versión (secuencial, OpenMP, MPI y CUDA).

- **include/**  
  Archivos de cabecera (headers) de los módulos del framework.

- **src/**  
  Código fuente del framework, incluídas las implementaciones de cada módulo y subcarpetas para las estrategias paralelas.

- **tests/**  
  Tests unitarios e de integración para validar cada componente de manera independiente.

- **CMakeLists.txt**  
  Archivo de configuración para construir el proyecto utilizando CMake.

- **README.md**  
  Este archivo que ofrece una visión global del proyecto y las instrucciones generales.

---

## Requisitos

- C++ (estándar C++14 o superior)
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

   *Si deseas configurar un build en modo Release, puedes usar:*
   ```
   cmake -DCMAKE_BUILD_TYPE=Release ..
   ```

4. Compila el proyecto:
   ```
   make
   ```

   El ejecutable se llamará `bioscreening`.

---

## Ejecución

Existen distintos scripts de ejecución en la carpeta `examples/` para ilustrar el uso de cada versión:

- **Versión Secuencial:**  
  Ejecutar:
  ```
  ./bioscreening sequential
  ```
  O bien, utilizar el script:
  ```
  cd examples
  ./run_sequential.sh
  ```

- **Versión OpenMP:**  
  Ejecutar:
  ```
  export OMP_NUM_THREADS=4
  ./bioscreening openmp
  ```
  O usar:
  ```
  cd examples
  ./run_openmp.sh
  ```

- **Versión MPI:**  
  Ejecutar (por ejemplo, con 4 procesos):
  ```
  mpirun -np 4 ./bioscreening mpi
  ```
  O usando:
  ```
  cd examples
  ./run_mpi.sh
  ```

- **Versión CUDA:**  
  Ejecutar:
  ```
  ./bioscreening cuda
  ```
  O usar:
  ```
  cd examples
  ./run_cuda.sh
  ```

---

## Ejecución de Tests

Dentro de la carpeta `tests/` se encuentran varios archivos fuente con tests unitarios para cada módulo y modalidad paralela. Puedes compilar y ejecutar estos tests individualmente mediante CMake (configurando targets para tests) o compilándolos manualmente.

Por ejemplo, para testear el DataManager:
```
g++ -std=c++14 tests/test_DataManager.cpp src/DataManager.cpp -Iinclude -o test_DataManager
./test_DataManager
```

---

## Documentación Adicional

La carpeta `doc/` contiene documentación exhaustiva sobre:
- Visión general y objetivos.
- Arquitectura y diseño del framework.
- Guía de instalación, uso y ejecución.
- Guía para desarrolladores y colaboradores.

Se recomienda revisar estos documentos para entender a fondo el proyecto y para orientarse en futuras extensiones o mejoras.

---

## Contribuciones

Se agradecen contribuciones que mejoren el framework. Si deseas colaborar:

- Haz un fork del repositorio.
- Crea una rama para tu feature o solución.
- Envía un pull request con los cambios propuestos.

---

## Licencia

Este proyecto está bajo la Licencia MIT.