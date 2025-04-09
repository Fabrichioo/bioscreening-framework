# Guía de Instalación y Compilación

## Requisitos Previos

- Compilador con soporte para C++14 o superior.
- [CMake](https://cmake.org/) versión 3.10 o superior.
- Bibliotecas y herramientas para:
  - **OpenMP:** Normalmente integrado en compiladores modernos.
  - **MPI:** Por ejemplo, OpenMPI o MPICH.
  - **CUDA:** SDK de NVIDIA para compilación y ejecución en GPU.

## Instrucciones

1. **Clonar el repositorio:**
   ```
   git clone https://tu-repositorio-url.git
   cd BioScreeningFramework
   ```

2. **Crear una carpeta de compilación:**
   ```
   mkdir build && cd build
   ```

3. **Configurar el proyecto con CMake:**
   ```
   cmake ..
   ```
   Para un build en modo Release:
   ```
   cmake -DCMAKE_BUILD_TYPE=Release ..
   ```

4. **Compilar el proyecto:**
   ```
   make
   ```

El ejecutable `bioscreening` se generará en el directorio de compilación.