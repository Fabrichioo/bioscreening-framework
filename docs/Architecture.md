# Arquitectura del Framework

El framework se compone de varios módulos interconectados:

## Módulos Principales

1. **DataManager:**  
   - Responsable de la carga y el parseo de archivos que contienen información de moléculas.
   - Permite la incorporación futura de formatos de archivo (por ejemplo, PDB, MOL2).

2. **Molecule:**  
   - Define la estructura para representar una molécula a través de átomos y sus propiedades.
   - Provee métodos para agregar y consultar átomos.

3. **Docking:**  
   - Implementa un algoritmo simplificado para simular el acoplamiento (docking) entre una proteína y un ligando.
   - Calcula un score de unión basado en parámetros dummy, sirviendo como base para futuras mejoras.

4. **Parallel:**  
   - Abstrae las distintas técnicas de paralelización.
   - Contiene submódulos (OpenMP, MPI, CUDA) que implementan la paralelización de la evaluación de docking.
  
5. **Utils:**  
   - Funciones auxiliares para tareas comunes como logging, temporización y manejo de errores.

## Organización de Directorios

- **include/**: Archivos de cabecera con las declaraciones de clases y funciones.
- **src/**: Implementaciones de los módulos y el programa principal.
- **tests/**: Tests unitarios y de integración para validar el funcionamiento de cada componente.
- **examples/**: Scripts de ejecución para cada modalidad de ejecución.
- **data/**: Archivos mockup para representar proteínas y ligandos.
- **doc/**: Documentación detallada (este documento, guías y especificaciones).

La modularidad y separación de responsabilidades permiten mejorar o extender componentes individualmente sin alterar todo el sistema.