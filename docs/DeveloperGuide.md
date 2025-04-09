# Guía para Desarrolladores

Esta guía está pensada para colaboradores que deseen ampliar o mejorar el framework.

## Organización del Código

- **include/** y **src/**:  
  Se mantiene una separación estricta entre las declaraciones (headers) y las implementaciones.  
  Cada módulo (DataManager, Molecule, Docking, Parallel, Utils) está estructurado para facilitar la extensión y el mantenimiento.

- **src/Parallel/**:  
  Contiene implementaciones específicas para cada técnica de paralelización. Es posible agregar nuevos módulos siguiendo el mismo patrón.

## Estándares y Buenas Prácticas

- Utilizar C++14 y aprovechar las características modernas del lenguaje.
- Documentar cada función y clase con comentarios claros.
- Mantener el código modular y desacoplado.
- Ejecutar tests unitarios y de integración después de cada cambio importante.
- Seguir la estructura de commits y ramas en Git para facilitar la revisión del código.

## Añadiendo Nuevas Funcionalidades

- **Nuevos Formatos de Datos:**  
  Amplía el módulo DataManager implementando nuevas funciones de parseo en archivos, siguiendo la interfaz definida en `DataManager.h`.

- **Nuevos Métodos de Docking:**  
  Si se requiere un algoritmo de docking más avanzado, se sugiere crear una nueva clase o funciones adicionales en el módulo Docking, sin romper la interfaz actual.

- **Optimizaciones HPC:**  
  Explora la posibilidad de incorporar nuevas técnicas de paralelización o mejora de los kernels CUDA. Asegúrate de que las nuevas implementaciones respeten la interfaz definida en `Parallel.h`.

## Pruebas y Validación

Cada modificación debe ser acompañada por pruebas unitarias en la carpeta `tests/`.  
Se recomienda documentar los casos de prueba para facilitar su mantenimiento y validación en futuras versiones.
