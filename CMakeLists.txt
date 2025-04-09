cmake_minimum_required(VERSION 3.10)
project(BioScreeningFramework LANGUAGES CXX CUDA)

# Establece el estándar de C++ y activa el uso de CUDA
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Incluir directorios de cabeceras
include_directories(${PROJECT_SOURCE_DIR}/include)

# Buscar y configurar OpenMP, si está disponible
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP encontrado")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(STATUS "OpenMP no encontrado")
endif()

# Buscar y configurar MPI (requerido)
find_package(MPI REQUIRED)
if(MPI_FOUND)
    message(STATUS "MPI encontrado")
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
endif()

# Buscar y configurar CUDA
find_package(CUDA)
if(CUDA_FOUND)
    message(STATUS "CUDA encontrado")
    enable_language(CUDA)
else()
    message(STATUS "CUDA no encontrado, los componentes CUDA no se compilarán")
endif()

# Lista de archivos fuente
set(SRC_FILES
    src/main.cpp
    src/DataManager.cpp
    src/Molecule.cpp
    src/Docking.cpp
    src/Utils.cpp
    src/Parallel/OpenMP.cpp
    src/Parallel/MPI.cpp
)

if(CUDA_FOUND)
    list(APPEND SRC_FILES src/Parallel/CUDA.cu)
endif()

# Crear el ejecutable
add_executable(bioscreening ${SRC_FILES})

# Vincular bibliotecas MPI
if(MPI_FOUND)
    target_link_libraries(bioscreening MPI::MPI_CXX)
endif()

# Vincular OpenMP si se encontró
if(OpenMP_CXX_FOUND)
    target_link_libraries(bioscreening OpenMP::OpenMP_CXX)
endif()

# Vincular CUDA si está disponible
if(CUDA_FOUND)
    target_link_libraries(bioscreening ${CUDA_LIBRARIES})
endif()

# Opcional: se pueden agregar tests u otros targets