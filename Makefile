# Makefile para compilar la versión secuencial de BioScreeningFramework

# Compilador y flags
CXX       = g++
CXXFLAGS  = -std=c++14 -O2 -Wall -Iinclude -lm -fopenmp

# Lista de archivos fuente a compilar (solo la parte secuencial)
SRC       = src/main.cpp src/DataManager.cpp src/Molecule.cpp src/Docking.cpp src/Utils.cpp src/Parallel/OpenMP.cpp
OBJ       = $(SRC:.cpp=.o)

# Nombre del ejecutable
EXEC      = bioscreening_seq

# Regla principal: compilar el ejecutable
all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXEC)

# Regla para compilar cada archivo .cpp a .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Objetivo para limpiar archivos objeto y el ejecutable compilado
clean:
	rm -f $(OBJ) $(EXEC)

.PHONY: all clean