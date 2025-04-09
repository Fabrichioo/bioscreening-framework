/* src/DataManager.cpp */
#include "DataManager.h"
#include "Molecule.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <dirent.h>      // Para la iteración de directorios
#include <sys/stat.h>    // Para stat() y verificar tipos de archivo
#include <cstring>       // Para strcmp

// Función auxiliar para parsear archivos PDB
static Molecule parsePDB(const std::string& filename) {
    Molecule mol;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error abriendo archivo: " << filename << std::endl;
        return mol;
    }
    std::string line;
    while (std::getline(file, line)) {
        // Verificar que la línea comience con "ATOM" o "HETATM"
        if (line.size() >= 4 && (line.substr(0, 4) == "ATOM" || 
           (line.size() >= 6 && line.substr(0, 6) == "HETATM"))) {
            try {
                // Extraer coordenadas (columnas fijas según el estándar PDB)
                std::string xStr = line.substr(30, 8);
                std::string yStr = line.substr(38, 8);
                std::string zStr = line.substr(46, 8);
                float x = std::stof(xStr);
                float y = std::stof(yStr);
                float z = std::stof(zStr);

                // Extraer el símbolo del elemento (col. 77-78)
                std::string element;
                if (line.size() >= 78) {
                    element = line.substr(76, 2);
                    element.erase(std::remove_if(element.begin(), element.end(), ::isspace), element.end());
                } else {
                    element = "X";  // Valor por defecto
                }
                Atom atom = { x, y, z, element };
                mol.addAtom(atom);
            } catch (const std::exception& e) {
                std::cerr << "Error parseando línea en " << filename << ": " << e.what() << std::endl;
            }
        }
    }
    file.close();
    return mol;
}

// Función auxiliar para parsear archivos SDF
static Molecule parseSDF(const std::string& filename) {
    Molecule mol;
    std::ifstream file(filename);
    if (!file.is_open()){
        std::cerr << "Error abriendo archivo: " << filename << std::endl;
        return mol;
    }
    std::string line;
    // Omitir las primeras 3 líneas de cabecera
    for (int i = 0; i < 3 && std::getline(file, line); ++i) {}

    // La cuarta línea es la línea de conteo
    if (!std::getline(file, line)) {
        std::cerr << "Error leyendo línea de conteo en " << filename << std::endl;
        return mol;
    }
    int numAtoms = 0;
    try {
        numAtoms = std::stoi(line.substr(0, 3));
    } catch (const std::exception& e) {
        std::cerr << "Error parseando número de átomos en " << filename << ": " << e.what() << std::endl;
        return mol;
    }
    // Parsear las siguientes numAtoms líneas (sección de átomos)
    for (int i = 0; i < numAtoms; ++i) {
        if (!std::getline(file, line))
            break;
        try {
            // Para SDF (V2000):  
            // X: columnas 1–10, Y: columnas 11–20, Z: columnas 21–30, Elemento: columnas 32–34
            std::string xStr = line.substr(0, 10);
            std::string yStr = line.substr(10, 10);
            std::string zStr = line.substr(20, 10);
            float x = std::stof(xStr);
            float y = std::stof(yStr);
            float z = std::stof(zStr);

            std::string element;
            if (line.size() >= 34) {
                element = line.substr(31, 3);
                element.erase(std::remove_if(element.begin(), element.end(), ::isspace), element.end());
            } else {
                element = "X";
            }
            Atom atom = { x, y, z, element };
            mol.addAtom(atom);
        } catch (const std::exception & e) {
            std::cerr << "Error parseando átomo en " << filename << ": " << e.what() << std::endl;
        }
    }
    file.close();
    return mol;
}

// Constructor y Destructor
DataManager::DataManager() {
    // Inicializaciones necesarias.
}

DataManager::~DataManager() {
    // Liberar recursos si es necesario.
}

// Función para obtener la extensión en minúsculas de un nombre de archivo
static std::string getExtension(const std::string& filename) {
    std::size_t pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    std::string ext = filename.substr(pos);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

// Función para cargar proteínas desde archivos PDB en un directorio
bool DataManager::loadProteins(const std::string& path, std::vector<Molecule>& proteins) {
    DIR *dir = opendir(path.c_str());
    if (dir == nullptr) {
        std::cerr << "Error abriendo el directorio: " << path << std::endl;
        return false;
    }
    struct dirent *entry;
    int count = 0;
    while ((entry = readdir(dir)) != nullptr) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
            continue;

        std::string filepath = path + "/" + entry->d_name;
        struct stat s;
        if (stat(filepath.c_str(), &s) == 0 && S_ISREG(s.st_mode)) {
            std::string ext = getExtension(entry->d_name);
            if (ext == ".pdb") {
                Molecule mol = parsePDB(filepath);
                if (mol.getAtoms().empty()) {
                    std::cerr << "No se parsearon átomos en " << filepath << std::endl;
                } else {
                    proteins.push_back(mol);
                    count++;
                }
            }
        }
    }
    closedir(dir);
    return (count > 0);
}

// Función para cargar ligandos desde archivos SDF o PDB en un directorio
bool DataManager::loadLigands(const std::string& path, std::vector<Molecule>& ligands) {
    DIR *dir = opendir(path.c_str());
    if (dir == nullptr) {
        std::cerr << "Error abriendo el directorio: " << path << std::endl;
        return false;
    }
    struct dirent *entry;
    int count = 0;
    while ((entry = readdir(dir)) != nullptr) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
            continue;

        std::string filepath = path + "/" + entry->d_name;
        struct stat s;
        if (stat(filepath.c_str(), &s) == 0 && S_ISREG(s.st_mode)) {
            std::string ext = getExtension(entry->d_name);
            if (ext == ".pdb") {
                Molecule mol = parsePDB(filepath);
                if (mol.getAtoms().empty())
                    std::cerr << "No se parsearon átomos en " << filepath << std::endl;
                else {
                    ligands.push_back(mol);
                    count++;
                }
            } else if (ext == ".sdf") {
                Molecule mol = parseSDF(filepath);
                if (mol.getAtoms().empty())
                    std::cerr << "No se parsearon átomos en " << filepath << std::endl;
                else {
                    ligands.push_back(mol);
                    count++;
                }
            }
        }
    }
    closedir(dir);
    return (count > 0);
}