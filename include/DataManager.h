#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <string>
#include <vector>
#include "Molecule.h"

class DataManager {
public:
    DataManager();
    ~DataManager();

    // Carga proteínas desde un directorio o fichero
    bool loadProteins(const std::string& path, std::vector<Molecule>& proteins);

    // Carga ligandos desde un directorio o fichero
    bool loadLigands(const std::string& path, std::vector<Molecule>& ligands);
};

#endif // DATAMANAGER_H