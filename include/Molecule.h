#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <vector>

// Estructura para representar un átomo (muy simplificada)
struct Atom {
    float x, y, z;
    std::string element;
};

class Molecule {
public:
    Molecule();
    ~Molecule();

    void addAtom(const Atom& atom);
    const std::vector<Atom>& getAtoms() const;

private:
    std::vector<Atom> atoms;
};

#endif // MOLECULE_H