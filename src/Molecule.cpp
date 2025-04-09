#include "Molecule.h"

Molecule::Molecule() {
    // Inicialización del objeto molecule.
}

Molecule::~Molecule() {
    // Destructor.
}

void Molecule::addAtom(const Atom& atom) {
    atoms.push_back(atom);
}

const std::vector<Atom>& Molecule::getAtoms() const {
    return atoms;
}