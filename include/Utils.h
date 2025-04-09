#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include "Docking.h"  // Para que se conozca la definición de DockingResult

class Timer {
public:
    Timer();
    void start();
    void stop();
    double elapsedMilliseconds() const;
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
};

void logMessage(const std::string& msg);

// Función para analizar y mostrar los resultados del docking:
// Recibe el vector de scores, el número de proteínas y el número de ligandos.
void analyzeDockingResults(const std::vector<float>& scores, int numProteins, int numLigands);

#endif // UTILS_H