#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <iostream>
#include <string>

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

#endif // UTILS_H
