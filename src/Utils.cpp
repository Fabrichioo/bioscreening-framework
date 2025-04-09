#include "Utils.h"

Timer::Timer() {}

void Timer::start() {
    m_start = std::chrono::high_resolution_clock::now();
}

void Timer::stop() {
    m_end = std::chrono::high_resolution_clock::now();
}

double Timer::elapsedMilliseconds() const {
    return std::chrono::duration<double, std::milli>(m_end - m_start).count();
}

void logMessage(const std::string& msg) {
    std::cout << msg << std::endl;
}