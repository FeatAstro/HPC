#pragma once

#include <chrono>

// Measures wall-clock time. steady_clock is a clock that only goes forward, made for durations.
class Stopwatch {
public:
    // Seconds since the stopwatch was created or since the previous lap, then starts a new lap.
    double lap_seconds() {
        const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = now - lap_start;
        lap_start = now;
        return elapsed.count();
    }

private:
    std::chrono::steady_clock::time_point lap_start = std::chrono::steady_clock::now();
};
