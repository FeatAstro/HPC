#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>

inline void check_close(double actual, double expected, double tolerance, const char* what) {
    if (std::fabs(actual - expected) > tolerance) {
        std::cerr << "FAIL: " << what << " expected " << expected << " got " << actual << " (tolerance "
                  << tolerance << ")\n";
        std::exit(1);
    }
}

inline void check_true(bool condition, const char* what) {
    if (!condition) {
        std::cerr << "FAIL: " << what << "\n";
        std::exit(1);
    }
}
