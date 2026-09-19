#pragma once

#include <cstddef>
#include <vector>

#include "vec3.hpp"

using ScalarField = std::vector<double>;
using VectorField = std::vector<Vec3>;

inline VectorField scaled(const VectorField& field, double factor) {
    VectorField result(field.size());
    for (std::size_t k = 0; k < field.size(); ++k) {
        result[k] = factor * field[k];
    }
    return result;
}

// base + factor * increment
inline VectorField added(const VectorField& base, double factor, const VectorField& increment) {
    VectorField result(base.size());
    for (std::size_t k = 0; k < base.size(); ++k) {
        result[k] = base[k] + factor * increment[k];
    }
    return result;
}

inline VectorField averaged(const VectorField& first, const VectorField& second) {
    return added(scaled(first, 0.5), 0.5, second);
}
