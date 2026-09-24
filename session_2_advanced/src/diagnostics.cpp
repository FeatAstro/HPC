#include "diagnostics.hpp"

#include <algorithm>
#include <cmath>

#include "finite_differences.hpp"

double ion_kinetic_energy(const Grid2D& g, const std::vector<Particle>& ions) {
    double energy = 0.0;
    for (const Particle& ion : ions) {
        energy += 0.5 * ion.m * dot(ion.v, ion.v) * ion.w * g.dx * g.dy;
    }
    return energy;
}

double magnetic_energy(const Grid2D& g, const PlasmaParameters& plasma) {
    double energy = 0.0;
    for (const Vec3& B : g.magnetic_field) {
        energy += dot(B, B) / (2.0 * plasma.vacuum_permeability) * g.dx * g.dy;
    }
    return energy;
}

double max_divergence_of_magnetic_field(const Grid2D& g) {
    const auto Bx = component_of(g.magnetic_field, axis::x);
    const auto By = component_of(g.magnetic_field, axis::y);

    double largest = 0.0;
    g.for_each_node([&](int i, int j, int) {
        const double divergence = forward_difference_x(g, Bx, i, j) + forward_difference_y(g, By, i, j);
        largest = std::max(largest, std::fabs(divergence));
    });
    return largest;
}
