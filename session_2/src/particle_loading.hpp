#pragma once

#include <vector>

#include "grid.hpp"
#include "particle.hpp"
#include "plasma_parameters.hpp"

struct IonLoading {
    int particles_per_cell_x = 2;
    int particles_per_cell_y = 2;
    double thermal_speed = 0.0; // standard deviation of each velocity component; 0 = cold
    unsigned random_seed = 1;
};

// Fills the domain with ions on a regular lattice (noise-free density) at
// plasma.reference_density, each with statistical weight density / particles_per_cell.
std::vector<Particle> load_uniform_ions(const Grid2D& g, const PlasmaParameters& plasma, const IonLoading& loading);
