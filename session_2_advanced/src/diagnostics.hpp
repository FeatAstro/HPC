#pragma once

#include <vector>

#include "grid.hpp"
#include "particle.hpp"
#include "plasma_parameters.hpp"

// A macro-particle stands for w_p * dx * dy real ions (w_p is a density contribution).
double ion_kinetic_energy(const Grid2D& g, const std::vector<Particle>& ions);

double magnetic_energy(const Grid2D& g, const PlasmaParameters& plasma);

// div B at the cell centres; stays at round-off level with the Yee scheme.
double max_divergence_of_magnetic_field(const Grid2D& g);
