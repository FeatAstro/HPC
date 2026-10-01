#pragma once

#include <functional>

#include "grid.hpp"
#include "population.hpp"
#include "vec3.hpp"

// A profile is any function of x: a C++ function, a lambda, or a Python function.
using Profile = std::function<double(double)>;
using VectorProfile = std::function<Vec3(double)>;

// Replaces the particles of the population by a Maxwellian plasma with the given profiles.
// In every cell: particles_per_cell particles at uniform random positions, all with the density,
// bulk velocity and temperature of the cell centre (as in PHARE), so w = n / particles_per_cell
// and each velocity component is u + sqrt(T/m) * (standard normal number).
// The same seed gives the same particles.
void load_maxwellian(Population& population, const Grid& grid, int particles_per_cell,
                     const Profile& density, const VectorProfile& bulk_velocity, const Profile& temperature,
                     unsigned seed);
