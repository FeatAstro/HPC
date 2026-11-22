#pragma once

#include <functional>
#include <vector>

#include "grid.hpp"

// A population: one species of ions (here protons) as macro-particles, one array per quantity:
// particle i is (x[i], vx[i], vy[i], vz[i], w[i]).
// Position in 1D, velocity in 3D because v x B needs all three components.
struct Population {
    double q = 1.0;
    double m = 1.0;

    std::vector<double> x;
    std::vector<double> vx, vy, vz;
    std::vector<double> w;  // statistical weight: the density this macro-particle deposits
};

// The ion fluid seen from the grid, on the nodes (slide 4: n_i = sum_p S(x_p - x_i) w_p).
struct Moments {
    ScalarField density;
    VectorField bulk_velocity;  // u = sum_p S w_p v_p / n
    ScalarField temperature;    // T = (m / 3) (sum_p S w_p |v_p|^2 / n - |u|^2): the spread of v around u
};

// A profile is any function of x given by the user: a C++ function, a lambda, or a Python function.
using Profile = std::function<double(double)>;

// Protons (q = m = 1) loaded from profiles of density n(x), bulk velocity (ux, uy, uz)(x) and
// temperature T(x). In every cell, particles_per_cell protons at random positions; each one takes
// the profiles at its own position: w = n(x) / particles_per_cell, and each velocity component is
// u(x) + sqrt(T(x) / m) times a standard normal number (a Maxwellian). The same seed gives the
// same protons.
Population maxwellian_protons(const Grid& grid, int particles_per_cell, const Profile& density, const Profile& ux,
                        const Profile& uy, const Profile& uz, const Profile& temperature, unsigned seed);

// Boris step (slide 7), with E and B taken at the mid-step position; x stays in [0, L).
void push(Population& ions, const Fields& fields, const Grid& grid, double dt);

Moments deposit(const Population& ions, const Grid& grid);
