#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "grid.hpp"
#include "vec3.hpp"

// One ion species and its macro-particles, stored as one array per quantity
// (structure of arrays): particle i is (x[i], vx[i], vy[i], vz[i], w[i]).
struct Population {
    std::string name;
    double q = 1.0;
    double m = 1.0;

    std::vector<double> x;
    std::vector<double> vx, vy, vz;
    std::vector<double> w;  // statistical weight: the density this macro-particle deposits

    // Moments on the grid nodes, filled by deposit() (slide 4: n_i = sum_p S(x_p - x_i) w_p).
    ScalarField density;
    VectorField flux;                    // sum_p S w_p v_p
    ScalarField kinetic_energy_density;  // sum_p S w_p m |v_p|^2 / 2

    Population(const std::string& name, double q, double m) : name(name), q(q), m(m) {}

    std::size_t size() const { return x.size(); }

    // Boris step (slide 7) with E and B gathered at the mid-step position.
    // Also wraps x periodically into [0, grid.length()).
    void push(const Grid& grid, const Fields& fields, double dt);

    void deposit(const Grid& grid);

    // From the deposited moments: T = (2/3) (kinetic energy density / n - m |u|^2 / 2), u = flux / n.
    // Zero on a node without particles.
    ScalarField temperature() const;

    void check_sizes() const;
    void check_deposited_on(const Grid& grid) const;
};

// Sums over all populations (slide 5).
ScalarField total_density(const std::vector<const Population*>& populations, const Grid& grid);
VectorField bulk_velocity(const std::vector<const Population*>& populations, const Grid& grid);

// Totals weighted by w; index 0 is the initial state.
struct PopulationHistory {
    std::vector<double> t;
    std::vector<double> kinetic_energy;
    std::vector<double> mean_vx, mean_vy, mean_vz;
};

PopulationHistory run_population(Population& population, const Grid& grid, const Fields& fields,
                                 double dt, int n_steps);
