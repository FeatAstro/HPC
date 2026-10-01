#include "population.hpp"

#include <cmath>
#include <stdexcept>

#include "boris.hpp"
#include "grid.hpp"

// A tiny negative x can round to exactly `length` after the subtraction, hence the second line.
static double wrap_periodic(double x, double length) {
    x -= length * std::floor(x / length);
    if (x >= length) x -= length;
    return x;
}

void Population::check_sizes() const {
    const std::size_t n = size();
    if (vx.size() != n || vy.size() != n || vz.size() != n || w.size() != n) {
        throw std::invalid_argument("Population '" + name + "': x, vx, vy, vz, w must have the same size");
    }
}

void Population::check_deposited_on(const Grid& grid) const {
    const std::size_t nx = grid.nx;
    if (density.size() != nx || flux.x.size() != nx || kinetic_energy_density.size() != nx) {
        throw std::invalid_argument("Population '" + name + "' has not been deposited on this grid");
    }
}

static void check_fields_match_grid(const Fields& fields, const Grid& grid) {
    const std::size_t nx = grid.nx;
    for (const VectorField* field : {&fields.E, &fields.B}) {
        if (field->x.size() != nx || field->y.size() != nx || field->z.size() != nx) {
            throw std::invalid_argument("Fields must have one value per grid point");
        }
    }
}

void Population::push(const Grid& grid, const Fields& fields, double dt) {
    check_sizes();
    check_fields_match_grid(fields, grid);
    const double length = grid.length();
    for (std::size_t i = 0; i < size(); ++i) {
        const double x_half = wrap_periodic(x[i] + 0.5 * dt * vx[i], length);

        const ParticleShape shape = particle_shape(grid, x_half);
        const Vec3 E = gather(fields.E, shape);
        const Vec3 B = gather(fields.B, shape);
        const Vec3 v = boris_velocity_update({vx[i], vy[i], vz[i]}, E, B, q, m, dt);
        vx[i] = v.x;
        vy[i] = v.y;
        vz[i] = v.z;

        x[i] = wrap_periodic(x_half + 0.5 * dt * v.x, length);
    }
}

void Population::deposit(const Grid& grid) {
    check_sizes();
    density.assign(grid.nx, 0.0);
    flux = VectorField(grid);
    kinetic_energy_density.assign(grid.nx, 0.0);
    for (std::size_t i = 0; i < size(); ++i) {
        const ShapeWeights shape = shape_weights(grid, x[i]);
        const Vec3 v = {vx[i], vy[i], vz[i]};
        deposit_onto(density, shape, w[i]);
        deposit_onto(flux.x, shape, w[i] * v.x);
        deposit_onto(flux.y, shape, w[i] * v.y);
        deposit_onto(flux.z, shape, w[i] * v.z);
        deposit_onto(kinetic_energy_density, shape, 0.5 * m * w[i] * dot(v, v));
    }
}

ScalarField Population::temperature() const {
    ScalarField T(density.size(), 0.0);
    for (std::size_t i = 0; i < density.size(); ++i) {
        if (density[i] > 0.0) {
            const Vec3 u = {flux.x[i] / density[i], flux.y[i] / density[i], flux.z[i] / density[i]};
            T[i] = (2.0 / 3.0) * (kinetic_energy_density[i] / density[i] - 0.5 * m * dot(u, u));
        }
    }
    return T;
}

ScalarField total_density(const std::vector<const Population*>& populations, const Grid& grid) {
    ScalarField density(grid.nx, 0.0);
    for (const Population* population : populations) {
        population->check_deposited_on(grid);
        for (int i = 0; i < grid.nx; ++i) {
            density[i] += population->density[i];
        }
    }
    return density;
}

VectorField bulk_velocity(const std::vector<const Population*>& populations, const Grid& grid) {
    const ScalarField density = total_density(populations, grid);
    VectorField velocity(grid);
    for (const Population* population : populations) {
        for (int i = 0; i < grid.nx; ++i) {
            velocity.x[i] += population->flux.x[i];
            velocity.y[i] += population->flux.y[i];
            velocity.z[i] += population->flux.z[i];
        }
    }
    for (int i = 0; i < grid.nx; ++i) {
        if (density[i] > 0.0) {
            velocity.x[i] /= density[i];
            velocity.y[i] /= density[i];
            velocity.z[i] /= density[i];
        }
    }
    return velocity;
}

static void record_totals(PopulationHistory& history, const Population& population, double t) {
    double total_weight = 0.0, kinetic_energy = 0.0;
    Vec3 weighted_velocity_sum;
    for (std::size_t i = 0; i < population.size(); ++i) {
        const double w = population.w[i];
        const Vec3 v = {population.vx[i], population.vy[i], population.vz[i]};
        total_weight += w;
        kinetic_energy += 0.5 * population.m * w * dot(v, v);
        weighted_velocity_sum = weighted_velocity_sum + w * v;
    }
    history.t.push_back(t);
    history.kinetic_energy.push_back(kinetic_energy);
    history.mean_vx.push_back(weighted_velocity_sum.x / total_weight);
    history.mean_vy.push_back(weighted_velocity_sum.y / total_weight);
    history.mean_vz.push_back(weighted_velocity_sum.z / total_weight);
}

PopulationHistory run_population(Population& population, const Grid& grid, const Fields& fields,
                                 double dt, int n_steps) {
    population.check_sizes();
    PopulationHistory history;
    record_totals(history, population, 0.0);
    for (int step = 1; step <= n_steps; ++step) {
        population.push(grid, fields, dt);
        record_totals(history, population, step * dt);
    }
    return history;
}
