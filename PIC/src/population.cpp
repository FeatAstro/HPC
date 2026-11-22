#include "population.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

#include <omp.h>

#include "vec3.hpp"

Population maxwellian_protons(const Grid& grid, int particles_per_cell, const Profile& density, const Profile& ux,
                        const Profile& uy, const Profile& uz, const Profile& temperature, unsigned seed) {
    std::mt19937_64 random_numbers(seed);
    std::uniform_real_distribution<double> uniform_in_cell(0.0, 1.0);
    std::normal_distribution<double> standard_normal(0.0, 1.0);

    Population protons;
    for (int cell = 0; cell < grid.nx; ++cell) {
        for (int k = 0; k < particles_per_cell; ++k) {
            const double x = (cell + uniform_in_cell(random_numbers)) * grid.dx;
            const double thermal_speed = std::sqrt(temperature(x) / protons.m);
            protons.x.push_back(x);
            protons.vx.push_back(ux(x) + thermal_speed * standard_normal(random_numbers));
            protons.vy.push_back(uy(x) + thermal_speed * standard_normal(random_numbers));
            protons.vz.push_back(uz(x) + thermal_speed * standard_normal(random_numbers));
            protons.w.push_back(density(x) / particles_per_cell);
        }
    }
    return protons;
}

// Velocity part of the Boris step (slide 7), with the names of the slide.
static Vec3 boris_velocity_update(const Vec3& v, const Vec3& E, const Vec3& B, double q, double m, double dt) {
    const double half_q_dt_over_m = 0.5 * q * dt / m;

    // first half of the push by E
    const Vec3 v_minus = v + half_q_dt_over_m * E;

    // rotation around B by the angle q B dt / m, which keeps the speed exactly
    const Vec3 t = half_q_dt_over_m * B;
    const Vec3 s = (2.0 / (1.0 + dot(t, t))) * t;
    const Vec3 v_prime = v_minus + cross(v_minus, t);
    const Vec3 v_plus = v_minus + cross(v_prime, s);

    // second half of the push by E
    return v_plus + half_q_dt_over_m * E;
}

// Moves x at speed vx during a time t. The box is periodic: leaving it on one side means coming
// back on the other.
static double drift(double x, double vx, double t, const Grid& grid) {
    const double L = grid.length();
    x = std::fmod(x + vx * t, L);  // the remainder of the division by L, in (-L, L)
    if (x < 0.0) x += L;
    if (x >= L) x -= L;  // a tiny negative x can round to exactly L
    return x;
}

// Ex at the cell centres, Ey and Ez on the nodes.
static Vec3 electric_field_at(double x, const VectorField& E, const Grid& grid) {
    return {interpolate(E.x, weights_on_cell_centres(x, grid)),
            interpolate(E.y, weights_on_nodes(x, grid)),
            interpolate(E.z, weights_on_nodes(x, grid))};
}

// Bx on the nodes, By and Bz at the cell centres.
static Vec3 magnetic_field_at(double x, const VectorField& B, const Grid& grid) {
    return {interpolate(B.x, weights_on_nodes(x, grid)),
            interpolate(B.y, weights_on_cell_centres(x, grid)),
            interpolate(B.z, weights_on_cell_centres(x, grid))};
}

// A run that is unstable (time step too large) grows until a velocity overflows to inf or NaN;
// stopping here gives a message instead of a crash when the ions are deposited.
static void check_velocities_are_finite(const Population& ions) {
    for (std::size_t i = 0; i < ions.x.size(); ++i) {
        if (!std::isfinite(ions.vx[i])) {
            throw std::runtime_error("an ion has an infinite velocity: the run has blown up (is dt too large?)");
        }
    }
}

// Each ion drifts half a step, is accelerated by the fields where it is, then drifts the other half.
void push(Population& ions, const Fields& fields, const Grid& grid, double dt) {
    // OpenMP: the particles are shared out between the threads (the cores), each thread pushing its
    // own block of particles. No two threads write to the same particle, so they cannot collide.
    #pragma omp parallel for
    for (std::size_t i = 0; i < ions.x.size(); ++i) {
        const double x_mid_step = drift(ions.x[i], ions.vx[i], 0.5 * dt, grid);

        const Vec3 E = electric_field_at(x_mid_step, fields.E, grid);
        const Vec3 B = magnetic_field_at(x_mid_step, fields.B, grid);
        const Vec3 v_before = {ions.vx[i], ions.vy[i], ions.vz[i]};
        const Vec3 v_after = boris_velocity_update(v_before, E, B, ions.q, ions.m, dt);
        ions.vx[i] = v_after.x;
        ions.vy[i] = v_after.y;
        ions.vz[i] = v_after.z;

        ions.x[i] = drift(x_mid_step, ions.vx[i], 0.5 * dt, grid);
    }
    // After the threads, not in them: an error thrown inside a thread would stop the program at once.
    check_velocities_are_finite(ions);
}

// The sums the moments are made of. Unlike the moments themselves (which are divisions), sums made
// by different threads can simply be added.
struct DepositSums {
    ScalarField weight;         // sum_p S w_p = density
    VectorField flux;           // sum_p S w_p v_p
    ScalarField speed_squared;  // sum_p S w_p |v_p|^2

    DepositSums(const Grid& grid) : weight(grid.nx, 0.0), flux(grid), speed_squared(grid.nx, 0.0) {}

    void add_particle(const ShapeWeights& shape, double w, const Vec3& v) {
        deposit_onto(weight, shape, w);
        deposit_onto(flux.x, shape, w * v.x);
        deposit_onto(flux.y, shape, w * v.y);
        deposit_onto(flux.z, shape, w * v.z);
        deposit_onto(speed_squared, shape, w * dot(v, v));
    }

    void add(const DepositSums& other) {
        for (std::size_t i = 0; i < weight.size(); ++i) {
            weight[i] += other.weight[i];
            flux.x[i] += other.flux.x[i];
            flux.y[i] += other.flux.y[i];
            flux.z[i] += other.flux.z[i];
            speed_squared[i] += other.speed_squared[i];
        }
    }
};

Moments deposit(const Population& ions, const Grid& grid) {
    // OpenMP: two threads depositing particles of the same cell would write to the same node at the
    // same time, and one of the two additions could be lost (a race condition). So each thread
    // deposits onto its own copy of the sums, and the copies are added afterwards (the reduction).
    // schedule(static): each thread always gets the same block of particles, and the copies are added
    // in the order of the threads, so a run gives the same result every time.
    std::vector<DepositSums> sums_of_thread(omp_get_max_threads(), DepositSums(grid));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < ions.x.size(); ++i) {
        const ShapeWeights shape = weights_on_nodes(ions.x[i], grid);
        const Vec3 v = {ions.vx[i], ions.vy[i], ions.vz[i]};
        sums_of_thread[omp_get_thread_num()].add_particle(shape, ions.w[i], v);
    }

    DepositSums sums(grid);
    for (const DepositSums& sums_of_one_thread : sums_of_thread) {
        sums.add(sums_of_one_thread);
    }
    const ScalarField& density = sums.weight;
    const VectorField& flux = sums.flux;
    const ScalarField& speed_squared_sum = sums.speed_squared;

    VectorField bulk_velocity(grid);
    ScalarField temperature(grid.nx);
    for (int i = 0; i < grid.nx; ++i) {
        const Vec3 u = {flux.x[i] / density[i], flux.y[i] / density[i], flux.z[i] / density[i]};
        bulk_velocity.x[i] = u.x;
        bulk_velocity.y[i] = u.y;
        bulk_velocity.z[i] = u.z;
        temperature[i] = (ions.m / 3.0) * (speed_squared_sum[i] / density[i] - dot(u, u));
    }

    return {density, bulk_velocity, temperature};
}
