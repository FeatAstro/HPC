#pragma once

#include <vector>

#include "vec3.hpp"

// 1D3V: position x only, velocity in 3D because v x B needs all three components.
struct Particle {
    double x = 0.0;
    Vec3 v;
    double q = 1.0;
    double m = 1.0;
};

// Velocity part of the Boris step (slide 7); E and B are taken at the mid-step position.
inline Vec3 boris_velocity_update(const Vec3& v, const Vec3& E, const Vec3& B, double q, double m, double dt) {
    const double half_q_dt_over_m = 0.5 * q * dt / m;

    const Vec3 v_minus = v + half_q_dt_over_m * E;

    // rotation around B, |v_plus| = |v_minus| exactly
    const Vec3 t = half_q_dt_over_m * B;
    const Vec3 s = (2.0 / (1.0 + dot(t, t))) * t;
    const Vec3 v_prime = v_minus + cross(v_minus, t);
    const Vec3 v_plus = v_minus + cross(v_prime, s);

    return v_plus + half_q_dt_over_m * E;
}

// Full Boris step in uniform fields (slide 7): half drift, velocity update, half drift.
void boris_push(Particle& particle, const Vec3& E, const Vec3& B, double dt);

// Index 0 is the initial state.
struct Trajectory {
    std::vector<double> t, x, vx, vy, vz;
};

Trajectory run_particle(Particle& particle, const Vec3& E, const Vec3& B, double dt, int n_steps);
