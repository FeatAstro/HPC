#include "boris.hpp"

void boris_push(Particle& particle, const Vec3& E, const Vec3& B, double dt) {
    particle.x += 0.5 * dt * particle.v.x;
    particle.v = boris_velocity_update(particle.v, E, B, particle.q, particle.m, dt);
    particle.x += 0.5 * dt * particle.v.x;
}

static void record_state(Trajectory& trajectory, const Particle& particle, double t) {
    trajectory.t.push_back(t);
    trajectory.x.push_back(particle.x);
    trajectory.vx.push_back(particle.v.x);
    trajectory.vy.push_back(particle.v.y);
    trajectory.vz.push_back(particle.v.z);
}

Trajectory run_particle(Particle& particle, const Vec3& E, const Vec3& B, double dt, int n_steps) {
    Trajectory trajectory;
    record_state(trajectory, particle, 0.0);
    for (int step = 1; step <= n_steps; ++step) {
        boris_push(particle, E, B, dt);
        record_state(trajectory, particle, step * dt);
    }
    return trajectory;
}
