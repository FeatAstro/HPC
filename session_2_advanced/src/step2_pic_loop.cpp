#include "step2_pic_loop.hpp"

#include "particle_push.hpp"

void advance_particles(std::vector<Particle>& particles, const Grid2D& g, double dt) {
    drift_particles(particles, g, dt / 2.0);
    accelerate_particles(particles, g, dt);
    drift_particles(particles, g, dt / 2.0);
}
