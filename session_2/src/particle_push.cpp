#include "particle_push.hpp"

#include "boris_pusher.hpp"
#include "interpolation.hpp"

void drift_particles(std::vector<Particle>& particles, const Grid2D& g, double duration) {
    for (Particle& p : particles) {
        p.x = g.wrap_x(p.x + p.v.x * duration);
        p.y = g.wrap_y(p.y + p.v.y * duration);
    }
}

void accelerate_particles(std::vector<Particle>& particles, const Grid2D& g, double dt) {
    for (Particle& p : particles) {
        const Vec3 E = gather(g.electric_field, yee::electric_field, g, p.x, p.y);
        const Vec3 B = gather(g.magnetic_field, yee::magnetic_field, g, p.x, p.y);
        boris_push(p, E, B, dt);
    }
}
