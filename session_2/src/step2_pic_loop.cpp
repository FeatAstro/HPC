#include "step2_pic_loop.hpp"

#include "boris_pusher.hpp"
#include "interpolation.hpp"

void advance_particles(std::vector<Particle>& particles, const Grid2D& g, double dt) {
    for (Particle& p : particles) {
        p.x += p.v.x * dt / 2.0;
        p.y += p.v.y * dt / 2.0;

        const Vec3 E = gather(g.E, g, p.x, p.y);
        const Vec3 B = gather(g.B, g, p.x, p.y);

        boris_push(p, E, B, dt);

        p.x += p.v.x * dt / 2.0;
        p.y += p.v.y * dt / 2.0;
    }
}
