#include "step1_single_particle.hpp"

#include "boris_pusher.hpp"

void advance_particle(Particle& p, const FieldFunc& E_of, const FieldFunc& B_of, double dt) {
    p.x += p.v.x * dt / 2.0;
    p.y += p.v.y * dt / 2.0;

    const Vec3 E = E_of(p.x, p.y);
    const Vec3 B = B_of(p.x, p.y);

    boris_push(p, E, B, dt);

    p.x += p.v.x * dt / 2.0;
    p.y += p.v.y * dt / 2.0;
}
