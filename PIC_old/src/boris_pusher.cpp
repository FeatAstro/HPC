#include "boris_pusher.hpp"

void boris_push(Particle& p, const Vec3& E, const Vec3& B, double dt) {
    const double qmdt2 = 0.5 * p.q * dt / p.m;

    const Vec3 v_minus = p.v + qmdt2 * E;

    const Vec3 t = qmdt2 * B;
    const Vec3 v_prime = v_minus + cross(v_minus, t);
    const Vec3 s = (2.0 / (1.0 + dot(t, t))) * t;
    const Vec3 v_plus = v_minus + cross(v_prime, s);

    p.v = v_plus + qmdt2 * E;
}
