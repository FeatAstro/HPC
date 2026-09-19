// Step 1 test: single particle, uniform B, no E -> pure gyration.
// Checks: |v| conservation (exact, Boris rotation is orthogonal), constant
// gyroradius about the analytic guiding center, and return-to-start after
// one full period.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>

#include "particle.hpp"
#include "step1_single_particle.hpp"
#include "vec3.hpp"

namespace {
void check_close(double a, double b, double tol, const char* what) {
    if (std::fabs(a - b) > tol) {
        std::cerr << "FAIL: " << what << " expected " << b << " got " << a << " (tol " << tol << ")\n";
        std::exit(1);
    }
}
} // namespace

int main() {
    const double B0 = 1.0;
    const auto E_of = [](double, double) { return Vec3{0.0, 0.0, 0.0}; };
    const auto B_of = [B0](double, double) { return Vec3{0.0, 0.0, B0}; };

    Particle p;
    p.q = 1.0;
    p.m = 1.0;
    p.x = 0.0;
    p.y = 0.0;
    p.v = Vec3{1.0, 0.0, 0.0};

    const double x0 = p.x;
    const double y0 = p.y;
    const double v0 = norm(p.v);
    const double omega_c = p.q * B0 / p.m;
    const double r_L = v0 / omega_c;
    const double yc = y0 - r_L; // analytic guiding-center y (see derivation in plan)

    const double speed0 = norm(p.v);

    const int steps_per_period = 2000;
    const double period = 2.0 * std::numbers::pi / omega_c;
    const double dt = period / steps_per_period;

    double max_radius_error = 0.0;
    for (int step = 0; step < steps_per_period; ++step) {
        advance_particle(p, E_of, B_of, dt);

        check_close(norm(p.v), speed0, 1e-9, "|v| conservation");

        const double r = std::sqrt((p.x - x0) * (p.x - x0) + (p.y - yc) * (p.y - yc));
        max_radius_error = std::max(max_radius_error, std::fabs(r - r_L));
    }

    check_close(max_radius_error, 0.0, 1e-3, "gyroradius about guiding center");
    check_close(p.x, x0, 1e-3, "return-to-start x after one period");
    check_close(p.y, y0, 1e-3, "return-to-start y after one period");

    std::cout << "test_boris_pusher: all checks passed\n";
    return 0;
}
