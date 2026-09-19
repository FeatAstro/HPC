// Step 2 tests: bilinear deposit/gather correctness.
//  1. weight conservation: total deposited weight on the grid == total particle weight
//  2. exactness on linear fields: gather must reproduce a linear f(x,y) exactly
//  3. consistency with Step 1: N=1 particle, uniform prescribed field -> same
//     trajectory as the Step 1 Boris pusher test.
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "grid.hpp"
#include "interpolation.hpp"
#include "particle.hpp"
#include "step1_single_particle.hpp"
#include "step2_pic_loop.hpp"
#include "vec3.hpp"
#include "check.hpp"
#include "particle_push.hpp"

namespace {

void test_weight_conservation() {
    Grid2D grid(11, 11, 1.0, 1.0, 0.0, 0.0);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> pos_dist(0.0, 10.0);

    std::vector<Particle> particles(500);
    double total_weight = 0.0;
    for (Particle& p : particles) {
        p.x = pos_dist(rng);
        p.y = pos_dist(rng);
        p.w = 1.0;
        total_weight += p.w;
    }
    // include an exact-corner and exact-edge particle
    particles.push_back(Particle{0.0, 0.0, Vec3{}, 1.0, 1.0, 1.0});
    particles.push_back(Particle{10.0, 10.0, Vec3{}, 1.0, 1.0, 1.0});
    total_weight += 2.0;

    deposit_moments(grid, particles);

    double deposited_weight = 0.0;
    for (double n_ij : grid.ion_density) {
        deposited_weight += n_ij;
    }

    check_close(deposited_weight, total_weight, 1e-9, "weight conservation");
    std::cout << "test_weight_conservation: passed\n";
}

void test_nonunit_spacing() {
    // dx = 0.5, dy = 4.0: domain x in [0,5], y in [0,20]
    Grid2D grid(11, 6, 0.5, 4.0, 0.0, 0.0);

    // (a) weight conservation with non-unit spacing
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> xdist(0.0, 5.0);
    std::uniform_real_distribution<double> ydist(0.0, 20.0);
    std::vector<Particle> particles(300);
    double total_weight = 0.0;
    for (Particle& p : particles) {
        p.x = xdist(rng);
        p.y = ydist(rng);
        p.w = 2.0;
        total_weight += p.w;
    }
    deposit_moments(grid, particles);
    double deposited_weight = 0.0;
    for (double n_ij : grid.ion_density) {
        deposited_weight += n_ij;
    }
    check_close(deposited_weight, total_weight, 1e-9, "weight conservation, dx!=1");

    // (b) exact node value and bulk velocity at a node
    // two particles sitting exactly on node (2,2) = (1.0, 8.0), weights 1 and 3
    Particle a{1.0, 8.0, Vec3{1.0, 0.0, 0.0}, 1.0, 1.0, 1.0};
    Particle b{1.0, 8.0, Vec3{5.0, 0.0, 0.0}, 1.0, 1.0, 3.0};
    deposit_moments(grid, {a, b});
    const int k = grid.index(2, 2);
    check_close(grid.ion_density[k], 1.0 + 3.0, 1e-12, "n = summed weight, independent of spacing");
    check_close(grid.ion_velocity[k].x, (1.0 * 1.0 + 3.0 * 5.0) / 4.0, 1e-12, "weighted mean velocity");

    std::cout << "test_nonunit_spacing: passed\n";
}

void test_linear_field_exactness() {
    Grid2D grid(6, 6, 1.0, 1.0, 0.0, 0.0);

    const double a = 2.0, b = 0.7, c = -1.3;
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            const double val = a + b * grid.node_x(i) + c * grid.node_y(j);
            grid.electric_field[grid.index(i, j)] = Vec3{val, 0.0, 0.0};
        }
    }

    std::mt19937 rng(11);
    std::uniform_real_distribution<double> pos_dist(0.0, 5.0);
    for (int trial = 0; trial < 200; ++trial) {
        const double x = pos_dist(rng);
        const double y = pos_dist(rng);
        const Vec3 gathered = gather(grid.electric_field, yee::nodes, grid, x, y);
        const double expected = a + b * x + c * y;
        check_close(gathered.x, expected, 1e-9, "linear field exactness");
    }
    std::cout << "test_linear_field_exactness: passed\n";
}

void test_consistency_with_step1() {
    const double B0 = 1.0;
    const auto E_of = [](double, double) { return Vec3{0.0, 0.0, 0.0}; };
    const auto B_of = [B0](double, double) { return Vec3{0.0, 0.0, B0}; };

    Particle p_step1;
    p_step1.q = 1.0;
    p_step1.m = 1.0;
    p_step1.x = 5.0;
    p_step1.y = 5.0;
    p_step1.v = Vec3{1.0, 0.0, 0.0};

    Grid2D grid(11, 11, 1.0, 1.0, 0.0, 0.0);
    for (int k = 0; k < grid.nx * grid.ny; ++k) {
        grid.electric_field[k] = Vec3{0.0, 0.0, 0.0};
        grid.magnetic_field[k] = Vec3{0.0, 0.0, B0};
    }

    std::vector<Particle> particles(1, p_step1);

    const double dt = 0.001;
    const int nsteps = 500;
    for (int step = 0; step < nsteps; ++step) {
        advance_particle(p_step1, E_of, B_of, dt);
        advance_particles(particles, grid, dt);
    }

    check_close(particles[0].x, p_step1.x, 1e-9, "step2 x matches step1");
    check_close(particles[0].y, p_step1.y, 1e-9, "step2 y matches step1");
    check_close(particles[0].v.x, p_step1.v.x, 1e-9, "step2 vx matches step1");
    check_close(particles[0].v.y, p_step1.v.y, 1e-9, "step2 vy matches step1");

    std::cout << "test_consistency_with_step1: passed\n";
}

} // namespace

int main() {
    test_weight_conservation();
    test_nonunit_spacing();
    test_linear_field_exactness();
    test_consistency_with_step1();
    std::cout << "test_deposit_gather: all checks passed\n";
    return 0;
}
