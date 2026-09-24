// Periodic boundary tests: index and position wrapping, particles crossing the seam,
// and deposit/gather continuity across it.
#include <cmath>
#include <iostream>
#include <numbers>
#include <vector>

#include "check.hpp"
#include "grid.hpp"
#include "interpolation.hpp"
#include "particle.hpp"
#include "particle_push.hpp"

namespace {

void test_index_wrapping() {
    Grid2D grid(8, 6, 1.0, 1.0);
    check_true(grid.index(-1, 0) == grid.index(7, 0), "index wraps to the right edge");
    check_true(grid.index(8, 0) == grid.index(0, 0), "index wraps to the left edge");
    check_true(grid.index(0, -1) == grid.index(0, 5), "index wraps to the top edge");
    check_true(grid.index(0, 6) == grid.index(0, 0), "index wraps to the bottom edge");
    check_true(grid.index(-9, 13) == grid.index(7, 1), "index wraps by more than one period");
    std::cout << "test_index_wrapping: passed\n";
}

void test_position_wrapping() {
    Grid2D grid(10, 4, 0.5, 2.0, 1.0, -1.0); // x in [1, 6), y in [-1, 7)
    check_close(grid.length_x(), 5.0, 1e-12, "domain length along x");
    check_close(grid.wrap_x(6.5), 1.5, 1e-12, "position beyond the right edge");
    check_close(grid.wrap_x(0.5), 5.5, 1e-12, "position beyond the left edge");
    check_close(grid.wrap_x(3.0), 3.0, 1e-12, "position inside is unchanged");
    check_close(grid.wrap_y(7.5), -0.5, 1e-12, "position beyond the top edge");
    std::cout << "test_position_wrapping: passed\n";
}

void test_particle_crosses_seam() {
    Grid2D grid(10, 10, 1.0, 1.0);
    std::vector<Particle> particles(1);
    particles[0].x = grid.length_x() - 0.1;
    particles[0].y = 0.1;
    particles[0].v = Vec3{1.0, -1.0, 0.0};

    drift_particles(particles, grid, 0.3);

    check_close(particles[0].x, 0.2, 1e-12, "particle reappears on the left");
    check_close(particles[0].y, grid.length_y() - 0.2, 1e-12, "particle reappears at the top");
    std::cout << "test_particle_crosses_seam: passed\n";
}

void test_deposit_across_seam() {
    Grid2D grid(8, 8, 1.0, 1.0);
    Particle particle;
    particle.x = grid.length_x() - 0.5; // halfway between the last node and node 0
    particle.y = 3.0;                   // exactly on a node row
    particle.w = 2.0;

    deposit_moments(grid, {particle});

    check_close(grid.ion_density[grid.index(7, 3)], 1.0, 1e-12, "half the weight on the last node");
    check_close(grid.ion_density[grid.index(0, 3)], 1.0, 1e-12, "half the weight wraps to node 0");
    std::cout << "test_deposit_across_seam: passed\n";
}

void test_gather_is_periodic() {
    Grid2D grid(16, 8, 0.5, 0.5);
    grid.fill_field(grid.electric_field, [&](double x, double y) {
        return Vec3{std::sin(2.0 * std::numbers::pi * x / grid.length_x()),
                    std::cos(2.0 * std::numbers::pi * y / grid.length_y()), 1.0};
    });

    const double x = grid.length_x() - 0.15;
    const double y = 0.3;
    const Vec3 inside = gather(grid.electric_field, grid, x, y);
    const Vec3 shifted = gather(grid.electric_field, grid, x + grid.length_x(), y - grid.length_y());

    check_close(shifted.x, inside.x, 1e-12, "gather x-component periodic");
    check_close(shifted.y, inside.y, 1e-12, "gather y-component periodic");

    const Vec3 at_seam = gather(grid.electric_field, grid, grid.length_x() - 0.25, 0.0);
    const double mean_of_neighbours = 0.5 * (grid.electric_field[grid.index(15, 0)].x + grid.electric_field[grid.index(0, 0)].x);
    check_close(at_seam.x, mean_of_neighbours, 1e-12, "gather interpolates between last node and node 0");
    std::cout << "test_gather_is_periodic: passed\n";
}

}

int main() {
    test_index_wrapping();
    test_position_wrapping();
    test_particle_crosses_seam();
    test_deposit_across_seam();
    test_gather_is_periodic();
    std::cout << "test_periodic_grid: all checks passed\n";
    return 0;
}
