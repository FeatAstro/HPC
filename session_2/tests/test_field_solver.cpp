// Field solver tests: curl operators (second-order accuracy on the Yee layout),
// Ohm's law terms one at a time, and div B preservation by the Faraday update.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>

#include "check.hpp"
#include "diagnostics.hpp"
#include "field_solver.hpp"
#include "grid.hpp"

namespace {

constexpr double two_pi = 2.0 * std::numbers::pi;

Grid2D periodic_square_grid(int cells) {
    return Grid2D(cells, cells, two_pi / cells, two_pi / cells);
}

Vec3 smooth_vector_field(double x, double y) {
    return Vec3{std::sin(y), std::sin(x) * std::cos(y), std::cos(x + y)};
}

Vec3 curl_of_smooth_vector_field(double x, double y) {
    return Vec3{-std::sin(x + y), std::sin(x + y), std::cos(x) * std::cos(y) - std::cos(y)};
}

using CurlOperator = VectorField (*)(const Grid2D&, const VectorField&);

double max_curl_error(int cells, const Staggering& input_layout, const Staggering& output_layout,
                      CurlOperator curl_operator) {
    const Grid2D grid = periodic_square_grid(cells);
    VectorField field(grid.ion_density.size());
    grid.fill_field(field, input_layout, smooth_vector_field);

    const VectorField numerical_curl = curl_operator(grid, field);

    double largest_error = 0.0;
    grid.for_each_node([&](int i, int j, int node) {
        for (int component = 0; component < axis::count; ++component) {
            const double x = grid.node_x(i, output_layout[component].x);
            const double y = grid.node_y(j, output_layout[component].y);
            const double exact = curl_of_smooth_vector_field(x, y)[component];
            largest_error = std::max(largest_error, std::fabs(numerical_curl[node][component] - exact));
        }
    });
    return largest_error;
}

void check_second_order(const Staggering& input_layout, const Staggering& output_layout, CurlOperator curl_operator,
                        const char* what) {
    const double coarse_error = max_curl_error(16, input_layout, output_layout, curl_operator);
    const double fine_error = max_curl_error(32, input_layout, output_layout, curl_operator);
    const double convergence_ratio = coarse_error / fine_error;

    check_true(fine_error < 5e-3, what);
    check_true(convergence_ratio > 3.5 && convergence_ratio < 4.5, what);
}

void test_curls_are_second_order() {
    check_second_order(yee::electric_field, yee::magnetic_field, curl_of_electric_field, "curl E convergence");
    check_second_order(yee::magnetic_field, yee::current_density, curl_of_magnetic_field, "curl B convergence");
    std::cout << "test_curls_are_second_order: passed\n";
}

void test_ohm_convective_term() {
    Grid2D grid(16, 8, 0.5, 0.5);
    PlasmaParameters plasma;
    plasma.electron_temperature = 0.7; // must not matter for a uniform density

    const Vec3 ion_velocity{0.2, -0.1, 0.3};
    const Vec3 magnetic_field{0.3, 0.4, 0.5};
    std::fill(grid.ion_density.begin(), grid.ion_density.end(), 1.0);
    std::fill(grid.ion_velocity.begin(), grid.ion_velocity.end(), ion_velocity);
    std::fill(grid.magnetic_field.begin(), grid.magnetic_field.end(), magnetic_field);

    update_electric_field(grid, plasma);

    const Vec3 expected = -1.0 * cross(ion_velocity, magnetic_field);
    for (const Vec3& E : grid.electric_field) {
        check_close(E.x, expected.x, 1e-12, "E = -v_i x B, x");
        check_close(E.y, expected.y, 1e-12, "E = -v_i x B, y");
        check_close(E.z, expected.z, 1e-12, "E = -v_i x B, z");
    }
    std::cout << "test_ohm_convective_term: passed\n";
}

void test_ohm_electron_pressure_term() {
    Grid2D grid(64, 8, two_pi / 64, two_pi / 64);
    PlasmaParameters plasma;
    plasma.electron_temperature = 0.5;

    const double amplitude = 0.2;
    grid.for_each_node([&](int i, int, int node) { grid.ion_density[node] = 1.0 + amplitude * std::sin(grid.node_x(i)); });

    update_electric_field(grid, plasma);

    // E_x = -T_e (dn/dx) / (n e), evaluated where E_x lives (half a cell to the right of node i)
    grid.for_each_node([&](int i, int, int node) {
        const double x = grid.node_x(i, yee::electric_field[axis::x].x);
        const double density = 1.0 + amplitude * std::sin(x);
        const double expected = -plasma.electron_temperature * amplitude * std::cos(x) / density;
        check_close(grid.electric_field[node].x, expected, 2e-4, "E_x from the electron pressure gradient");
        check_close(grid.electric_field[node].y, 0.0, 1e-12, "no pressure gradient along y");
        check_close(grid.electric_field[node].z, 0.0, 1e-12, "no pressure gradient along z");
    });
    std::cout << "test_ohm_electron_pressure_term: passed\n";
}

void test_ohm_hall_term() {
    Grid2D grid(64, 8, two_pi / 64, two_pi / 64);
    const PlasmaParameters plasma;

    const double background = 1.0;
    const double perturbation = 0.3;
    std::fill(grid.ion_density.begin(), grid.ion_density.end(), 1.0);
    grid.fill_field(grid.magnetic_field, yee::magnetic_field, [&](double x, double) {
        return Vec3{background, 0.0, perturbation * std::sin(x)};
    });

    update_electric_field(grid, plasma);

    // Ions at rest, so E = (j x B) / (n e) with j = (0, -perturbation cos x, 0).
    grid.for_each_node([&](int i, int, int node) {
        const double x_of_ex = grid.node_x(i, yee::electric_field[axis::x].x);
        const double x_of_ez = grid.node_x(i, yee::electric_field[axis::z].x);
        check_close(grid.electric_field[node].x, -perturbation * perturbation * std::cos(x_of_ex) * std::sin(x_of_ex),
                    5e-4, "Hall term, E_x");
        check_close(grid.electric_field[node].y, 0.0, 1e-12, "Hall term, E_y");
        check_close(grid.electric_field[node].z, perturbation * std::cos(x_of_ez) * background, 5e-4, "Hall term, E_z");
    });
    std::cout << "test_ohm_hall_term: passed\n";
}

void test_faraday_update_preserves_div_b() {
    Grid2D grid = periodic_square_grid(32);
    PlasmaParameters plasma;
    plasma.electron_temperature = 0.3;

    // B = curl(A_z z_hat) with A_z on the nodes is divergence-free by construction on the Yee grid.
    VectorField vector_potential(grid.ion_density.size());
    grid.for_each_node([&](int i, int j, int node) {
        vector_potential[node].z = std::cos(grid.node_x(i)) * std::cos(2.0 * grid.node_y(j));
    });
    grid.magnetic_field = curl_of_electric_field(grid, vector_potential);
    check_true(max_divergence_of_magnetic_field(grid) < 1e-10, "initial B is divergence-free");

    // arbitrary non-uniform ion moments so that E has every kind of term
    grid.for_each_node([&](int i, int j, int node) {
        const double x = grid.node_x(i);
        const double y = grid.node_y(j);
        grid.ion_density[node] = 1.0 + 0.3 * std::sin(x) * std::cos(y);
        grid.ion_velocity[node] = Vec3{std::sin(y), std::cos(x), 0.5};
    });

    const VectorField before = grid.magnetic_field;
    for (int step = 0; step < 20; ++step) {
        advance_magnetic_field(grid, plasma, 0.01);
    }

    double largest_change = 0.0;
    for (std::size_t node = 0; node < before.size(); ++node) {
        largest_change = std::max(largest_change, norm(grid.magnetic_field[node] - before[node]));
    }
    check_true(largest_change > 1e-4, "B must actually evolve");
    check_true(max_divergence_of_magnetic_field(grid) < 1e-10, "div B stays at round-off after Faraday updates");
    std::cout << "test_faraday_update_preserves_div_b: passed\n";
}

}

int main() {
    test_curls_are_second_order();
    test_ohm_convective_term();
    test_ohm_electron_pressure_term();
    test_ohm_hall_term();
    test_faraday_update_preserves_div_b();
    std::cout << "test_field_solver: all checks passed\n";
    return 0;
}
