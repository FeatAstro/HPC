// Yee layout tests: every component of E, B and j is interpolated from the position
// where it lives, so a linear function of position is gathered exactly.
#include <iostream>
#include <random>

#include "check.hpp"
#include "grid.hpp"
#include "interpolation.hpp"

namespace {

Vec3 linear_field(double x, double y) {
    return Vec3{1.0 + 2.0 * x - 0.5 * y, -3.0 + 0.7 * x + 1.5 * y, 0.25 - 1.1 * x - 0.3 * y};
}

void check_layout_gathers_linear_field(const Grid2D& grid, const VectorField& field, const Staggering& layout,
                                       const char* what) {
    std::mt19937 random_generator(5);
    // stay away from the periodic seam, where a linear function is not periodic
    std::uniform_real_distribution<double> x_position(1.0 * grid.dx, (grid.nx - 2) * grid.dx);
    std::uniform_real_distribution<double> y_position(1.0 * grid.dy, (grid.ny - 2) * grid.dy);

    for (int trial = 0; trial < 200; ++trial) {
        const double x = x_position(random_generator);
        const double y = y_position(random_generator);
        const Vec3 gathered = gather(field, layout, grid, x, y);
        const Vec3 expected = linear_field(x, y);
        check_close(gathered.x, expected.x, 1e-10, what);
        check_close(gathered.y, expected.y, 1e-10, what);
        check_close(gathered.z, expected.z, 1e-10, what);
    }
}

void test_staggered_gather_is_exact_for_linear_fields() {
    Grid2D grid(12, 9, 0.5, 0.8);

    grid.fill_field(grid.electric_field, yee::electric_field, linear_field);
    grid.fill_field(grid.magnetic_field, yee::magnetic_field, linear_field);
    grid.fill_field(grid.current_density, yee::current_density, linear_field);

    check_layout_gathers_linear_field(grid, grid.electric_field, yee::electric_field, "E gather exact");
    check_layout_gathers_linear_field(grid, grid.magnetic_field, yee::magnetic_field, "B gather exact");
    check_layout_gathers_linear_field(grid, grid.current_density, yee::current_density, "j gather exact");
    std::cout << "test_staggered_gather_is_exact_for_linear_fields: passed\n";
}

void test_wrong_layout_is_detected() {
    Grid2D grid(12, 9, 0.5, 0.8);
    grid.fill_field(grid.electric_field, yee::electric_field, linear_field);

    // Reading E as if it were node-centred moves Ex by half a cell (0.25 in x), i.e. 0.5 off for a slope of 2.
    const double x = 3.0;
    const double y = 3.0;
    const Vec3 wrong = gather(grid.electric_field, yee::nodes, grid, x, y);
    check_true(std::fabs(wrong.x - linear_field(x, y).x) > 0.1, "ignoring the staggering must give a visibly wrong Ex");
    std::cout << "test_wrong_layout_is_detected: passed\n";
}

}

int main() {
    test_staggered_gather_is_exact_for_linear_fields();
    test_wrong_layout_is_detected();
    std::cout << "test_yee_grid: all checks passed\n";
    return 0;
}
