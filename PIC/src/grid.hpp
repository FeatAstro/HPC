#pragma once

#include <cmath>
#include <vector>

// Periodic 1D grid: nx cells of size dx, so the domain is [0, nx dx).
// Yee layout: a quantity lives either on the nodes, x_i = i dx, or at the cell centres,
// x_i = (i + 1/2) dx. Cell centre i lies between node i and node i + 1.
struct Grid {
    int nx;
    double dx;

    double length() const { return nx * dx; }
};

// One value per node, or one per cell centre.
using ScalarField = std::vector<double>;

struct VectorField {
    ScalarField x, y, z;

    VectorField() = default;
    VectorField(const Grid& grid) : x(grid.nx, 0.0), y(grid.nx, 0.0), z(grid.nx, 0.0) {}
};

// Bx, Ey, Ez live on the nodes (with the ion moments); Ex, By, Bz at the cell centres.
struct Fields {
    VectorField E, B;
};

// The grid repeats itself: point -1 is point nx - 1, and point nx is point 0.
inline int periodic_index(int i, const Grid& grid) {
    if (i < 0) return i + grid.nx;
    if (i >= grid.nx) return i - grid.nx;
    return i;
}

inline int next_point(int i, const Grid& grid) { return periodic_index(i + 1, grid); }
inline int previous_point(int i, const Grid& grid) { return periodic_index(i - 1, grid); }

// From the nodes to the cell centres and back; centred, so second order in dx.
ScalarField derivative_at_cell_centres(const ScalarField& on_nodes, const Grid& grid);
ScalarField derivative_at_nodes(const ScalarField& on_cell_centres, const Grid& grid);
ScalarField average_at_cell_centres(const ScalarField& on_nodes, const Grid& grid);
ScalarField average_at_nodes(const ScalarField& on_cell_centres, const Grid& grid);

// Order-1 shape function (slides 4-5): a particle between two grid points shares itself between
// them, in proportion to how close it is to each. The same weights deposit and gather.
struct ShapeWeights {
    int left_point, right_point;
    double left_weight, right_weight;
};

// distance_in_cells: how many cells the particle is to the right of point 0.
inline ShapeWeights shape_weights(double distance_in_cells, const Grid& grid) {
    // -1 left of the first cell centre; nx when x is just below L and the division rounds up.
    const int cells_passed = int(std::floor(distance_in_cells));
    const double distance_from_left_point = distance_in_cells - cells_passed;  // in [0, 1)

    const int left_point = periodic_index(cells_passed, grid);
    return {left_point, next_point(left_point, grid), 1.0 - distance_from_left_point, distance_from_left_point};
}

// x must be in [0, L).
inline ShapeWeights weights_on_nodes(double x, const Grid& grid) {
    return shape_weights(x / grid.dx, grid);
}

inline ShapeWeights weights_on_cell_centres(double x, const Grid& grid) {
    return shape_weights(x / grid.dx - 0.5, grid);
}

inline double interpolate(const ScalarField& values, const ShapeWeights& shape) {
    return shape.left_weight * values[shape.left_point] + shape.right_weight * values[shape.right_point];
}

inline void deposit_onto(ScalarField& moment, const ShapeWeights& shape, double amount) {
    moment[shape.left_point] += shape.left_weight * amount;
    moment[shape.right_point] += shape.right_weight * amount;
}
