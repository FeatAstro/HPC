#pragma once

#include <cmath>
#include <vector>

#include "vec3.hpp"

// Yee layout in 1D: a quantity lives either on the nodes, x_i = i dx, or at the cell centres,
// x_i = (i + 1/2) dx. A derivative taken from one set of points is centred on the other.
enum class Centering { node, cell_centre };

// Periodic 1D grid: nx cells of size dx, node nx is node 0 again, so the domain is [0, nx dx).
struct Grid {
    int nx;
    double dx;

    Grid(int nx, double dx);

    double length() const { return nx * dx; }
    double node_position(int i) const { return i * dx; }
    double cell_centre_position(int i) const { return (i + 0.5) * dx; }
};

// One value per node (or per cell centre).
using ScalarField = std::vector<double>;

// One value per node (or per cell centre) for each component.
struct VectorField {
    ScalarField x, y, z;
    Centering x_centering = Centering::node;
    Centering y_centering = Centering::node;
    Centering z_centering = Centering::node;

    VectorField() = default;
    explicit VectorField(const Grid& grid) : x(grid.nx, 0.0), y(grid.nx, 0.0), z(grid.nx, 0.0) {}
    VectorField(const Grid& grid, Centering x_centering, Centering y_centering, Centering z_centering)
        : x(grid.nx, 0.0), y(grid.nx, 0.0), z(grid.nx, 0.0),
          x_centering(x_centering), y_centering(y_centering), z_centering(z_centering) {}
};

// Yee layout: Bx, Ey, Ez on the nodes (with the moments); By, Bz, Ex at the cell centres.
// The current j = curl B lives where E does.
inline VectorField electric_field_layout(const Grid& grid) {
    return VectorField(grid, Centering::cell_centre, Centering::node, Centering::node);
}

inline VectorField magnetic_field_layout(const Grid& grid) {
    return VectorField(grid, Centering::node, Centering::cell_centre, Centering::cell_centre);
}

struct Fields {
    VectorField E, B;

    explicit Fields(const Grid& grid) : E(electric_field_layout(grid)), B(magnetic_field_layout(grid)) {}
};

// Centred differences and averages between the two sets of points, second order in dx.
ScalarField derivative_at_cell_centres(const ScalarField& node_values, const Grid& grid);
ScalarField derivative_at_nodes(const ScalarField& cell_centre_values, const Grid& grid);
ScalarField average_at_cell_centres(const ScalarField& node_values);
ScalarField average_at_nodes(const ScalarField& cell_centre_values);

// Order-1 shape function (slides 4-5): a particle between two points shares itself between them,
// linearly in its distance to each. The same weights are used to deposit and to gather.
struct ShapeWeights {
    int left_point, right_point;
    double left_weight, right_weight;
};

// x must be in [0, grid.length()).
inline ShapeWeights shape_weights(const Grid& grid, double x, Centering centering = Centering::node) {
    const double half_cell_shift = (centering == Centering::cell_centre) ? 0.5 : 0.0;
    const double position_in_cells = x / grid.dx - half_cell_shift;
    int left_point = static_cast<int>(std::floor(position_in_cells));
    const double right_weight = position_in_cells - left_point;
    if (left_point < 0) left_point += grid.nx;         // left of the first cell centre
    if (left_point >= grid.nx) left_point -= grid.nx;  // x just below L can round to exactly nx cells
    const int right_point = (left_point + 1 == grid.nx) ? 0 : left_point + 1;
    return {left_point, right_point, 1.0 - right_weight, right_weight};
}

// The weights of one particle for both centerings.
struct ParticleShape {
    ShapeWeights on_nodes, on_cell_centres;

    const ShapeWeights& on(Centering centering) const {
        return (centering == Centering::node) ? on_nodes : on_cell_centres;
    }
};

inline ParticleShape particle_shape(const Grid& grid, double x) {
    return {shape_weights(grid, x, Centering::node), shape_weights(grid, x, Centering::cell_centre)};
}

inline double interpolate(const ScalarField& values, const ShapeWeights& shape) {
    return shape.left_weight * values[shape.left_point] + shape.right_weight * values[shape.right_point];
}

inline Vec3 gather(const VectorField& field, const ParticleShape& shape) {
    return {interpolate(field.x, shape.on(field.x_centering)),
            interpolate(field.y, shape.on(field.y_centering)),
            interpolate(field.z, shape.on(field.z_centering))};
}

inline void deposit_onto(ScalarField& moment, const ShapeWeights& shape, double amount) {
    moment[shape.left_point] += shape.left_weight * amount;
    moment[shape.right_point] += shape.right_weight * amount;
}
