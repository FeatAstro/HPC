#pragma once

#include <cmath>
#include <cstddef>

#include "vector_field.hpp"
#include "yee_layout.hpp"

// Periodic 2D grid with nx-by-ny nodes at (x0 + i*dx, y0 + j*dy). The domain is
// [x0, x0 + nx*dx) x [y0, y0 + ny*dy): node nx is node 0 again.
// Moments live on the nodes; E, B and j follow the Yee layout in yee_layout.hpp.
struct Grid2D {
    int nx, ny;
    double dx, dy;
    double x0, y0;

    ScalarField ion_density;
    VectorField ion_velocity;
    VectorField electric_field;
    VectorField magnetic_field;
    VectorField current_density;

    Grid2D(int nx_, int ny_, double dx_, double dy_, double x0_ = 0.0, double y0_ = 0.0)
        : nx(nx_), ny(ny_), dx(dx_), dy(dy_), x0(x0_), y0(y0_),
          ion_density(static_cast<std::size_t>(nx_) * ny_, 0.0),
          ion_velocity(static_cast<std::size_t>(nx_) * ny_),
          electric_field(static_cast<std::size_t>(nx_) * ny_),
          magnetic_field(static_cast<std::size_t>(nx_) * ny_),
          current_density(static_cast<std::size_t>(nx_) * ny_) {}

    // Indices wrap around, so neighbours of edge nodes are always valid.
    int index(int i, int j) const { return wrap_index(j, ny) * nx + wrap_index(i, nx); }

    double node_x(int i, double offset_in_cells = 0.0) const { return x0 + (i + offset_in_cells) * dx; }
    double node_y(int j, double offset_in_cells = 0.0) const { return y0 + (j + offset_in_cells) * dy; }

    double length_x() const { return nx * dx; }
    double length_y() const { return ny * dy; }

    double wrap_x(double x) const { return wrap_position(x, x0, length_x()); }
    double wrap_y(double y) const { return wrap_position(y, y0, length_y()); }

    template <typename Visit>
    void for_each_node(Visit visit) const {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                visit(i, j, index(i, j));
            }
        }
    }

    // Sets every component of `field` from field_at(x, y), evaluated where that component lives.
    template <typename FieldFunction>
    void fill_field(VectorField& field, const Staggering& staggering, FieldFunction field_at) const {
        for_each_node([&](int i, int j, int node) {
            for (int component = 0; component < axis::count; ++component) {
                const double x = node_x(i, staggering[component].x);
                const double y = node_y(j, staggering[component].y);
                field[node][component] = field_at(x, y)[component];
            }
        });
    }

private:
    static int wrap_index(int index, int size) { return ((index % size) + size) % size; }

    static double wrap_position(double position, double origin, double length) {
        return position - length * std::floor((position - origin) / length);
    }
};
