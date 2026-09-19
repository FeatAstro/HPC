#pragma once

#include <vector>

#include "vec3.hpp"

// Simple node-centered 2D grid, nodes at (x0 + i*dx, y0 + j*dy), i in [0,nx), j in [0,ny).
// Not Yee-staggered: no curl operators are needed yet (no field solve at this stage),
// so E, B, and the moments all live on the same nodes.
struct Grid2D {
    int nx, ny;
    double dx, dy;
    double x0, y0;

    std::vector<double> n; // density moment, size nx*ny
    std::vector<Vec3> v;   // bulk velocity moment, size nx*ny
    std::vector<Vec3> E;   // prescribed field, set externally for now
    std::vector<Vec3> B;   // prescribed field, set externally for now

    Grid2D(int nx_, int ny_, double dx_, double dy_, double x0_ = 0.0, double y0_ = 0.0)
        : nx(nx_), ny(ny_), dx(dx_), dy(dy_), x0(x0_), y0(y0_),
          n(static_cast<std::size_t>(nx_) * ny_, 0.0),
          v(static_cast<std::size_t>(nx_) * ny_),
          E(static_cast<std::size_t>(nx_) * ny_),
          B(static_cast<std::size_t>(nx_) * ny_) {}

    int index(int i, int j) const { return j * nx + i; }

    double node_x(int i) const { return x0 + i * dx; }
    double node_y(int j) const { return y0 + j * dy; }
};
