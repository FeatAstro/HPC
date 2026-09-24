#include "interpolation.hpp"

#include <algorithm>
#include <cmath>

CellWeights bilinear_weights(const Grid2D& g, double x, double y) {
    const double cell_coordinate_x = (x - g.x0) / g.dx;
    const double cell_coordinate_y = (y - g.y0) / g.dy;

    const int i0 = static_cast<int>(std::floor(cell_coordinate_x));
    const int j0 = static_cast<int>(std::floor(cell_coordinate_y));
    const double fx = cell_coordinate_x - i0;
    const double fy = cell_coordinate_y - j0;

    CellWeights cw;
    cw.i0 = i0;
    cw.j0 = j0;
    cw.w00 = (1.0 - fx) * (1.0 - fy);
    cw.w10 = fx * (1.0 - fy);
    cw.w01 = (1.0 - fx) * fy;
    cw.w11 = fx * fy;
    return cw;
}

void deposit_moments(Grid2D& g, const std::vector<Particle>& particles) {
    std::fill(g.ion_density.begin(), g.ion_density.end(), 0.0);
    std::fill(g.ion_velocity.begin(), g.ion_velocity.end(), Vec3{});

    for (const Particle& p : particles) {
        const CellWeights cw = bilinear_weights(g, p.x, p.y);
        for_each_corner(g, cw, [&](int node, double weight) {
            g.ion_density[node] += weight * p.w;
            g.ion_velocity[node] += (weight * p.w) * p.v;
        });
    }

    // ion_velocity holds the summed weight*velocity; dividing by the summed weight
    // ion_density gives the weighted mean (bulk) velocity.
    for (std::size_t node = 0; node < g.ion_velocity.size(); ++node) {
        if (g.ion_density[node] > 0.0) {
            g.ion_velocity[node] = (1.0 / g.ion_density[node]) * g.ion_velocity[node];
        }
    }
}

Vec3 gather(const VectorField& field, const Grid2D& g, double x, double y) {
    const CellWeights cw = bilinear_weights(g, x, y);
    Vec3 result;
    for_each_corner(g, cw, [&](int node, double weight) { result += weight * field[node]; });
    return result;
}
