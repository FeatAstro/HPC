#include "interpolation.hpp"

#include <algorithm>
#include <cmath>

CellWeights bilinear_weights(const Grid2D& g, double x, double y) {
    int i0 = static_cast<int>(std::floor((x - g.x0) / g.dx));
    int j0 = static_cast<int>(std::floor((y - g.y0) / g.dy));
    i0 = std::clamp(i0, 0, g.nx - 2);
    j0 = std::clamp(j0, 0, g.ny - 2);

    const double fx = (x - g.node_x(i0)) / g.dx;
    const double fy = (y - g.node_y(j0)) / g.dy;

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
    std::fill(g.n.begin(), g.n.end(), 0.0);
    std::fill(g.v.begin(), g.v.end(), Vec3{});

    for (const Particle& p : particles) {
        const CellWeights cw = bilinear_weights(g, p.x, p.y);
        const int idx00 = g.index(cw.i0, cw.j0);
        const int idx10 = g.index(cw.i0 + 1, cw.j0);
        const int idx01 = g.index(cw.i0, cw.j0 + 1);
        const int idx11 = g.index(cw.i0 + 1, cw.j0 + 1);

        g.n[idx00] += cw.w00 * p.w;
        g.n[idx10] += cw.w10 * p.w;
        g.n[idx01] += cw.w01 * p.w;
        g.n[idx11] += cw.w11 * p.w;

        g.v[idx00] += (cw.w00 * p.w) * p.v;
        g.v[idx10] += (cw.w10 * p.w) * p.v;
        g.v[idx01] += (cw.w01 * p.w) * p.v;
        g.v[idx11] += (cw.w11 * p.w) * p.v;
    }

    // At this point g.n holds the summed weight per node and g.v the summed
    // weight*velocity. The bulk velocity divides by the summed weight; only
    // afterwards is g.n converted to a density (weight per unit area).
    const double inv_cell_area = 1.0 / (g.dx * g.dy);
    for (std::size_t k = 0; k < g.v.size(); ++k) {
        if (g.n[k] > 0.0) {
            g.v[k] = (1.0 / g.n[k]) * g.v[k];
        }
        g.n[k] *= inv_cell_area;
    }
}

Vec3 gather(const std::vector<Vec3>& field, const Grid2D& g, double x, double y) {
    const CellWeights cw = bilinear_weights(g, x, y);
    Vec3 result;
    result += cw.w00 * field[g.index(cw.i0, cw.j0)];
    result += cw.w10 * field[g.index(cw.i0 + 1, cw.j0)];
    result += cw.w01 * field[g.index(cw.i0, cw.j0 + 1)];
    result += cw.w11 * field[g.index(cw.i0 + 1, cw.j0 + 1)];
    return result;
}
