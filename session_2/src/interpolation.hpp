#pragma once

#include <vector>

#include "grid.hpp"
#include "particle.hpp"
#include "vec3.hpp"

// First-order (bilinear) shape function, slide 5: for a particle between nodes
// (i0,j0) and (i0+1,j0+1), the weight on each corner is the product of the two
// linear 1D weights (x_p-x_i)/dx and (y_p-y_j)/dy. Order 1 = exact linear interpolation.
struct CellWeights {
    int i0, j0;
    double w00, w10, w01, w11; // weights at (i0,j0),(i0+1,j0),(i0,j0+1),(i0+1,j0+1)
};

CellWeights bilinear_weights(const Grid2D& g, double x, double y);

// Deposit: accumulate particle weight/velocity onto the grid moments n_ij, v_ij
// (slide 4/5: n_ij = sum_p S(r_p - r_ij) w_p). Following the slides, w_p is
// already a density contribution, so no division by the cell area is done here
// and sum(n_ij) equals sum(w_p) for any grid spacing.
void deposit_moments(Grid2D& g, const std::vector<Particle>& particles);

// Gather: interpolate a node-centered field back to a particle position, using
// the same shape function as deposit (slide 6: Q(r) = sum Q_ij S(r_p - r_ij)).
Vec3 gather(const std::vector<Vec3>& field, const Grid2D& g, double x, double y);
