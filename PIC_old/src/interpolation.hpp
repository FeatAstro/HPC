#pragma once

#include <vector>

#include "grid.hpp"
#include "particle.hpp"
#include "vec3.hpp"
#include "yee_layout.hpp"

// First-order (bilinear) shape function, slide 5: for a particle between nodes
// (i0,j0) and (i0+1,j0+1), the weight on each corner is the product of the two
// linear 1D weights (x_p-x_i)/dx and (y_p-y_j)/dy. Order 1 = exact linear interpolation.
struct CellWeights {
    int i0, j0;
    double w00, w10, w01, w11; // weights at (i0,j0),(i0+1,j0),(i0,j0+1),(i0+1,j0+1)
};

// `offset` is where the quantity lives relative to the nodes (see yee_layout.hpp).
CellWeights bilinear_weights(const Grid2D& g, double x, double y, CellOffset offset = {});

// Calls visit(node_index, weight) for the four nodes around the particle.
template <typename Visit>
void for_each_corner(const Grid2D& g, const CellWeights& cw, Visit visit) {
    visit(g.index(cw.i0, cw.j0), cw.w00);
    visit(g.index(cw.i0 + 1, cw.j0), cw.w10);
    visit(g.index(cw.i0, cw.j0 + 1), cw.w01);
    visit(g.index(cw.i0 + 1, cw.j0 + 1), cw.w11);
}

template <typename ValueAtNode>
double weighted_sum(const Grid2D& g, const CellWeights& cw, ValueAtNode value_at_node) {
    double sum = 0.0;
    for_each_corner(g, cw, [&](int node, double weight) { sum += weight * value_at_node(node); });
    return sum;
}

// Deposit: accumulate particle weight/velocity onto the nodal moments ion_density,
// ion_velocity (slide 4/5: n_ij = sum_p S(r_p - r_ij) w_p). Following the slides, w_p is
// already a density contribution, so no division by the cell area is done here
// and sum(n_ij) equals sum(w_p) for any grid spacing.
void deposit_moments(Grid2D& g, const std::vector<Particle>& particles);

// Gather: interpolate a field back to a position, using the same shape function as
// deposit (slide 6: Q(r) = sum Q_ij S(r_p - r_ij)). Each component of a vector field
// is interpolated from where it lives, given by `staggering`.
Vec3 gather(const VectorField& field, const Staggering& staggering, const Grid2D& g, double x, double y);
double gather_scalar(const ScalarField& field, const Grid2D& g, double x, double y);
