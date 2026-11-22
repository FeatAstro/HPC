#include "field_solver.hpp"

#include <cstddef>
#include <stdexcept>

VectorField ampere(const VectorField& B, const Grid& grid) {
    const ScalarField dBy_dx = derivative_at_nodes(B.y, grid);
    const ScalarField dBz_dx = derivative_at_nodes(B.z, grid);

    VectorField j(grid);
    for (int i = 0; i < grid.nx; ++i) {
        j.y[i] = -dBz_dx[i];
        j.z[i] = dBy_dx[i];
    }
    return j;
}

// Ohm's law divides by the density. A production code would replace an empty node by a small
// density floor; here an empty node is treated as a mistake in the setup.
// Written !(n > 0) rather than n <= 0 so that it also stops a run that has blown up: every
// comparison with NaN is false.
static void check_density_is_positive(const ScalarField& density) {
    for (double n : density) {
        if (!(n > 0.0)) {
            throw std::invalid_argument("Ohm's law needs a positive density on every node (is the run unstable?)");
        }
    }
}

static ScalarField laplacian_at_nodes(const ScalarField& on_nodes, const Grid& grid) {
    return derivative_at_nodes(derivative_at_cell_centres(on_nodes, grid), grid);
}

VectorField ohm(const Moments& moments, const VectorField& B, const Grid& grid, const PlasmaParameters& parameters) {
    const ScalarField& n = moments.density;
    const VectorField& u = moments.bulk_velocity;
    check_density_is_positive(n);
    const VectorField j = ampere(B, grid);
    const double T_e = parameters.electron_temperature;
    const double eta = parameters.resistivity;
    const double nu = parameters.hyper_resistivity;

    // Ex lives at the cell centres, with By and Bz: the node quantities are averaged there.
    const ScalarField n_centre = average_at_cell_centres(n, grid);
    const ScalarField uy_centre = average_at_cell_centres(u.y, grid);
    const ScalarField uz_centre = average_at_cell_centres(u.z, grid);
    const ScalarField jy_centre = average_at_cell_centres(j.y, grid);
    const ScalarField jz_centre = average_at_cell_centres(j.z, grid);
    const ScalarField dn_dx = derivative_at_cell_centres(n, grid);

    // Ey and Ez live on the nodes, with the moments and Bx: By and Bz are averaged there.
    const ScalarField By_node = average_at_nodes(B.y, grid);
    const ScalarField Bz_node = average_at_nodes(B.z, grid);
    const ScalarField laplacian_jy = laplacian_at_nodes(j.y, grid);
    const ScalarField laplacian_jz = laplacian_at_nodes(j.z, grid);

    // In 1D nothing varies along y or z and jx = 0: the pressure term is only in Ex,
    // the resistive terms only in Ey and Ez.
    VectorField E(grid);
    for (int i = 0; i < grid.nx; ++i) {
        E.x[i] = -(uy_centre[i] * B.z[i] - uz_centre[i] * B.y[i])                 // -u x B
                 + (jy_centre[i] * B.z[i] - jz_centre[i] * B.y[i]) / n_centre[i]  // + j x B / n
                 - T_e * dn_dx[i] / n_centre[i];                                  // - grad(p_e) / n

        E.y[i] = -(u.z[i] * B.x[i] - u.x[i] * Bz_node[i])  // -u x B
                 + j.z[i] * B.x[i] / n[i]                  // + j x B / n
                 + eta * j.y[i]                            // + eta j
                 - nu * laplacian_jy[i];                   // - nu laplacian(j)

        E.z[i] = -(u.x[i] * By_node[i] - u.y[i] * B.x[i])  // -u x B
                 - j.y[i] * B.x[i] / n[i]                  // + j x B / n
                 + eta * j.z[i]                            // + eta j
                 - nu * laplacian_jz[i];                   // - nu laplacian(j)
    }
    return E;
}

VectorField faraday(const VectorField& B, const VectorField& E, const Grid& grid, double dt) {
    const ScalarField dEy_dx = derivative_at_cell_centres(E.y, grid);
    const ScalarField dEz_dx = derivative_at_cell_centres(E.z, grid);

    VectorField B_after_dt = B;
    for (int i = 0; i < grid.nx; ++i) {
        B_after_dt.y[i] += dt * dEz_dx[i];
        B_after_dt.z[i] -= dt * dEy_dx[i];
    }
    return B_after_dt;
}

VectorField average(const VectorField& a, const VectorField& b) {
    VectorField mean = a;
    for (std::size_t i = 0; i < a.x.size(); ++i) {
        mean.x[i] = 0.5 * (a.x[i] + b.x[i]);
        mean.y[i] = 0.5 * (a.y[i] + b.y[i]);
        mean.z[i] = 0.5 * (a.z[i] + b.z[i]);
    }
    return mean;
}
