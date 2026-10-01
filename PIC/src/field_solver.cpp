#include "field_solver.hpp"

#include <stdexcept>

ScalarField electron_pressure(const ScalarField& density, const PlasmaParameters& parameters) {
    ScalarField pressure(density.size());
    for (std::size_t i = 0; i < density.size(); ++i) {
        pressure[i] = density[i] * parameters.electron_temperature;
    }
    return pressure;
}

VectorField ampere(const VectorField& B, const Grid& grid) {
    const ScalarField dBy_dx = derivative_at_nodes(B.y, grid);
    const ScalarField dBz_dx = derivative_at_nodes(B.z, grid);

    VectorField j = electric_field_layout(grid);
    for (int i = 0; i < grid.nx; ++i) {
        j.y[i] = -dBz_dx[i];
        j.z[i] = dBy_dx[i];
    }
    return j;
}

// Ohm's law divides by the density. A production code would replace an empty node by a small
// density floor; here an empty node is treated as a mistake in the setup.
static void check_density_is_positive(const ScalarField& density) {
    for (double n : density) {
        if (n <= 0.0) {
            throw std::invalid_argument("Ohm's law needs a positive density on every node");
        }
    }
}

static ScalarField second_derivative_at_nodes(const ScalarField& node_values, const Grid& grid) {
    return derivative_at_nodes(derivative_at_cell_centres(node_values, grid), grid);
}

VectorField ohm(const ScalarField& density, const VectorField& bulk_velocity, const VectorField& B,
                const VectorField& j, const Grid& grid, const PlasmaParameters& parameters) {
    check_density_is_positive(density);
    const VectorField& u = bulk_velocity;
    const ScalarField& n = density;
    const double eta = parameters.resistivity;
    const double nu = parameters.hyper_resistivity;

    // Ex lives at the cell centres, with By and Bz: the node quantities are averaged there.
    const ScalarField n_centre = average_at_cell_centres(n);
    const ScalarField uy_centre = average_at_cell_centres(u.y);
    const ScalarField uz_centre = average_at_cell_centres(u.z);
    const ScalarField jy_centre = average_at_cell_centres(j.y);
    const ScalarField jz_centre = average_at_cell_centres(j.z);
    const ScalarField dpe_dx = derivative_at_cell_centres(electron_pressure(n, parameters), grid);

    // Ey and Ez live on the nodes, with the moments and Bx: By and Bz are averaged there.
    const ScalarField By_node = average_at_nodes(B.y);
    const ScalarField Bz_node = average_at_nodes(B.z);
    const ScalarField laplacian_jy = second_derivative_at_nodes(j.y, grid);
    const ScalarField laplacian_jz = second_derivative_at_nodes(j.z, grid);

    VectorField E = electric_field_layout(grid);
    for (int i = 0; i < grid.nx; ++i) {
        E.x[i] = -(uy_centre[i] * B.z[i] - uz_centre[i] * B.y[i])
                 + (jy_centre[i] * B.z[i] - jz_centre[i] * B.y[i]) / n_centre[i]
                 - dpe_dx[i] / n_centre[i];

        E.y[i] = -(u.z[i] * B.x[i] - u.x[i] * Bz_node[i])
                 + j.z[i] * B.x[i] / n[i]
                 + eta * j.y[i] - nu * laplacian_jy[i];

        E.z[i] = -(u.x[i] * By_node[i] - u.y[i] * B.x[i])
                 - j.y[i] * B.x[i] / n[i]
                 + eta * j.z[i] - nu * laplacian_jz[i];
    }
    return E;
}

VectorField faraday(const VectorField& B, const VectorField& E, const Grid& grid, double dt) {
    const ScalarField dEy_dx = derivative_at_cell_centres(E.y, grid);
    const ScalarField dEz_dx = derivative_at_cell_centres(E.z, grid);

    VectorField advanced = B;
    for (int i = 0; i < grid.nx; ++i) {
        advanced.y[i] += dt * dEz_dx[i];
        advanced.z[i] -= dt * dEy_dx[i];
    }
    return advanced;
}

static VectorField average(const VectorField& a, const VectorField& b) {
    VectorField mean = a;
    for (std::size_t i = 0; i < a.x.size(); ++i) {
        mean.x[i] = 0.5 * (a.x[i] + b.x[i]);
        mean.y[i] = 0.5 * (a.y[i] + b.y[i]);
        mean.z[i] = 0.5 * (a.z[i] + b.z[i]);
    }
    return mean;
}

static VectorField electric_field(const ScalarField& density, const VectorField& bulk_velocity, const VectorField& B,
                                  const Grid& grid, const PlasmaParameters& parameters) {
    return ohm(density, bulk_velocity, B, ampere(B, grid), grid, parameters);
}

void advance_fields_with_fixed_ions(Fields& fields, const ScalarField& density, const VectorField& bulk_velocity,
                                    const Grid& grid, const PlasmaParameters& parameters, double dt) {
    const VectorField E_start = electric_field(density, bulk_velocity, fields.B, grid, parameters);

    VectorField E_mid_step = E_start;  // first guess: E does not change during the step
    for (int prediction = 1; prediction <= 2; ++prediction) {
        const VectorField B_predicted = faraday(fields.B, E_mid_step, grid, dt);
        const VectorField E_predicted = electric_field(density, bulk_velocity, B_predicted, grid, parameters);
        E_mid_step = average(E_start, E_predicted);
    }

    fields.B = faraday(fields.B, E_mid_step, grid, dt);
    fields.E = electric_field(density, bulk_velocity, fields.B, grid, parameters);
}
