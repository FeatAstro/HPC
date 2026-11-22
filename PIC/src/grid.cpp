#include "grid.hpp"

// Cell centre i is between node i and node i + 1.
ScalarField derivative_at_cell_centres(const ScalarField& on_nodes, const Grid& grid) {
    ScalarField derivative(grid.nx);
    for (int i = 0; i < grid.nx; ++i) {
        derivative[i] = (on_nodes[next_point(i, grid)] - on_nodes[i]) / grid.dx;
    }
    return derivative;
}

// Node i is between cell centre i - 1 and cell centre i.
ScalarField derivative_at_nodes(const ScalarField& on_cell_centres, const Grid& grid) {
    ScalarField derivative(grid.nx);
    for (int i = 0; i < grid.nx; ++i) {
        derivative[i] = (on_cell_centres[i] - on_cell_centres[previous_point(i, grid)]) / grid.dx;
    }
    return derivative;
}

ScalarField average_at_cell_centres(const ScalarField& on_nodes, const Grid& grid) {
    ScalarField average(grid.nx);
    for (int i = 0; i < grid.nx; ++i) {
        average[i] = 0.5 * (on_nodes[i] + on_nodes[next_point(i, grid)]);
    }
    return average;
}

ScalarField average_at_nodes(const ScalarField& on_cell_centres, const Grid& grid) {
    ScalarField average(grid.nx);
    for (int i = 0; i < grid.nx; ++i) {
        average[i] = 0.5 * (on_cell_centres[previous_point(i, grid)] + on_cell_centres[i]);
    }
    return average;
}
