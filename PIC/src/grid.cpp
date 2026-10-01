#include "grid.hpp"

#include <cstddef>
#include <stdexcept>

Grid::Grid(int nx, double dx) : nx(nx), dx(dx) {
    if (nx < 2 || dx <= 0.0) {
        throw std::invalid_argument("Grid needs nx >= 2 and dx > 0");
    }
}

// Cell centre i lies between node i and node i + 1; node i lies between cell centres i - 1 and i.
// The neighbour beyond either end of the array is the one at the other end (periodic).

ScalarField derivative_at_cell_centres(const ScalarField& node_values, const Grid& grid) {
    const std::size_t n = node_values.size();
    ScalarField derivative(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t next = (i + 1 == n) ? 0 : i + 1;
        derivative[i] = (node_values[next] - node_values[i]) / grid.dx;
    }
    return derivative;
}

ScalarField derivative_at_nodes(const ScalarField& cell_centre_values, const Grid& grid) {
    const std::size_t n = cell_centre_values.size();
    ScalarField derivative(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t previous = (i == 0) ? n - 1 : i - 1;
        derivative[i] = (cell_centre_values[i] - cell_centre_values[previous]) / grid.dx;
    }
    return derivative;
}

ScalarField average_at_cell_centres(const ScalarField& node_values) {
    const std::size_t n = node_values.size();
    ScalarField average(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t next = (i + 1 == n) ? 0 : i + 1;
        average[i] = 0.5 * (node_values[i] + node_values[next]);
    }
    return average;
}

ScalarField average_at_nodes(const ScalarField& cell_centre_values) {
    const std::size_t n = cell_centre_values.size();
    ScalarField average(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t previous = (i == 0) ? n - 1 : i - 1;
        average[i] = 0.5 * (cell_centre_values[previous] + cell_centre_values[i]);
    }
    return average;
}
