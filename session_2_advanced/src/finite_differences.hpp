#pragma once

#include "grid.hpp"

// Differences of any quantity stored at grid nodes, given as value_at_node(node_index).
// A forward difference at (i, j) lands half a cell to the right/above node (i, j);
// a backward difference lands half a cell to the left/below.

template <typename ValueAtNode>
double forward_difference_x(const Grid2D& g, ValueAtNode value_at_node, int i, int j) {
    return (value_at_node(g.index(i + 1, j)) - value_at_node(g.index(i, j))) / g.dx;
}

template <typename ValueAtNode>
double forward_difference_y(const Grid2D& g, ValueAtNode value_at_node, int i, int j) {
    return (value_at_node(g.index(i, j + 1)) - value_at_node(g.index(i, j))) / g.dy;
}

template <typename ValueAtNode>
double backward_difference_x(const Grid2D& g, ValueAtNode value_at_node, int i, int j) {
    return (value_at_node(g.index(i, j)) - value_at_node(g.index(i - 1, j))) / g.dx;
}

template <typename ValueAtNode>
double backward_difference_y(const Grid2D& g, ValueAtNode value_at_node, int i, int j) {
    return (value_at_node(g.index(i, j)) - value_at_node(g.index(i, j - 1))) / g.dy;
}

inline auto component_of(const VectorField& field, int component) {
    return [&field, component](int node) { return field[node][component]; };
}
