#pragma once

#include <array>

// Where each vector component lives relative to node (i, j), in units of the cell size:
// component c of entry (i, j) sits at (x0 + (i + offset.x) * dx, y0 + (j + offset.y) * dy).
struct CellOffset {
    double x = 0.0;
    double y = 0.0;
};

using Staggering = std::array<CellOffset, 3>;

// Yee layout for 2D3V: E and j on cell edges, B on cell faces.
namespace yee {

inline constexpr Staggering nodes = {{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}};
inline constexpr Staggering electric_field = {{{0.5, 0.0}, {0.0, 0.5}, {0.0, 0.0}}};
inline constexpr Staggering current_density = electric_field;
inline constexpr Staggering magnetic_field = {{{0.0, 0.5}, {0.5, 0.0}, {0.5, 0.5}}};

} // namespace yee
