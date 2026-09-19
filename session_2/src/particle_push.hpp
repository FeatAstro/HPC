#pragma once

#include <vector>

#include "grid.hpp"
#include "particle.hpp"

// Moves every particle at its current velocity for `duration`, wrapping around the periodic domain.
void drift_particles(std::vector<Particle>& particles, const Grid2D& g, double duration);

// Boris-pushes every particle's velocity over `dt` with the E and B gathered at its position.
void accelerate_particles(std::vector<Particle>& particles, const Grid2D& g, double dt);
