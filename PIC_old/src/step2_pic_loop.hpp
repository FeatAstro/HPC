#pragma once

#include <vector>

#include "grid.hpp"
#include "particle.hpp"

// Leapfrog step against the E, B already on the grid (prescribed, not solved):
// drift dt/2, gather and Boris-push, drift dt/2.
void advance_particles(std::vector<Particle>& particles, const Grid2D& g, double dt);
