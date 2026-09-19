#pragma once

#include <vector>

#include "grid.hpp"
#include "particle.hpp"

// One PIC-loop step for a set of particles against a grid with prescribed
// (static) E, B fields (slide 6, minus the field-equations box, which is not
// implemented yet): for each particle, gather E,B at its mid-step position,
// Boris-push it, then deposit moments for diagnostics.
void advance_particles(std::vector<Particle>& particles, const Grid2D& g, double dt);
