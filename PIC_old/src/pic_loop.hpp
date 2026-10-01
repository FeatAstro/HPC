#pragma once

#include <vector>

#include "grid.hpp"
#include "particle.hpp"
#include "plasma_parameters.hpp"

// One self-consistent hybrid step (slide 6): kinetic ions from t^n to t^{n+1},
// fluid electrons through Ohm's law, B through Faraday's law.
//
// On entry: ion positions and velocities at t^n, grid.magnetic_field = B^n.
// On exit:  ions at t^{n+1}, grid.magnetic_field = B^{n+1},
//           grid.electric_field = E^{n+1/2}, grid.ion_density/ion_velocity = moments at t^{n+1/2}.
//
// A predictor pass estimates the half-step ion velocities, which the corrector pass
// uses to deposit the moments that drive the fields the ions are finally pushed with.
void advance_hybrid_step(Grid2D& g, std::vector<Particle>& ions, const PlasmaParameters& plasma, double dt);
