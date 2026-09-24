#pragma once

#include <functional>

#include "particle.hpp"
#include "vec3.hpp"

using FieldFunc = std::function<Vec3(double x, double y)>;

// Leapfrog step (slide 7): r^{n+1/2} = r^n + v^n dt/2, Boris push with E, B at r^{n+1/2},
// r^{n+1} = r^{n+1/2} + v^{n+1} dt/2.
void advance_particle(Particle& p, const FieldFunc& E_of, const FieldFunc& B_of, double dt);
