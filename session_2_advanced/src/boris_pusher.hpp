#pragma once

#include "particle.hpp"
#include "vec3.hpp"

// Rotates/accelerates p.v in place using the Boris scheme (slide 7):
//   v- = v + (q dt / 2m) E
//   t  = (q dt / 2m) B ;  v' = v- + v- x t ;  s = 2t / (1 + t.t) ;  v+ = v- + v' x s
//   v_new = v+ + (q dt / 2m) E
// E and B must already be the field values sampled at the particle's mid-step
// position r^{n+1/2}. Position updates are the caller's responsibility.
void boris_push(Particle& p, const Vec3& E, const Vec3& B, double dt);
