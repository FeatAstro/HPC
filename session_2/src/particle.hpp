#pragma once

#include "vec3.hpp"

struct Particle {
    double x = 0.0;
    double y = 0.0;
    Vec3 v;
    double q = 1.0;
    double m = 1.0;
    double w = 1.0; // statistical weight w_p
};
