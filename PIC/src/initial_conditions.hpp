#pragma once

#include "field_solver.hpp"
#include "grid.hpp"
#include "population.hpp"

// Everything a run needs at t = 0.
struct Setup {
    Grid grid;
    Population ions;
    Fields fields;
    PlasmaParameters parameters;
    double dt;
    int n_steps;
    int steps_between_outputs;
};

// The simplest case: cold protons at rest, density 1, in a uniform magnetic field B = (1, 0, 0).
// No current, no pressure, no motion, so E = 0: nothing should change.
Setup uniform_plasma_at_rest();

// An ion sound wave: cold ions along B = (1, 0, 0), warm electrons (T_e = 1), and a density ripple
// n = 1 + delta cos kx moving with the ions. The electron pressure pushes the ions out of the dense
// regions, and the ripple travels at the sound speed sqrt(T_e): omega = k sqrt(T_e).
Setup ion_sound_wave();

// Resistive decay: no Bx, ions at rest, and a field By = b (sin x + sin 2x). The field does not move,
// it fades: each mode k decays as exp(-gamma t), gamma = eta k^2 + nu k^4.
Setup resistive_decay();

// A case only to measure the speed of the code: many warm protons (T = 1) in B = (1, 0, 0). The physics
// is not checked; the ions move and cross cells, as in a real run.
Setup timing_run();
