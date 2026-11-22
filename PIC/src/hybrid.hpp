#pragma once

#include "field_solver.hpp"
#include "grid.hpp"
#include "population.hpp"

// Where the time of the run goes: push and deposit are split between the threads (OpenMP), the rest is not.
struct StepTimes {
    double particles = 0.0;  // push, deposit: in parallel
    double copies = 0.0;     // copies of the ions, to restart each prediction from the start of the step
    double fields = 0.0;     // Faraday, Ohm (with Ampere), averages
};

// One time step of the ions and the fields together, by the iterated Crank-Nicolson scheme
// (session 3, slide 7): two predictions, then the correction. Each prediction advances B, deduces E,
// pushes the ions from their state at the start of the step in the fields averaged over the step,
// and deposits their moments. The ions keep the result of the second prediction.
// Needs at the start: fields.E from Ohm's law. Adds the time spent to `times`.
void advance_one_step(Population& ions, Fields& fields, const Grid& grid, const PlasmaParameters& parameters,
                      double dt, StepTimes& times);
