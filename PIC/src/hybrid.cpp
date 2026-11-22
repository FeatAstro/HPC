#include "hybrid.hpp"

#include "stopwatch.hpp"

void advance_one_step(Population& ions, Fields& fields, const Grid& grid, const PlasmaParameters& parameters,
                      double dt, StepTimes& times) {
    Stopwatch stopwatch;
    const Population ions_at_start = ions;
    times.copies += stopwatch.lap_seconds();
    Moments moments = deposit(ions, grid);
    times.particles += stopwatch.lap_seconds();

    const Fields start = fields;
    Fields mid_step = start;  // first guess: the fields do not change during the step
    for (int prediction = 1; prediction <= 2; ++prediction) {
        const VectorField B_at_end = faraday(start.B, mid_step.E, grid, dt);
        const VectorField E_at_end = ohm(moments, B_at_end, grid, parameters);
        mid_step = {average(start.E, E_at_end), average(start.B, B_at_end)};
        times.fields += stopwatch.lap_seconds();

        ions = ions_at_start;
        times.copies += stopwatch.lap_seconds();
        push(ions, mid_step, grid, dt);
        moments = deposit(ions, grid);
        times.particles += stopwatch.lap_seconds();
    }

    // correction: the same fields at the end, from the second prediction
    fields.B = faraday(start.B, mid_step.E, grid, dt);
    fields.E = ohm(moments, fields.B, grid, parameters);
    times.fields += stopwatch.lap_seconds();
}
