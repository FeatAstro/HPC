#include "pic_loop.hpp"

#include "field_solver.hpp"
#include "interpolation.hpp"
#include "particle_push.hpp"

namespace {

// Advances B from magnetic_field_start over dt with the moments already on the grid.
// Leaves the half-step B and the E computed from it on the grid; returns B at the end of the step.
VectorField solve_fields_for_step(Grid2D& g, const PlasmaParameters& plasma,
                                  const VectorField& magnetic_field_start, double dt) {
    g.magnetic_field = magnetic_field_start;
    advance_magnetic_field(g, plasma, dt);
    const VectorField magnetic_field_end = g.magnetic_field;

    g.magnetic_field = averaged(magnetic_field_start, magnetic_field_end);
    update_electric_field(g, plasma);
    return magnetic_field_end;
}

} // namespace

void advance_hybrid_step(Grid2D& g, std::vector<Particle>& ions, const PlasmaParameters& plasma, double dt) {
    const VectorField magnetic_field_start = g.magnetic_field;

    drift_particles(ions, g, dt / 2.0);

    // Predictor: moments from the velocities at t^n give a first estimate of v^{n+1}.
    deposit_moments(g, ions);
    solve_fields_for_step(g, plasma, magnetic_field_start, dt);
    std::vector<Particle> predicted_ions = ions;
    accelerate_particles(predicted_ions, g, dt);

    // Corrector: moments from the half-step velocities.
    std::vector<Particle> half_step_ions = ions;
    for (std::size_t k = 0; k < ions.size(); ++k) {
        half_step_ions[k].v = 0.5 * (ions[k].v + predicted_ions[k].v);
    }
    deposit_moments(g, half_step_ions);
    const VectorField magnetic_field_end = solve_fields_for_step(g, plasma, magnetic_field_start, dt);

    accelerate_particles(ions, g, dt);
    drift_particles(ions, g, dt / 2.0);
    g.magnetic_field = magnetic_field_end;
}
