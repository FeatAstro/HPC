// Full hybrid PIC loop tests: a uniform plasma in a uniform field stays in equilibrium, and
// parallel circularly polarised waves propagate at the frequency given by the exact
// dispersion relation (kinetic ions, fluid electrons).
#include <cmath>
#include <iostream>

#include "check.hpp"
#include "diagnostics.hpp"
#include "field_solver.hpp"
#include "linear_wave.hpp"
#include "particle_loading.hpp"
#include "pic_loop.hpp"

namespace {

void test_uniform_drift_is_an_equilibrium() {
    Grid2D grid(16, 8, 0.5, 0.5);
    PlasmaParameters plasma;
    plasma.electron_temperature = 0.4;

    const Vec3 background_field{0.3, 0.4, 0.5};
    const Vec3 drift_velocity{0.2, -0.1, 0.3};
    std::fill(grid.magnetic_field.begin(), grid.magnetic_field.end(), background_field);

    std::vector<Particle> ions = load_uniform_ions(grid, plasma, IonLoading{});
    const std::vector<Particle> initial_ions = [&] {
        for (Particle& ion : ions) {
            ion.v = drift_velocity;
        }
        return ions;
    }();

    const double time_step = 0.05;
    const int number_of_steps = 50;
    for (int step = 0; step < number_of_steps; ++step) {
        advance_hybrid_step(grid, ions, plasma, time_step);
    }

    // E = -v x B cancels the magnetic force, so ions keep drifting in a straight line.
    const double elapsed = number_of_steps * time_step;
    for (std::size_t k = 0; k < ions.size(); ++k) {
        check_close(ions[k].v.x, drift_velocity.x, 1e-9, "drift velocity x");
        check_close(ions[k].v.y, drift_velocity.y, 1e-9, "drift velocity y");
        check_close(ions[k].v.z, drift_velocity.z, 1e-9, "drift velocity z");
        check_close(ions[k].x, grid.wrap_x(initial_ions[k].x + drift_velocity.x * elapsed), 1e-9, "drifted x");
        check_close(ions[k].y, grid.wrap_y(initial_ions[k].y + drift_velocity.y * elapsed), 1e-9, "drifted y");
    }
    for (const Vec3& B : grid.magnetic_field) {
        check_close(B.x, background_field.x, 1e-10, "B unchanged, x");
        check_close(B.y, background_field.y, 1e-10, "B unchanged, y");
        check_close(B.z, background_field.z, 1e-10, "B unchanged, z");
    }
    std::cout << "test_uniform_drift_is_an_equilibrium: passed\n";
}

void check_wave_matches_dispersion_relation(WaveBranch branch, const char* what) {
    WaveRunSettings settings;
    settings.wave.branch = branch;

    const WaveRunResult result = run_parallel_wave(settings);

    std::cout << "  " << what << ": theory " << result.theoretical_frequency << ", measured "
              << result.measured_frequency << ", amplitude ratio " << result.amplitude_ratio << ", energy change "
              << result.relative_energy_change << "\n";

    check_close(result.measured_frequency / result.theoretical_frequency, 1.0, 0.02, what);
    check_close(result.amplitude_ratio, 1.0, 0.03, what);
    check_close(result.relative_energy_change, 0.0, 0.02, what);
}

void test_waves_follow_the_dispersion_relation() {
    check_wave_matches_dispersion_relation(WaveBranch::ion_cyclotron, "ion-cyclotron branch");
    check_wave_matches_dispersion_relation(WaveBranch::whistler, "whistler branch");
    std::cout << "test_waves_follow_the_dispersion_relation: passed\n";
}

}

int main() {
    test_uniform_drift_is_an_equilibrium();
    test_waves_follow_the_dispersion_relation();
    std::cout << "test_hybrid_waves: all checks passed\n";
    return 0;
}
