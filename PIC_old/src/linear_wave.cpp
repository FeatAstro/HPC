#include "linear_wave.hpp"

#include <cmath>
#include <complex>
#include <numbers>

#include "diagnostics.hpp"
#include "pic_loop.hpp"

namespace {

// Mean of (b_y + i b_z) exp(-i k x): the complex amplitude of the wave.
std::complex<double> measure_wave_amplitude(const Grid2D& g, double wave_number) {
    const std::complex<double> imaginary_unit(0.0, 1.0);
    std::complex<double> sum = 0.0;
    g.for_each_node([&](int i, int, int node) {
        const double x_of_by = g.node_x(i, yee::magnetic_field[axis::y].x);
        const double x_of_bz = g.node_x(i, yee::magnetic_field[axis::z].x);
        sum += g.magnetic_field[node].y * std::polar(1.0, -wave_number * x_of_by);
        sum += imaginary_unit * g.magnetic_field[node].z * std::polar(1.0, -wave_number * x_of_bz);
    });
    return sum / static_cast<double>(g.nx * g.ny);
}

double background_magnetic_energy(const Grid2D& g, const PlasmaParameters& plasma, double background_field) {
    return 0.5 * background_field * background_field / plasma.vacuum_permeability * g.length_x() * g.length_y();
}

}

double parallel_wave_frequency(const PlasmaParameters& plasma, double background_field, double wave_number,
                               WaveBranch branch) {
    const double ion_cyclotron_frequency = plasma.elementary_charge * background_field / plasma.ion_mass;
    const double alfven_speed =
        background_field / std::sqrt(plasma.vacuum_permeability * plasma.reference_density * plasma.ion_mass);
    const double wave_number_in_ion_lengths = wave_number * alfven_speed / ion_cyclotron_frequency;

    const double squared = wave_number_in_ion_lengths * wave_number_in_ion_lengths;
    const double root = wave_number_in_ion_lengths * std::sqrt(squared + 4.0);
    const double normalized_frequency = branch == WaveBranch::ion_cyclotron ? 0.5 * (root - squared)
                                                                            : -0.5 * (root + squared);
    return normalized_frequency * ion_cyclotron_frequency;
}

WaveRunResult run_parallel_wave(const WaveRunSettings& settings) {
    const PlasmaParameters& plasma = settings.plasma;
    const ParallelWave& wave = settings.wave;

    Grid2D grid(settings.nx, settings.ny, settings.dx, settings.dy);
    std::vector<Particle> ions = load_uniform_ions(grid, plasma, settings.loading);

    const double wave_number = 2.0 * std::numbers::pi * wave.mode_number / grid.length_x();
    const double frequency = parallel_wave_frequency(plasma, settings.background_field, wave_number, wave.branch);

    grid.fill_field(grid.magnetic_field, yee::magnetic_field, [&](double x, double) {
        return Vec3{settings.background_field, wave.amplitude * std::cos(wave_number * x),
                    wave.amplitude * std::sin(wave_number * x)};
    });

    // Cold ions move with the wave: u = -(e/m) w / (k (Omega_i - w)) * b.
    const double ion_cyclotron_frequency = plasma.elementary_charge * settings.background_field / plasma.ion_mass;
    const double velocity_per_field = -(plasma.elementary_charge / plasma.ion_mass) * frequency /
                                      (wave_number * (ion_cyclotron_frequency - frequency));
    for (Particle& ion : ions) {
        ion.v = Vec3{0.0, velocity_per_field * wave.amplitude * std::cos(wave_number * ion.x),
                     velocity_per_field * wave.amplitude * std::sin(wave_number * ion.x)};
    }

    const double initial_energy = ion_kinetic_energy(grid, ions) + magnetic_energy(grid, plasma);
    const double wave_energy = initial_energy - background_magnetic_energy(grid, plasma, settings.background_field);

    const double wave_period = 2.0 * std::numbers::pi / std::fabs(frequency);
    const int number_of_steps =
        static_cast<int>(std::lround(settings.number_of_wave_periods * wave_period / settings.time_step));

    const std::complex<double> initial_amplitude = measure_wave_amplitude(grid, wave_number);
    std::complex<double> amplitude = initial_amplitude;
    double accumulated_phase = 0.0;

    for (int step = 0; step < number_of_steps; ++step) {
        advance_hybrid_step(grid, ions, plasma, settings.time_step);
        const std::complex<double> new_amplitude = measure_wave_amplitude(grid, wave_number);
        accumulated_phase += std::arg(new_amplitude / amplitude);
        amplitude = new_amplitude;
    }

    const double final_energy = ion_kinetic_energy(grid, ions) + magnetic_energy(grid, plasma);
    const double duration = number_of_steps * settings.time_step;

    WaveRunResult result;
    result.theoretical_frequency = frequency;
    result.measured_frequency = -accumulated_phase / duration;
    result.amplitude_ratio = std::abs(amplitude) / std::abs(initial_amplitude);
    result.relative_energy_change = (final_energy - initial_energy) / wave_energy;
    return result;
}
