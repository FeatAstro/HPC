#pragma once

#include "particle_loading.hpp"
#include "plasma_parameters.hpp"

// Circularly polarised waves travelling along the background field B0 (along x), with cold
// ions and massless electrons. With b = b_y + i b_z varying as exp(i(kx - wt)), the exact
// dispersion relation is  w^2 = k^2 v_A^2 (1 - w/Omega_i), with two roots:
//   ion_cyclotron: w > 0, left-handed, w -> Omega_i for large k
//   whistler:      w < 0, right-handed, w grows like k^2
enum class WaveBranch { ion_cyclotron, whistler };

struct ParallelWave {
    int mode_number = 2;      // wavelengths that fit along x
    double amplitude = 0.01;  // magnitude of the transverse magnetic perturbation
    WaveBranch branch = WaveBranch::ion_cyclotron;
};

struct WaveRunSettings {
    int nx = 64;
    int ny = 4;
    double dx = 0.5;
    double dy = 0.5;
    PlasmaParameters plasma;
    IonLoading loading;
    double background_field = 1.0;
    ParallelWave wave;
    double time_step = 0.02;
    double number_of_wave_periods = 1.0;
};

struct WaveRunResult {
    double theoretical_frequency;
    double measured_frequency;
    double amplitude_ratio;         // final / initial wave amplitude
    double relative_energy_change;  // total energy change divided by the initial wave energy
};

// Signed frequency of the branch (negative for the whistler, by the convention above).
double parallel_wave_frequency(const PlasmaParameters& plasma, double background_field, double wave_number,
                               WaveBranch branch);

WaveRunResult run_parallel_wave(const WaveRunSettings& settings);
