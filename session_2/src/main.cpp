#include <cmath>
#include <iomanip>
#include <numbers>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include "grid.hpp"
#include "interpolation.hpp"
#include "linear_wave.hpp"
#include "particle.hpp"
#include "step1_single_particle.hpp"
#include "step2_pic_loop.hpp"

namespace {

void run_step1_demo() {
    std::cout << "--- Step 1: single-particle Boris pusher (uniform B) ---\n";

    const double B0 = 1.0;
    const auto E_of = [](double, double) { return Vec3{0.0, 0.0, 0.0}; };
    const auto B_of = [B0](double, double) { return Vec3{0.0, 0.0, B0}; };

    Particle p;
    p.q = 1.0;
    p.m = 1.0;
    p.v = Vec3{1.0, 0.0, 0.0};

    const double omega_c = p.q * B0 / p.m;
    const double period = 2.0 * std::numbers::pi / omega_c;
    const int steps_per_period = 1000;
    const double dt = period / steps_per_period;

    const double speed0 = norm(p.v);
    for (int step = 0; step < steps_per_period; ++step) {
        advance_particle(p, E_of, B_of, dt);
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  initial speed = " << speed0 << ", final speed = " << norm(p.v) << "\n";
    std::cout << "  final position after one gyration (should be near start): (" << p.x << ", " << p.y
               << ")\n\n";
}

void run_step2_demo() {
    std::cout << "--- Step 2: N particles depositing to / gathering from mesh ---\n";

    Grid2D grid(21, 21, 1.0, 1.0, 0.0, 0.0);
    for (Vec3& B : grid.magnetic_field) {
        B = Vec3{0.0, 0.0, 1.0};
    }

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> pos_dist(2.0, 18.0);
    std::uniform_real_distribution<double> vel_dist(-0.5, 0.5);

    const int N = 200;
    std::vector<Particle> particles(N);
    double total_weight = 0.0;
    for (Particle& p : particles) {
        p.x = pos_dist(rng);
        p.y = pos_dist(rng);
        p.v = Vec3{vel_dist(rng), vel_dist(rng), 0.0};
        p.q = 1.0;
        p.m = 1.0;
        p.w = 1.0;
        total_weight += p.w;
    }

    const double dt = 0.01;
    for (int step = 0; step < 100; ++step) {
        advance_particles(particles, grid, dt);
    }

    deposit_moments(grid, particles);
    double deposited_weight = 0.0;
    for (double n_ij : grid.ion_density) {
        deposited_weight += n_ij;
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  total particle weight   = " << total_weight << "\n";
    std::cout << "  weight recovered on grid = " << deposited_weight << "\n\n";
}

void run_hybrid_wave_demo() {
    std::cout << "--- Hybrid PIC: parallel waves (kinetic ions, fluid electrons) ---\n";

    const std::pair<WaveBranch, const char*> branches[] = {{WaveBranch::ion_cyclotron, "ion-cyclotron"},
                                                           {WaveBranch::whistler, "whistler"}};
    for (const auto& [branch, name] : branches) {
        WaveRunSettings settings;
        settings.wave.branch = branch;
        const WaveRunResult result = run_parallel_wave(settings);

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  " << name << " branch: theoretical frequency = " << result.theoretical_frequency
                  << ", measured = " << result.measured_frequency << ", amplitude after one period = "
                  << result.amplitude_ratio << " x initial, energy change = " << std::scientific
                  << std::setprecision(2) << result.relative_energy_change << " x wave energy\n";
    }
    std::cout << "\n";
}

} // namespace

int main() {
    std::cout << "kinetic_fisher: session 2 - hybrid-kinetic PIC (kinetic ions, fluid electrons)\n\n";
    run_step1_demo();
    run_step2_demo();
    run_hybrid_wave_demo();
    return 0;
}
