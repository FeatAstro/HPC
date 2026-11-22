// Runs one case of src/initial_conditions.hpp and writes its magnetic field and density to <case>.csv.
// Usage:  ./hybrid uniform   (or sound, decay, timing)
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <omp.h>

#include "field_solver.hpp"
#include "hybrid.hpp"
#include "initial_conditions.hpp"
#include "population.hpp"
#include "stopwatch.hpp"

// One line: t, then By and Bz at every cell centre, then the ion density at every node.
static void write_output(std::ofstream& file, double t, const VectorField& B, const ScalarField& density) {
    file << t;
    for (double By : B.y) file << ',' << By;
    for (double Bz : B.z) file << ',' << Bz;
    for (double n : density) file << ',' << n;
    file << '\n';
}

// The wall-clock time of the time loop (not the loading of the ions), and how it splits between
// the particles and the fields. Amdahl's law: if a fraction p of the time can be split between
// cores, the run can never be more than 1 / (1 - p) times faster, however many cores are used.
static void print_timing(double total_seconds, const StepTimes& times, int n_steps, std::size_t n_particles) {
    const double seconds_per_step = total_seconds / n_steps;
    const double rest = total_seconds - times.particles - times.copies - times.fields;
    const double parallel_fraction = times.particles / total_seconds;
    std::cout << std::setprecision(3)  // three significant digits
              << "  threads:               " << omp_get_max_threads() << '\n'
              << "  particles:             " << n_particles << '\n'
              << "  total time:            " << total_seconds << " s\n"
              << "  per step:              " << 1e3 * seconds_per_step << " ms\n"
              << "  per particle per step: " << 1e9 * seconds_per_step / n_particles << " ns\n"
              << "  of which push, deposit (parallel):    " << times.particles << " s\n"
              << "           copies of the ions:          " << times.copies << " s\n"
              << "           fields (Faraday, Ohm, ...):  " << times.fields << " s\n"
              << "           rest (output):               " << rest << " s\n"
              << "  parallel fraction p = " << parallel_fraction
              << ", so at most " << 1.0 / (1.0 - parallel_fraction) << " times faster (Amdahl)\n";
}

int main(int argc, char* argv[]) {
    const std::string case_name = (argc > 1) ? argv[1] : "uniform";
    Setup setup;
    if (case_name == "uniform") {
        setup = uniform_plasma_at_rest();
    } else if (case_name == "sound") {
        setup = ion_sound_wave();
    } else if (case_name == "decay") {
        setup = resistive_decay();
    } else if (case_name == "timing") {
        setup = timing_run();
    } else {
        std::cerr << "unknown case '" << case_name << "': choose uniform, sound, decay or timing\n";
        return 1;
    }
    Population& ions = setup.ions;
    Fields& fields = setup.fields;

    const Grid& grid = setup.grid;

    fields.E = ohm(deposit(ions, grid), fields.B, grid, setup.parameters);

    std::ofstream file(case_name + ".csv");
    write_output(file, 0.0, fields.B, deposit(ions, grid).density);

    Stopwatch stopwatch;
    StepTimes times;
    for (int step = 1; step <= setup.n_steps; ++step) {
        advance_one_step(ions, fields, grid, setup.parameters, setup.dt, times);
        if (step % setup.steps_between_outputs == 0) {
            write_output(file, step * setup.dt, fields.B, deposit(ions, grid).density);
        }
    }
    const double total_seconds = stopwatch.lap_seconds();

    std::cout << case_name << ": " << setup.n_steps << " steps written to " << case_name << ".csv\n";
    print_timing(total_seconds, times, setup.n_steps, ions.x.size());
}
