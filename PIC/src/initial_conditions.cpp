#include "initial_conditions.hpp"

#include <cmath>
#include <cstddef>

const double pi = 3.14159265358979323846;

// A profile that has the same value everywhere.
static Profile constant(double value) {
    return [value](double) { return value; };
}

Setup uniform_plasma_at_rest() {
    const int particles_per_cell = 100;
    const double density = 1.0;
    const double temperature = 0.0;
    const unsigned seed = 0;
    const double Bx = 1.0;

    Setup setup;
    setup.grid = {32, 2.0 * pi / 32};
    setup.dt = 0.005;
    setup.n_steps = 1000;
    setup.steps_between_outputs = 20;
    setup.parameters = PlasmaParameters();  // cold electrons, no resistivity
    setup.ions = maxwellian_protons(setup.grid, particles_per_cell, constant(density), constant(0.0),
                                    constant(0.0), constant(0.0), constant(temperature), seed);
    setup.fields = {VectorField(setup.grid), VectorField(setup.grid)};

    for (int i = 0; i < setup.grid.nx; ++i) {
        setup.fields.B.x[i] = Bx;
    }
    return setup;
}

Setup ion_sound_wave() {
    const int particles_per_cell = 50;
    const double density = 1.0;
    const double temperature = 0.0;
    const unsigned seed = 0;
    const double Bx = 1.0;
    const double electron_temperature = 1.0;

    Setup setup;
    setup.grid = {32, 2.0 * pi / 32};
    setup.dt = 0.01;

    const double k = 2.0 * pi / setup.grid.length();  // one ripple in the box
    const double delta = 0.01;                        // small: a large ripple steepens into a shock
    const double sound_speed = std::sqrt(electron_temperature);
    const double two_periods = 2.0 * 2.0 * pi / (k * sound_speed);

    setup.n_steps = int(two_periods / setup.dt);
    setup.steps_between_outputs = 10;
    setup.parameters = PlasmaParameters();
    setup.parameters.electron_temperature = electron_temperature;
    setup.ions = maxwellian_protons(setup.grid, particles_per_cell, constant(density), constant(0.0),
                                    constant(0.0), constant(0.0), constant(temperature), seed);
    setup.fields = {VectorField(setup.grid), VectorField(setup.grid)};

    for (int i = 0; i < setup.grid.nx; ++i) {
        setup.fields.B.x[i] = Bx;
    }

    // The ions are evenly spaced, so that the only ripple of the density is the one we put: random
    // positions would add a noise of 1/sqrt(6 N), larger than delta. The ripple is made by heavier
    // and lighter particles. For it to travel in one direction only, the ions move with it:
    // u = sound_speed * (n - 1).
    Population& ions = setup.ions;
    for (std::size_t i = 0; i < ions.x.size(); ++i) {
        ions.x[i] = (i + 0.5) * setup.grid.length() / ions.x.size();
        const double ripple = delta * std::cos(k * ions.x[i]);
        ions.w[i] *= 1.0 + ripple;
        ions.vx[i] += sound_speed * ripple;
    }
    return setup;
}

Setup resistive_decay() {
    const int particles_per_cell = 50;
    const double density = 1.0;
    const double temperature = 0.0;
    const unsigned seed = 0;

    Setup setup;
    setup.grid = {64, 2.0 * pi / 64};
    setup.dt = 0.004;  // hyper-resistivity needs dt < dx^4 / (8 nu) = 0.0058
    setup.n_steps = 2500;
    setup.steps_between_outputs = 25;
    setup.parameters = PlasmaParameters();
    setup.parameters.resistivity = 0.01;
    setup.parameters.hyper_resistivity = 0.002;
    setup.ions = maxwellian_protons(setup.grid, particles_per_cell, constant(density), constant(0.0),
                                    constant(0.0), constant(0.0), constant(temperature), seed);
    setup.fields = {VectorField(setup.grid), VectorField(setup.grid)};

    const double b = 0.001;  // small: the ions are barely pushed and stay at rest
    for (int i = 0; i < setup.grid.nx; ++i) {
        const double x = (i + 0.5) * setup.grid.dx;  // By lives at the cell centres
        setup.fields.B.y[i] = b * (std::sin(x) + std::sin(2.0 * x));
    }
    return setup;
}

Setup timing_run() {
    const int particles_per_cell = 2000;
    const double density = 1.0;
    const double temperature = 1.0;
    const unsigned seed = 0;
    const double Bx = 1.0;

    Setup setup;
    setup.grid = {64, 2.0 * pi / 64};
    setup.dt = 0.004;  // the Hall term needs dt < n dx^2 / (2 B) = 0.0048
    setup.n_steps = 200;
    setup.steps_between_outputs = 50;
    setup.parameters = PlasmaParameters();
    setup.parameters.hyper_resistivity = 0.001;  // damps the grid-scale noise of the warm ions
    setup.ions = maxwellian_protons(setup.grid, particles_per_cell, constant(density), constant(0.0),
                                    constant(0.0), constant(0.0), constant(temperature), seed);
    setup.fields = {VectorField(setup.grid), VectorField(setup.grid)};

    for (int i = 0; i < setup.grid.nx; ++i) {
        setup.fields.B.x[i] = Bx;
    }
    return setup;
}
