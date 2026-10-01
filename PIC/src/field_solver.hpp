#pragma once

#include "grid.hpp"

// Normalised units: mu0 = e = m_i = 1.
struct PlasmaParameters {
    double electron_temperature = 0.0;
    double resistivity = 0.0;        // eta
    double hyper_resistivity = 0.0;  // nu
};

// Isothermal electrons: p_e = n T_e, on the nodes like the density.
ScalarField electron_pressure(const ScalarField& density, const PlasmaParameters& parameters);

// Ampere's law without displacement current: j = curl B.
// In 1D: jx = 0, jy = -dBz/dx, jz = dBy/dx.
VectorField ampere(const VectorField& B, const Grid& grid);

// Generalised Ohm's law (session 3, slide 7):
// E = -u x B + (j x B) / n - grad(p_e) / n + eta j - nu laplacian(j).
// density and bulk_velocity are the ion moments on the nodes.
VectorField ohm(const ScalarField& density, const VectorField& bulk_velocity, const VectorField& B,
                const VectorField& j, const Grid& grid, const PlasmaParameters& parameters);

// Faraday's law over one time step: returns B - dt curl E.
// In 1D: Bx is constant, dBy/dt = dEz/dx, dBz/dt = -dEy/dx.
VectorField faraday(const VectorField& B, const VectorField& E, const Grid& grid, double dt);

// One time step of B and E with the ion moments held fixed, by the iterated Crank-Nicolson scheme
// (session 3, slides 6-7): two predictions of B from the time-centred E, then the correction.
// Stable for dt below about n dx^2 / (2 B), the period of the grid-scale whistler.
void advance_fields_with_fixed_ions(Fields& fields, const ScalarField& density, const VectorField& bulk_velocity,
                                    const Grid& grid, const PlasmaParameters& parameters, double dt);
