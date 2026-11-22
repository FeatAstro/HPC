#pragma once

#include "grid.hpp"
#include "population.hpp"

// Normalised units: mu0 = e = m_i = 1.
struct PlasmaParameters {
    double electron_temperature = 0.0;
    double resistivity = 0.0;        // eta: collisions, lets field lines reconnect
    double hyper_resistivity = 0.0;  // nu: damps the noise at the scale of a cell
};

// Ampere's law without displacement current: j = curl B.
// In 1D: jx = 0, jy = -dBz/dx, jz = dBy/dx. j lives where E does.
VectorField ampere(const VectorField& B, const Grid& grid);

// Generalised Ohm's law (session 3, slide 7), with j from Ampere's law and isothermal electrons:
// E = -u x B + (j x B) / n - grad(p_e) / n + eta j - nu laplacian(j),  p_e = n T_e.
// n and u are the deposited ion moments.
VectorField ohm(const Moments& moments, const VectorField& B, const Grid& grid, const PlasmaParameters& parameters);

// Faraday's law over one time step: returns B - dt curl E.
// In 1D: Bx is constant, dBy/dt = dEz/dx, dBz/dt = -dEy/dx.
VectorField faraday(const VectorField& B, const VectorField& E, const Grid& grid, double dt);

VectorField average(const VectorField& a, const VectorField& b);
