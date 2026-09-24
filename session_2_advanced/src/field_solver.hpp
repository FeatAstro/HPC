#pragma once

#include "grid.hpp"
#include "plasma_parameters.hpp"

// Curl of E, landing where B lives (Faraday: dB/dt = -curl E).
VectorField curl_of_electric_field(const Grid2D& g, const VectorField& electric_field);

// Curl of B, landing where E and j live.
VectorField curl_of_magnetic_field(const Grid2D& g, const VectorField& magnetic_field);

// Ampere's law without displacement current: mu0 j = curl B.
VectorField current_density_from_ampere(const Grid2D& g, const VectorField& magnetic_field,
                                        const PlasmaParameters& plasma);

// Generalized Ohm's law for massless, quasi-neutral electrons (slide 1), with the ion
// moments currently on the grid:
//   v_e = v_i - j / (n e),   E = -v_e x B - (div P_e) / (n e)
VectorField electric_field_from_ohms_law(const Grid2D& g, const VectorField& magnetic_field,
                                         const PlasmaParameters& plasma);

// Sets grid.current_density and grid.electric_field from grid.magnetic_field.
void update_electric_field(Grid2D& g, const PlasmaParameters& plasma);

// Advances grid.magnetic_field by dt with Faraday's law (fourth-order Runge-Kutta),
// keeping the ion moments on the grid frozen.
void advance_magnetic_field(Grid2D& g, const PlasmaParameters& plasma, double dt);
