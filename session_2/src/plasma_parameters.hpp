#pragma once

// Physical constants and plasma properties of the hybrid model. The defaults are
// normalised units (mu0 = e = m_i = n0 = 1) where B is in units of B0, velocity in
// units of the Alfven speed and time in units of the inverse ion cyclotron frequency.
struct PlasmaParameters {
    double vacuum_permeability = 1.0;
    double elementary_charge = 1.0; // ions carry +elementary_charge, electrons -elementary_charge
    double ion_mass = 1.0;
    double reference_density = 1.0; // density the ions are loaded with
    double electron_temperature = 0.0; // isothermal electron fluid
    double minimum_density = 1e-3;     // floor for divisions by the density in Ohm's law
};
