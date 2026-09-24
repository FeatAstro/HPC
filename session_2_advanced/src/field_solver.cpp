#include "field_solver.hpp"

#include <algorithm>

#include "finite_differences.hpp"
#include "interpolation.hpp"

namespace {

// Closure for the electron fluid: change this one function to change the equation of state.
double electron_pressure(double density, const PlasmaParameters& plasma) {
    return density * plasma.electron_temperature;
}

ScalarField electron_pressure_at_nodes(const Grid2D& g, const PlasmaParameters& plasma) {
    ScalarField pressure(g.ion_density.size());
    for (std::size_t node = 0; node < pressure.size(); ++node) {
        pressure[node] = electron_pressure(g.ion_density[node], plasma);
    }
    return pressure;
}

// Component `component` of grad P_e, at the position where that component of E lives.
double pressure_gradient(const Grid2D& g, const ScalarField& pressure, int component, int i, int j) {
    const auto pressure_at_node = [&pressure](int node) { return pressure[node]; };
    switch (component) {
    case axis::x:
        return forward_difference_x(g, pressure_at_node, i, j);
    case axis::y:
        return forward_difference_y(g, pressure_at_node, i, j);
    default:
        return 0.0; // nothing varies along z
    }
}

}

VectorField curl_of_electric_field(const Grid2D& g, const VectorField& electric_field) {
    const auto Ex = component_of(electric_field, axis::x);
    const auto Ey = component_of(electric_field, axis::y);
    const auto Ez = component_of(electric_field, axis::z);

    VectorField curl(electric_field.size());
    g.for_each_node([&](int i, int j, int node) {
        curl[node] = Vec3{forward_difference_y(g, Ez, i, j),
                          -forward_difference_x(g, Ez, i, j),
                          forward_difference_x(g, Ey, i, j) - forward_difference_y(g, Ex, i, j)};
    });
    return curl;
}

VectorField curl_of_magnetic_field(const Grid2D& g, const VectorField& magnetic_field) {
    const auto Bx = component_of(magnetic_field, axis::x);
    const auto By = component_of(magnetic_field, axis::y);
    const auto Bz = component_of(magnetic_field, axis::z);

    VectorField curl(magnetic_field.size());
    g.for_each_node([&](int i, int j, int node) {
        curl[node] = Vec3{backward_difference_y(g, Bz, i, j),
                          -backward_difference_x(g, Bz, i, j),
                          backward_difference_x(g, By, i, j) - backward_difference_y(g, Bx, i, j)};
    });
    return curl;
}

VectorField current_density_from_ampere(const Grid2D& g, const VectorField& magnetic_field,
                                        const PlasmaParameters& plasma) {
    return scaled(curl_of_magnetic_field(g, magnetic_field), 1.0 / plasma.vacuum_permeability);
}

VectorField electric_field_from_ohms_law(const Grid2D& g, const VectorField& magnetic_field,
                                         const PlasmaParameters& plasma) {
    const VectorField current_density = current_density_from_ampere(g, magnetic_field, plasma);
    const ScalarField pressure = electron_pressure_at_nodes(g, plasma);
    const double e = plasma.elementary_charge;

    VectorField electric_field(magnetic_field.size());
    for (int component = 0; component < axis::count; ++component) {
        const CellOffset offset = yee::electric_field[component];

        g.for_each_node([&](int i, int j, int node) {
            const double x = g.node_x(i, offset.x);
            const double y = g.node_y(j, offset.y);

            const double density = std::max(gather_scalar(g.ion_density, g, x, y), plasma.minimum_density);
            const Vec3 ion_velocity = gather(g.ion_velocity, yee::nodes, g, x, y);
            const Vec3 current = gather(current_density, yee::current_density, g, x, y);
            const Vec3 B = gather(magnetic_field, yee::magnetic_field, g, x, y);

            const Vec3 electron_velocity = ion_velocity - (1.0 / (density * e)) * current;
            const double convective_term = -cross(electron_velocity, B)[component];
            const double pressure_term = -pressure_gradient(g, pressure, component, i, j) / (density * e);

            electric_field[node][component] = convective_term + pressure_term;
        });
    }
    return electric_field;
}

void update_electric_field(Grid2D& g, const PlasmaParameters& plasma) {
    g.current_density = current_density_from_ampere(g, g.magnetic_field, plasma);
    g.electric_field = electric_field_from_ohms_law(g, g.magnetic_field, plasma);
}

void advance_magnetic_field(Grid2D& g, const PlasmaParameters& plasma, double dt) {
    const auto faraday_rate_of_change = [&](const VectorField& magnetic_field) {
        return scaled(curl_of_electric_field(g, electric_field_from_ohms_law(g, magnetic_field, plasma)), -1.0);
    };

    const VectorField start = g.magnetic_field;
    const VectorField k1 = faraday_rate_of_change(start);
    const VectorField k2 = faraday_rate_of_change(added(start, dt / 2.0, k1));
    const VectorField k3 = faraday_rate_of_change(added(start, dt / 2.0, k2));
    const VectorField k4 = faraday_rate_of_change(added(start, dt, k3));

    for (std::size_t node = 0; node < start.size(); ++node) {
        g.magnetic_field[node] = start[node] + (dt / 6.0) * (k1[node] + 2.0 * k2[node] + 2.0 * k3[node] + k4[node]);
    }
}
