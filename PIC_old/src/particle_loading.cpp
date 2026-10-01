#include "particle_loading.hpp"

#include <random>

std::vector<Particle> load_uniform_ions(const Grid2D& g, const PlasmaParameters& plasma, const IonLoading& loading) {
    const int particles_per_cell = loading.particles_per_cell_x * loading.particles_per_cell_y;

    std::mt19937 random_generator(loading.random_seed);
    std::normal_distribution<double> velocity_component(0.0, loading.thermal_speed > 0.0 ? loading.thermal_speed : 1.0);

    std::vector<Particle> ions;
    ions.reserve(static_cast<std::size_t>(g.nx) * g.ny * particles_per_cell);

    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            for (int b = 0; b < loading.particles_per_cell_y; ++b) {
                for (int a = 0; a < loading.particles_per_cell_x; ++a) {
                    Particle ion;
                    ion.x = g.node_x(i, (a + 0.5) / loading.particles_per_cell_x);
                    ion.y = g.node_y(j, (b + 0.5) / loading.particles_per_cell_y);
                    ion.q = plasma.elementary_charge;
                    ion.m = plasma.ion_mass;
                    ion.w = plasma.reference_density / particles_per_cell;
                    if (loading.thermal_speed > 0.0) {
                        ion.v = Vec3{velocity_component(random_generator), velocity_component(random_generator),
                                     velocity_component(random_generator)};
                    }
                    ions.push_back(ion);
                }
            }
        }
    }
    return ions;
}
