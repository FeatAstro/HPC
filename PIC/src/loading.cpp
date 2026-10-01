#include "loading.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>

void load_maxwellian(Population& population, const Grid& grid, int particles_per_cell,
                     const Profile& density, const VectorProfile& bulk_velocity, const Profile& temperature,
                     unsigned seed) {
    if (particles_per_cell < 1) {
        throw std::invalid_argument("load_maxwellian needs at least one particle per cell");
    }
    std::mt19937_64 random_generator(seed);
    std::uniform_real_distribution<double> uniform_in_cell(0.0, 1.0);
    std::normal_distribution<double> standard_normal(0.0, 1.0);

    const std::size_t n_particles = static_cast<std::size_t>(grid.nx) * particles_per_cell;
    for (std::vector<double>* array : {&population.x, &population.vx, &population.vy, &population.vz, &population.w}) {
        array->clear();
        array->reserve(n_particles);
    }

    for (int cell = 0; cell < grid.nx; ++cell) {
        const double cell_centre = grid.cell_centre_position(cell);
        const double weight = density(cell_centre) / particles_per_cell;
        const Vec3 u = bulk_velocity(cell_centre);
        const double thermal_speed = std::sqrt(temperature(cell_centre) / population.m);

        for (int k = 0; k < particles_per_cell; ++k) {
            population.x.push_back((cell + uniform_in_cell(random_generator)) * grid.dx);
            population.vx.push_back(u.x + thermal_speed * standard_normal(random_generator));
            population.vy.push_back(u.y + thermal_speed * standard_normal(random_generator));
            population.vz.push_back(u.z + thermal_speed * standard_normal(random_generator));
            population.w.push_back(weight);
        }
    }
}
