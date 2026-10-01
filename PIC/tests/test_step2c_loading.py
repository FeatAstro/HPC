"""Step 2c: Maxwellian loading from user profiles n(x), u(x), T(x), and the 1/sqrt(N) convergence
of the deposited moments (see TESTS.md).

The deposited moments are compared with the exact average of the loader (loader_expectation), so
the only error left is the statistical noise of the particles, which must fall as 1/sqrt(N).
"""
import unittest

import numpy as np
from ddt import data, ddt

from common import example_profiles, loader_expectation, pic

GRID = pic.Grid(32, 0.25)
PARTICLES_PER_CELL = np.array([10, 40, 160, 640, 2560, 10240])
SEEDS = range(4)


def deposited_moments(particles_per_cell, seed):
    population = pic.Population("protons", 1.0, 1.0)
    pic.load_maxwellian(population, GRID, particles_per_cell, *example_profiles(GRID.length()), seed)
    population.deposit(GRID)
    return {
        "density": population.density,
        "bulk velocity": pic.bulk_velocity([population], GRID).x,
        "temperature": population.temperature(),
    }


def expected_moments():
    density, velocity, temperature = loader_expectation(GRID, 1.0, *example_profiles(GRID.length()))
    return {"density": density, "bulk velocity": velocity, "temperature": temperature}


def rms_error(moment, particles_per_cell):
    """Root mean square over the nodes and over several seeds."""
    expected = expected_moments()[moment]
    squared = [np.mean((deposited_moments(particles_per_cell, seed)[moment] - expected) ** 2) for seed in SEEDS]
    return np.sqrt(np.mean(squared))


@ddt
class TestLoading(unittest.TestCase):

    @data("density", "bulk velocity", "temperature")
    def test_error_decreases_as_one_over_sqrt_N(self, moment):
        errors = np.array([rms_error(moment, ppc) for ppc in PARTICLES_PER_CELL])
        total_particles = PARTICLES_PER_CELL * GRID.nx
        slope = np.polyfit(np.log(total_particles), np.log(errors), 1)[0]
        print(f"\n  {moment}: RMS error {np.array2string(errors, precision=4)}"
              f"\n  for N = {total_particles}: slope {slope:.3f} (expected -0.5)")
        self.assertAlmostEqual(slope, -0.5, delta=0.1)


if __name__ == "__main__":
    unittest.main()
