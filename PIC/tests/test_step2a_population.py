"""Step 2a: a population of N particles in uniform constant fields, periodic domain [0, L)
(see TESTS.md). The uniform fields are given on the grid of Step 2b, identical at every node."""
import unittest

import numpy as np

from common import NO_FIELD, field_along_z, pic, temperature, thermal_population, uniform_fields


class TestPopulation(unittest.TestCase):

    def test_thermal_population_keeps_energy_and_temperature(self):
        """A magnetic field does no work: the kinetic energy and the temperature of a Maxwellian
        population stay constant to round-off, over ~16 cyclotron periods."""
        grid = pic.Grid(10, 1.0)
        population = thermal_population(100_000, grid.length(), temperature=0.5)
        initial_temperature = temperature(population)
        fields = uniform_fields(grid, NO_FIELD, field_along_z(1.0))
        history = pic.run_population(population, grid, fields, 0.1, 1000)
        final_temperature = temperature(population)

        energy = history.kinetic_energy
        energy_change = np.max(np.abs(energy - energy[0])) / energy[0]
        print(f"\n  T before {initial_temperature:.12f}, after {final_temperature:.12f} (loaded T = 0.5)"
              f"\n  max relative change of the kinetic energy: {energy_change:.1e}")
        self.assertLess(energy_change, 1e-12)
        self.assertAlmostEqual(final_temperature / initial_temperature, 1.0, delta=1e-12)

    def test_periodic_wrap(self):
        """In a domain of length 1, fast particles wrap many times: every x stays in [0, L) and
        equals, modulo L, the position of the same particle in a domain where it never wraps."""
        dt, n_steps = 0.05, 400
        small_grid, large_grid = pic.Grid(10, 0.1), pic.Grid(20, 100.0)
        E = pic.Vec3(0.0, 0.3, 0.0)  # E x B drift along x: particles keep crossing the boundary
        wrapped = thermal_population(2000, small_grid.length(), temperature=4.0)
        unwrapped = thermal_population(2000, small_grid.length(), temperature=4.0)  # same seed, same particles
        offset = 1000.0
        unwrapped.x = unwrapped.x + offset  # middle of [0, 2000): no boundary reached in 400 steps

        pic.run_population(wrapped, small_grid, uniform_fields(small_grid, E, field_along_z(1.0)), dt, n_steps)
        pic.run_population(unwrapped, large_grid, uniform_fields(large_grid, E, field_along_z(1.0)), dt, n_steps)

        length = small_grid.length()
        self.assertTrue(np.all((wrapped.x >= 0.0) & (wrapped.x < length)))
        print(f"\n  unwrapped positions span [{unwrapped.x.min() - offset:.1f}, "
              f"{unwrapped.x.max() - offset:.1f}], L = {length}")
        # distance around the periodic domain, so that 0.999999 and 0.000001 are close
        difference = (wrapped.x - (unwrapped.x - offset)) % length
        distance = np.minimum(difference, length - difference)
        self.assertLess(distance.max(), 1e-10)
        np.testing.assert_allclose(wrapped.vx, unwrapped.vx, rtol=0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
