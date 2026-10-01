"""Step 2b: periodic grid, order-1 deposit and gather (see TESTS.md).

Both directions use the same shape function: a particle at distance d from node i gives it the
fraction max(0, 1 - d/dx) of its weight (the "hat" function), and feels the field of node i with
that same fraction.
"""
import unittest

import numpy as np
from ddt import data, ddt, unpack

from common import pic, thermal_population

GRIDS = [
    # nx, dx
    (8, 0.5),
    (13, 0.37),
]


def periodic_distance(x, node_position, length):
    distance = np.abs(x - node_position)
    return np.minimum(distance, length - distance)


def hat(distance, dx):
    return np.maximum(0.0, 1.0 - distance / dx)


def single_particle(x, w, v=(0.3, -0.2, 0.1)):
    population = pic.Population("protons", 1.0, 1.0)
    population.x, population.w = [x], [w]
    population.vx, population.vy, population.vz = [v[0]], [v[1]], [v[2]]
    return population


@ddt
class TestGrid(unittest.TestCase):

    @data(*GRIDS)
    @unpack
    def test_deposit_is_the_hat_function(self, nx, dx):
        """One particle of weight w swept through the whole domain, across the periodic boundary:
        the density on every node is w * hat(distance to the node)."""
        grid = pic.Grid(nx, dx)
        length, w = grid.length(), 0.7
        node_positions = np.array([grid.node_position(i) for i in range(nx)])
        worst = 0.0
        for x in np.linspace(0.0, length, 20 * nx, endpoint=False):
            population = single_particle(x, w)
            population.deposit(grid)
            expected = w * hat(periodic_distance(x, node_positions, length), dx)
            worst = max(worst, np.max(np.abs(population.density - expected)))
            np.testing.assert_allclose(population.flux.x, population.density * population.vx[0], rtol=1e-15)
        print(f"\n  nx={nx} dx={dx}: max |density - w hat| = {worst:.1e}")
        self.assertLess(worst, 1e-14)

    @data(*[(nx, dx, centering) for nx, dx in GRIDS for centering in ("node", "cell_centre")])
    @unpack
    def test_gather_is_linear_interpolation(self, nx, dx, centering):
        """Arbitrary values on the nodes, or at the cell centres (Yee layout, Step 3a): the gathered
        value is the straight line between the two points around the particle, across the periodic
        boundary too."""
        grid = pic.Grid(nx, dx)
        on_nodes = centering == "node"
        component_centering = pic.Centering.node if on_nodes else pic.Centering.cell_centre
        field = pic.VectorField(grid, component_centering, pic.Centering.node, pic.Centering.node)
        field.x = np.random.default_rng(1).uniform(-1.0, 1.0, nx)
        values = field.x
        first_point = 0.0 if on_nodes else 0.5 * dx

        worst = 0.0
        for x in np.linspace(0.0, grid.length(), 20 * nx, endpoint=False):
            position_in_points = (x - first_point) / dx
            left = int(np.floor(position_in_points))
            fraction = position_in_points - left
            expected = values[left % nx] + fraction * (values[(left + 1) % nx] - values[left % nx])
            worst = max(worst, abs(pic.gather(field, grid, x).x - expected))
        print(f"\n  nx={nx} dx={dx}, values on {centering}s: max |gathered - straight line| = {worst:.1e}")
        self.assertLess(worst, 1e-14)

    def test_deposit_conserves_weight_momentum_and_energy(self):
        """Many particles anywhere, two species: the deposit only shares out what the particles
        carry. Summed over the nodes, density = sum w, flux = sum w v, energy = sum w m v^2 / 2,
        and total density x bulk velocity = total flux of both species."""
        grid = pic.Grid(16, 0.25)
        protons = thermal_population(5000, grid.length(), temperature=1.0, m=1.0, seed=2)
        alphas = thermal_population(3000, grid.length(), temperature=1.0, q=2.0, m=4.0, seed=3)
        alphas.vx = alphas.vx + 0.5
        for population in (protons, alphas):
            population.deposit(grid)
            speed_squared = population.vx**2 + population.vy**2 + population.vz**2
            self.assertAlmostEqual(population.density.sum(), population.w.sum(), delta=1e-13)
            self.assertAlmostEqual(population.flux.x.sum(), np.sum(population.w * population.vx), delta=1e-13)
            self.assertAlmostEqual(population.kinetic_energy_density.sum(),
                                   np.sum(0.5 * population.m * population.w * speed_squared), delta=1e-13)

        total_density = pic.total_density([protons, alphas], grid)
        bulk_velocity = pic.bulk_velocity([protons, alphas], grid)
        np.testing.assert_allclose(total_density, protons.density + alphas.density, rtol=1e-15)
        total_flux = np.sum(protons.w * protons.vx) + np.sum(alphas.w * alphas.vx)
        print(f"\n  sum of n u_x = {np.sum(total_density * bulk_velocity.x):.15f}, "
              f"sum of w v_x = {total_flux:.15f}")
        self.assertAlmostEqual(np.sum(total_density * bulk_velocity.x), total_flux, delta=1e-13)


if __name__ == "__main__":
    unittest.main()
