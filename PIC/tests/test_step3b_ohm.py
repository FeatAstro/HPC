"""Step 3b: Ampere's law and the generalised Ohm's law on the Yee grid (see TESTS.md)."""
import math
import unittest

import numpy as np
from ddt import data, ddt, unpack

from common import SmoothPlasma, magnetic_field, pic, plasma_parameters, vector_field_on_nodes

LENGTH = 2.0 * math.pi
PARAMETERS = plasma_parameters(electron_temperature=0.5, resistivity=0.1, hyper_resistivity=0.05)


@ddt
class TestOhm(unittest.TestCase):

    def test_uniform_plasma_gives_minus_u_cross_B(self):
        """Uniform n, u and B: no current and no gradient, so E = -u x B exactly. A plasma drifting
        at u sees no electric field in its own frame."""
        grid = pic.Grid(16, 0.5)
        u, B = (0.3, -0.2, 0.5), (0.7, 0.4, -0.6)
        magnetic = magnetic_field(grid, *B)

        j = pic.ampere(magnetic, grid)
        E = pic.ohm(np.full(grid.nx, 1.3), vector_field_on_nodes(grid, *u), magnetic, j, grid, PARAMETERS)

        minus_u_cross_B = -np.cross(u, B)
        print(f"\n  E = ({E.x[0]:+.3f}, {E.y[0]:+.3f}, {E.z[0]:+.3f}), "
              f"-u x B = ({minus_u_cross_B[0]:+.3f}, {minus_u_cross_B[1]:+.3f}, {minus_u_cross_B[2]:+.3f})")
        for current_component in (j.x, j.y, j.z):
            np.testing.assert_array_equal(current_component, 0.0)
        for E_component, expected in zip((E.x, E.y, E.z), minus_u_cross_B):
            np.testing.assert_allclose(E_component, expected, rtol=0, atol=1e-15)

    @data(("Ex", 0), ("Ey", 1), ("Ez", 2))
    @unpack
    def test_error_is_second_order_in_dx(self, name, component):
        """Smooth plasma with every term of Ohm's law non-zero (flow, Hall, electron pressure,
        resistivity, hyper-resistivity): the error against the exact E must fall as dx^2."""
        plasma = SmoothPlasma(LENGTH, PARAMETERS)
        cells = np.array([32, 64, 128, 256, 512, 1024])
        errors = []
        for nx in cells:
            grid = pic.Grid(nx, LENGTH / nx)
            E = plasma.solved_E(grid)
            solved = (E.x, E.y, E.z)[component]
            errors.append(np.max(np.abs(solved - plasma.exact_E_on_the_yee_grid(grid)[component])))

        order = np.polyfit(np.log(LENGTH / cells), np.log(errors), 1)[0]
        print(f"\n  {name}: errors {np.array2string(np.array(errors), precision=2)}"
              f"\n  order {order:.3f} (expected 2)")
        self.assertAlmostEqual(order, 2.0, delta=0.05)


if __name__ == "__main__":
    unittest.main()
