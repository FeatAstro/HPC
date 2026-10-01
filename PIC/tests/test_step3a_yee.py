"""Step 3a: Yee layout. Derivatives and averages between the nodes and the cell centres
(see TESTS.md)."""
import math
import unittest

import numpy as np
from ddt import data, ddt

from common import cell_centre_positions, node_positions, pic

LENGTH = 2.0 * math.pi
WAVENUMBER = 2.0


def largest_error(operation, nx):
    """Applies the operation to sin(kx) and compares with the exact result where it lands."""
    grid = pic.Grid(nx, LENGTH / nx)
    wave = lambda x: np.sin(WAVENUMBER * x)
    slope = lambda x: WAVENUMBER * np.cos(WAVENUMBER * x)
    nodes, cell_centres = node_positions(grid), cell_centre_positions(grid)
    measured, expected = {
        "derivative at cell centres": (pic.derivative_at_cell_centres(wave(nodes), grid), slope(cell_centres)),
        "derivative at nodes": (pic.derivative_at_nodes(wave(cell_centres), grid), slope(nodes)),
        "average at cell centres": (pic.average_at_cell_centres(wave(nodes)), wave(cell_centres)),
        "average at nodes": (pic.average_at_nodes(wave(cell_centres)), wave(nodes)),
    }[operation]
    return np.max(np.abs(measured - expected))


@ddt
class TestYeeLayout(unittest.TestCase):

    @data("derivative at cell centres", "derivative at nodes", "average at cell centres", "average at nodes")
    def test_operations_are_second_order_in_dx(self, operation):
        """A sine wave on one set of points, differentiated or averaged onto the other set: the error
        against the exact value there must fall as dx^2."""
        cells = np.array([16, 32, 64, 128, 256, 512])
        errors = np.array([largest_error(operation, nx) for nx in cells])
        order = np.polyfit(np.log(LENGTH / cells), np.log(errors), 1)[0]
        print(f"\n  {operation}: errors {np.array2string(errors, precision=2)}\n  order {order:.3f} (expected 2)")
        self.assertAlmostEqual(order, 2.0, delta=0.05)


if __name__ == "__main__":
    unittest.main()
