"""Step 3c: Faraday's law advanced with the iterated Crank-Nicolson scheme, ions held fixed
(see TESTS.md)."""
import math
import unittest

import numpy as np
from ddt import data, ddt

from common import follow_amplitude, pic, plasma_parameters, transverse_wave

GRID = pic.Grid(128, 2.0 * math.pi / 128)
WAVE_AMPLITUDE = 1e-3


@ddt
class TestFaraday(unittest.TestCase):

    @data(1, 2, 4)
    def test_whistler_frequency(self, mode):
        """Ions at rest, Bx = 1, a small transverse wave and the Hall term only: the wave pattern
        turns at omega = k^2 Bx / n. Measured over two periods."""
        Bx, dt = 1.0, 5e-4
        fields, k = transverse_wave(GRID, mode, WAVE_AMPLITUDE, Bx)
        expected_frequency = k**2 * Bx
        steps = round(2 * 2.0 * math.pi / expected_frequency / dt)
        times, amplitudes = follow_amplitude(fields, GRID, plasma_parameters(), mode, dt, steps // 40, 40)

        measured_frequency = abs(np.polyfit(times, np.unwrap(np.angle(amplitudes)), 1)[0])
        print(f"\n  k = {k:g}: omega measured {measured_frequency:.5f}, expected {expected_frequency:.5f}; "
              f"amplitude kept to {abs(amplitudes[-1]) / abs(amplitudes[0]):.6f}")
        self.assertAlmostEqual(measured_frequency / expected_frequency, 1.0, delta=0.005)
        self.assertAlmostEqual(abs(amplitudes[-1]) / abs(amplitudes[0]), 1.0, delta=1e-3)

    def test_resistive_decay(self):
        """No Bx (no Hall effect), resistivity and hyper-resistivity: Faraday's law becomes a
        diffusion equation and the wave decays as exp(-(eta k^2 + nu k^4) t)."""
        mode, dt = 2, 2e-4  # the hyper-resistive term limits dt to about dx^4 / (8 nu)
        parameters = plasma_parameters(resistivity=0.05, hyper_resistivity=0.001)
        fields, k = transverse_wave(GRID, mode, WAVE_AMPLITUDE, Bx=0.0)
        expected_rate = parameters.resistivity * k**2 + parameters.hyper_resistivity * k**4
        times, amplitudes = follow_amplitude(fields, GRID, parameters, mode, dt, 700, 50)

        measured_rate = -np.polyfit(times, np.log(np.abs(amplitudes)), 1)[0]
        print(f"\n  decay rate measured {measured_rate:.5f}, expected {expected_rate:.5f}; "
              f"amplitude divided by {abs(amplitudes[0]) / abs(amplitudes[-1]):.2f}")
        self.assertAlmostEqual(measured_rate / expected_rate, 1.0, delta=0.005)


if __name__ == "__main__":
    unittest.main()
