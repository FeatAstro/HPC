"""Step 1: one particle, Boris pusher, uniform constant fields (see TESTS.md).

B is along +z: the particle turns in the (x, y) plane. In 1D3V we follow x, which oscillates
between -r_L and +r_L, while the velocity (vx, vy) turns on a circle.
"""
import unittest

import numpy as np
from ddt import data, ddt, unpack

from common import NO_FIELD, cyclotron_period, field_along_z, larmor_radius, make_particle, pic


@ddt
class TestBorisPusher(unittest.TestCase):

    @data(
        # name,        q,    m
        ("ion m=1", 1.0, 1.0),
        ("ion m=4", 1.0, 4.0),
        ("ion m=16", 1.0, 16.0),
        ("ion m=1836", 1.0, 1836.0),
        ("electron", -1.0, 1.0 / 1836.0),
    )
    @unpack
    def test_larmor_radius_and_rotation(self, name, q, m):
        """x oscillates with amplitude r_L = m v / |q| B; the velocity turns clockwise for
        q > 0 and anticlockwise for q < 0."""
        B, v = 1.0, 1.0
        steps = 200  # multiple of 4: the extremes of x (quarter periods) fall on a step
        particle = make_particle(v=(v, 0.0, 0.0), q=q, m=m)
        trajectory = pic.run_particle(particle, NO_FIELD, field_along_z(B), cyclotron_period(q, m, B) / steps, steps)

        measured_radius = (trajectory.x.max() - trajectory.x.min()) / 2.0
        expected_radius = larmor_radius(v, q, m, B)
        print(f"\n  {name}: r_L = {measured_radius:.10g}, expected {expected_radius:.10g}")
        # the extremes are sampled with Boris' tiny phase lag: ~1e-8 error
        self.assertAlmostEqual(measured_radius / expected_radius, 1.0, delta=1e-6)

        # z component of v(n) x v(n+1): negative means clockwise
        turn = trajectory.vx[:-1] * trajectory.vy[1:] - trajectory.vy[:-1] * trajectory.vx[1:]
        self.assertTrue(np.all(np.sign(turn) == -np.sign(q)))

    def test_error_is_second_order_in_dt(self):
        """After exactly one period the particle should be back at x = 0: |x| is the error,
        and it must fall as dt^2."""
        period = cyclotron_period(1.0, 1.0, 1.0)
        steps_per_period = np.array([16, 32, 64, 128, 256, 512])
        errors = []
        for steps in steps_per_period:
            particle = make_particle(v=(1.0, 0.0, 0.0))
            pic.run_particle(particle, NO_FIELD, field_along_z(1.0), period / steps, int(steps))
            errors.append(abs(particle.x))

        order = np.polyfit(np.log(period / steps_per_period), np.log(errors), 1)[0]
        print(f"\n  errors {np.array(errors)}\n  order {order:.3f} (expected 2)")
        self.assertAlmostEqual(order, 2.0, delta=0.05)

    @data(
        # Ex,   Ey,   B,    q,     m
        (0.0, 0.1, 1.0, 1.0, 1.0),
        (0.3, -0.2, 2.0, 1.0, 1.0),
        (0.0, 0.1, 1.0, -1.0, 1.0),
        (0.0, 0.1, 1.0, 1.0, 16.0),
    )
    @unpack
    def test_ExB_drift(self, Ex, Ey, B, q, m):
        """A particle started at v_E = E x B / B^2 feels no force (E + v x B = 0): it keeps
        the velocity v_E, whatever q and m."""
        drift_vx, drift_vy = Ey / B, -Ex / B
        particle = make_particle(v=(drift_vx, drift_vy, 0.0), q=q, m=m)
        trajectory = pic.run_particle(particle, pic.Vec3(Ex, Ey, 0.0), field_along_z(B), 0.1, 500)
        print(f"\n  E=({Ex}, {Ey}) B={B} q={q:+} m={m:g}: v_E = ({drift_vx:+.2f}, {drift_vy:+.2f}), "
              f"final v = ({trajectory.vx[-1]:+.15f}, {trajectory.vy[-1]:+.15f})")
        np.testing.assert_allclose(trajectory.vx, drift_vx, atol=1e-14)
        np.testing.assert_allclose(trajectory.vy, drift_vy, atol=1e-14)
        np.testing.assert_allclose(trajectory.x, drift_vx * trajectory.t, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
