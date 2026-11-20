"""Regression tests for arbitrary poses, unforced motion and local-only control.
Run: python -m unittest test_physics_regression test_pid
"""
import unittest
import numpy as np
from physics import derivatives, total_energy, energy_scale
from integrator import rk4_step
from controller import PIDController, suggested_gains

P = dict(m1=1., m2=1., L1=1., L2=1., g=9.81, c1=0., c2=0.)


class TestPhysicsAllAngles(unittest.TestCase):
    def test_energy_conservation_across_positions(self):
        for a, b in [(0, 0), (90, 0), (120, -35), (180, 180),
                     (270, 90), (-160, 20), (25, 330)]:
            with self.subTest(angles=(a, b)):
                y = np.deg2rad([a, 0, b, 0])
                start_E = float(total_energy(y, P)[0])
                dt = .002
                for k in range(3000):
                    y = rk4_step(lambda state, t: derivatives(state, t, P),
                                 y, k*dt, dt)
                self.assertTrue(np.all(np.isfinite(y)))
                relative_drift = abs(float(total_energy(y, P)[0]) - start_E) / energy_scale(P)
                self.assertLess(relative_drift, 1e-5)

    def test_downward_equilibrium_unforced(self):
        y = np.zeros(4)
        np.testing.assert_allclose(derivatives(y, 0., P), np.zeros(4), atol=1e-12)

    def test_pivot_torque_changes_acceleration_not_geometry(self):
        y = np.deg2rad([75, 12, -30, -9])
        delta = derivatives(y, 0., P, tau=2.) - derivatives(y, 0., P, tau=0.)
        self.assertAlmostEqual(delta[0], 0)
        self.assertAlmostEqual(delta[2], 0)
        self.assertGreater(abs(delta[1]), 0)


class TestLocalController(unittest.TestCase):
    def test_pid_applies_no_torque_far_from_upright(self):
        K = suggested_gains(P)
        for a, b in [(0, 0), (90, 90), (180, 0), (120, -35), (270, 180)]:
            with self.subTest(angles=(a, b)):
                controller = PIDController(K)
                y = np.deg2rad([a, 0, b, 0])
                self.assertEqual(controller.output(y), 0.)
                self.assertTrue(controller.disengaged)

    def test_pid_near_upper_can_act(self):
        controller = PIDController(suggested_gains(P))
        y = np.deg2rad([185, 0, 178, 0])
        self.assertNotEqual(controller.output(y), 0.)
        self.assertFalse(controller.disengaged)

    def test_disengagement_is_latched(self):
        controller = PIDController(suggested_gains(P))
        self.assertEqual(controller.output(np.zeros(4)), 0.)
        self.assertEqual(controller.output(np.deg2rad([185, 0, 175, 0])), 0.)
        controller.reset()
        self.assertNotEqual(controller.output(np.deg2rad([185, 0, 175, 0])), 0.)

    def test_speed_safety_gate(self):
        controller = PIDController(suggested_gains(P))
        y = np.deg2rad([180, 100, 180, 0])
        self.assertEqual(controller.output(y), 0.)


if __name__ == '__main__':
    unittest.main()
