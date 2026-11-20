"""Run: python -m unittest test_pid.py"""
import unittest
import numpy as np
from controller import PIDController, suggested_gains, upper_linearization, wrap_error
from integrator import rk4_step
from physics import derivatives

P = dict(m1=1., m2=1., L1=1., L2=1., g=9.81, c1=0.1, c2=0.1)

class TestUpperPID(unittest.TestCase):
    def test_upper_equilibrium(self):
        self.assertAlmostEqual(wrap_error(np.pi), 0.)
        self.assertLess(np.linalg.norm(derivatives(np.array([np.pi,0,np.pi,0]),0,P)),1e-10)

    def test_linearization(self):
        A,B = upper_linearization(P)
        self.assertEqual(A.shape, (4,4))
        self.assertEqual(B.shape, (4,))

    def test_small_perturbations(self):
        K = suggested_gains(P)
        for degrees in (1, 5, 10):
            y = np.radians([180+degrees,0,180-degrees,0])
            pid = PIDController(K, ki=0., limit=100.)
            for i in range(2000):
                torque = pid.output(y)
                y = rk4_step(lambda s,t: derivatives(s,t,P,tau=torque), y, i*0.005, 0.005)
            self.assertLess(abs(wrap_error(y[0])), np.radians(0.01))
            self.assertLess(abs(wrap_error(y[2])), np.radians(0.01))

if __name__ == '__main__':
    unittest.main()
