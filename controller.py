"""Local feedback stabilization of the *upper* double-pendulum equilibrium.

Angles are absolute and measured from downward vertical. The actuator applies
ONLY a pivot torque to the first link. This controller is NOT a swing-up controller:
for large deflections it must yield to the uncontrolled nonlinear physics.
"""
import numpy as np
from scipy.signal import place_poles
from physics import numeric_jacobian


def wrap_error(angle):
    """Shortest signed angular displacement from the UPWARD vertical (pi)."""
    return (float(angle) - np.pi + np.pi) % (2.0 * np.pi) - np.pi


def upper_linearization(params):
    return numeric_jacobian(params, np.array([np.pi, 0., np.pi, 0.]))


def suggested_gains(params):
    """Full-state PD gains from pole placement, *not* LQR optimization."""
    A, B = upper_linearization(params)
    gains = place_poles(A, B[:, None], [-2.0, -2.6, -3.3, -4.1]).gain_matrix[0]
    return np.asarray(gains, dtype=float)


class PIDController:
    """Locally valid feedback using first-link PID and second-link PD.

    An automatic safety gate latches control OFF if either link is too far
    from the upright position or is moving too fast. Re-enable explicitly
    (reset/start) after a latch; outside the gate, applied torque is exactly 0.
    This avoids applying a LOCAL linear regulator to arbitrary pendulum states.
    """

    def __init__(self, gains, ki=0.0, limit=100.0, integral_limit=10.0,
                 angle_limit_deg=25.0, speed_limit_deg_s=90.0):
        self.gains = np.asarray(gains, dtype=float)
        if self.gains.shape != (4,) or not np.all(np.isfinite(self.gains)):
            raise ValueError('Нужны четыре конечных коэффициента P1, D1, P2, D2')
        self.ki = float(ki)
        self.limit = float(limit)
        self.integral_limit = float(integral_limit)
        self.angle_limit = float(np.deg2rad(angle_limit_deg))
        self.speed_limit = float(np.deg2rad(speed_limit_deg_s))
        if not np.isfinite(self.limit) or self.limit <= 0 or self.angle_limit <= 0 or self.speed_limit <= 0:
            raise ValueError('Неверные ограничения регулятора')
        self.reset()

    def reset(self):
        self.integral = 0.0
        self.disengaged = False

    def error_state(self, state):
        return np.array([wrap_error(state[0]), state[1],
                         wrap_error(state[2]), state[3]], dtype=float)

    def is_local(self, state):
        e = self.error_state(state)
        return (np.all(np.isfinite(e)) and
                max(abs(e[0]), abs(e[2])) <= self.angle_limit and
                max(abs(e[1]), abs(e[3])) <= self.speed_limit)

    def output(self, state):
        if self.disengaged:
            return 0.0
        if not self.is_local(state):
            self.disengaged = True
            self.integral = 0.0
            return 0.0
        raw = -float(self.gains @ self.error_state(state)) - self.ki*self.integral
        return float(np.clip(raw, -self.limit, self.limit))

    def advance(self, state, dt):
        if self.disengaged:
            return
        if not self.is_local(state):
            self.disengaged = True
            self.integral = 0.0
            return
        e = self.error_state(state)
        raw = -float(self.gains @ e) - self.ki*self.integral
        # Integrate only when the output is not saturated or the integral
        # increment would REDUCE saturation (basic anti-windup).
        integral_delta = self.ki * e[0] * float(dt)
        if abs(raw) < self.limit or (raw > self.limit and integral_delta > 0) or (raw < -self.limit and integral_delta < 0):
            self.integral = float(np.clip(self.integral + e[0]*dt,
                                          -self.integral_limit, self.integral_limit))
