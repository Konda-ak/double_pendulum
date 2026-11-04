"""Интегратор RK4. Ничего не знает про маятник: ему нужна только функция
правой части f(y, t), которая возвращает dy/dt."""
import numpy as np


def rk4_step(f, y, t, h):
    """Один шаг метода Рунге-Кутты 4-го порядка для y' = f(y, t)."""
    k1 = f(y, t)
    k2 = f(y + 0.5 * h * k1, t + 0.5 * h)
    k3 = f(y + 0.5 * h * k2, t + 0.5 * h)
    k4 = f(y + h * k3, t + h)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def integrate(f, y0, t_max, dt):
    """Интегрирует от 0 до t_max с постоянным шагом dt.
    Возвращает массив времени t и массив состояний y формы (N, len(y0))."""
    n = int(round(t_max / dt))
    t = dt * np.arange(n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for i in range(n):
        y[i + 1] = rk4_step(f, y[i], t[i], dt)
    return t, y
