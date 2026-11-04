"""Физика двойного маятника. Вывод формул - в DERIVATION.md.

Состояние: state = [theta1, omega1, theta2, omega2],
углы отсчитываются от вертикали вниз.
Параметры params - словарь из конфига (m1, m2, L1, L2, g, c1, c2).
"""
import numpy as np
from scipy.linalg import eigh


def mass_matrix(th1, th2, p):
    c = p["m2"] * p["L1"] * p["L2"] * np.cos(th1 - th2)
    return np.array([[(p["m1"] + p["m2"]) * p["L1"] ** 2, c],
                     [c, p["m2"] * p["L2"] ** 2]])


def derivatives(state, t, params, tau=0.0):
    """Правая часть системы первого порядка: возвращает d(state)/dt.

    tau - внешний момент в точке опоры (пока 0, понадобится для LQR).
    """
    p = params
    th1, w1, th2, w2 = state
    s = np.sin(th1 - th2)
    k = p["m2"] * p["L1"] * p["L2"]

    # обобщённые силы трения (из функции Рэлея) и момент в опоре
    Q1 = -p["c1"] * w1 + p["c2"] * (w2 - w1) + tau
    Q2 = -p["c2"] * (w2 - w1)

    # правая часть: M(q) * q'' = f
    f1 = -k * w2 ** 2 * s - (p["m1"] + p["m2"]) * p["g"] * p["L1"] * np.sin(th1) + Q1
    f2 = k * w1 ** 2 * s - p["m2"] * p["g"] * p["L2"] * np.sin(th2) + Q2

    a1, a2 = np.linalg.solve(mass_matrix(th1, th2, p), [f1, f2])
    return np.array([w1, a1, w2, a2])


def make_rhs(params):
    """Превращает derivatives(state, t, params) в функцию f(y, t) для интегратора."""
    return lambda y, t: derivatives(y, t, params)


def initial_state(cfg):
    """Начальное состояние из конфига (градусы -> радианы)."""
    return np.radians([cfg["theta1_deg"], cfg["omega1_deg"],
                       cfg["theta2_deg"], cfg["omega2_deg"]])


# ---------- энергия ----------
def kinetic_energy(y, p):
    y = np.atleast_2d(y)
    th1, w1, th2, w2 = y.T
    return (0.5 * (p["m1"] + p["m2"]) * p["L1"] ** 2 * w1 ** 2
            + 0.5 * p["m2"] * p["L2"] ** 2 * w2 ** 2
            + p["m2"] * p["L1"] * p["L2"] * w1 * w2 * np.cos(th1 - th2))


def potential_energy(y, p):
    y = np.atleast_2d(y)
    th1, w1, th2, w2 = y.T
    return (-(p["m1"] + p["m2"]) * p["g"] * p["L1"] * np.cos(th1)
            - p["m2"] * p["g"] * p["L2"] * np.cos(th2))


def total_energy(y, p):
    return kinetic_energy(y, p) + potential_energy(y, p)


def energy_scale(p):
    """Характерная энергия системы: работа подъёма обоих грузов из нижнего
    положения на уровень опоры, g*(m1*L1 + m2*(L1+L2)). Не зависит от выбора
    нуля потенциала и от начального состояния."""
    return p["g"] * (p["m1"] * p["L1"] + p["m2"] * (p["L1"] + p["L2"]))


def friction_power(y, p):
    """Мощность сил трения P = -2R = -(c1*w1^2 + c2*(w2-w1)^2) <= 0."""
    y = np.atleast_2d(y)
    th1, w1, th2, w2 = y.T
    return -(p["c1"] * w1 ** 2 + p["c2"] * (w2 - w1) ** 2)


# ---------- координаты ----------
def positions(y, p):
    y = np.atleast_2d(y)
    th1, w1, th2, w2 = y.T
    x1 = p["L1"] * np.sin(th1)
    y1 = -p["L1"] * np.cos(th1)
    x2 = x1 + p["L2"] * np.sin(th2)
    y2 = y1 - p["L2"] * np.cos(th2)
    return x1, y1, x2, y2


# ---------- линеаризация около нижнего положения ----------
def linear_matrices(p):
    """M0 и K: при малых углах M0 * q'' + K * q = 0."""
    M0 = mass_matrix(0.0, 0.0, p)
    K = np.diag([(p["m1"] + p["m2"]) * p["g"] * p["L1"], p["m2"] * p["g"] * p["L2"]])
    return M0, K


def normal_modes(p):
    """Собственные частоты omega_i (рад/с) и формы колебаний (столбцы).
    Решается обобщённая задача K v = omega^2 M0 v."""
    M0, K = linear_matrices(p)
    lam, vecs = eigh(K, M0)
    return np.sqrt(lam), vecs


def linearized_system(p):
    """Аналитические A и B для x' = A x + B tau около нижнего положения
    x = (theta1, omega1, theta2, omega2) (см. DERIVATION.md, раздел 9)."""
    M0, K = linear_matrices(p)
    C = np.array([[p["c1"] + p["c2"], -p["c2"]],
                  [-p["c2"], p["c2"]]])
    Minv = np.linalg.inv(M0)
    Aq = -Minv @ K            # ускорения от углов
    Av = -Minv @ C            # ускорения от скоростей
    A = np.array([[0, 1, 0, 0],
                  [Aq[0, 0], Av[0, 0], Aq[0, 1], Av[0, 1]],
                  [0, 0, 0, 1],
                  [Aq[1, 0], Av[1, 0], Aq[1, 1], Av[1, 1]]])
    b = Minv @ np.array([1.0, 0.0])          # момент действует на theta1
    B = np.array([0.0, b[0], 0.0, b[1]])
    return A, B


def numeric_jacobian(p, x0=None, h=1e-6):
    """Якобиан правой части по состоянию и по tau, центральными разностями.
    Работает в любой точке x0 (нужно будет для верхнего положения)."""
    x0 = np.zeros(4) if x0 is None else np.asarray(x0, dtype=float)
    A = np.zeros((4, 4))
    for j in range(4):
        d = np.zeros(4)
        d[j] = h
        A[:, j] = (derivatives(x0 + d, 0.0, p) - derivatives(x0 - d, 0.0, p)) / (2 * h)
    B = (derivatives(x0, 0.0, p, tau=h) - derivatives(x0, 0.0, p, tau=-h)) / (2 * h)
    return A, B
