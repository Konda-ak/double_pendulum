"""Механика двойного маятника: параметры, уравнения движения, RK4, проверка."""
import numpy as np

# ---------- Все параметры в одном месте ----------
config = {
    "m1": 1.0,          # масса 1-го груза, кг
    "m2": 1.0,          # масса 2-го груза, кг
    "L1": 1.0,          # длина 1-го стержня, м
    "L2": 1.0,          # длина 2-го стержня, м
    "g": 9.81,          # ускорение свободного падения, м/с^2
    "c1": 0.0,          # трение в точке опоры
    "c2": 0.0,          # трение в среднем шарнире
    "theta1": 120.0,    # начальный угол 1 (градусы)
    "theta2": -10.0,    # начальный угол 2 (градусы)
    "omega1": 0.0,      # начальная угловая скорость 1 (град/с)
    "omega2": 0.0,      # начальная угловая скорость 2 (град/с)
    "dt": 0.01,         # шаг по времени, с
    "t_max": 20.0,      # время симуляции, с
}


# ---------- Механика ----------
def derivatives(state, t, params):
    """Возвращает производную состояния [theta1, omega1, theta2, omega2].

    Из уравнений Лагранжа получается система
        M * [a1, a2] = [f1, f2],
    где a1, a2 - угловые ускорения. Её решаем через np.linalg.solve.
    """
    th1, w1, th2, w2 = state
    m1, m2 = params["m1"], params["m2"]
    L1, L2 = params["L1"], params["L2"]
    g = params["g"]
    c1, c2 = params["c1"], params["c2"]

    delta = th1 - th2

    # матрица масс
    M = np.array([
        [(m1 + m2) * L1**2, m2 * L1 * L2 * np.cos(delta)],
        [m2 * L1 * L2 * np.cos(delta), m2 * L2**2],
    ])

    # правая часть (центробежные силы, гравитация, трение)
    f1 = (-m2 * L1 * L2 * w2**2 * np.sin(delta)
          - (m1 + m2) * g * L1 * np.sin(th1)
          - c1 * w1 + c2 * (w2 - w1))
    f2 = (m2 * L1 * L2 * w1**2 * np.sin(delta)
          - m2 * g * L2 * np.sin(th2)
          - c2 * (w2 - w1))

    a1, a2 = np.linalg.solve(M, [f1, f2])
    return np.array([w1, a1, w2, a2])


def energy(state, params):
    """Полная энергия T + V."""
    th1, w1, th2, w2 = state
    m1, m2 = params["m1"], params["m2"]
    L1, L2 = params["L1"], params["L2"]
    g = params["g"]

    T = (0.5 * (m1 + m2) * L1**2 * w1**2
         + 0.5 * m2 * L2**2 * w2**2
         + m2 * L1 * L2 * w1 * w2 * np.cos(th1 - th2))
    V = -(m1 + m2) * g * L1 * np.cos(th1) - m2 * g * L2 * np.cos(th2)
    return T + V


# ---------- Интегратор RK4 ----------
def rk4_step(state, t, h, params):
    k1 = derivatives(state, t, params)
    k2 = derivatives(state + 0.5 * h * k1, t + 0.5 * h, params)
    k3 = derivatives(state + 0.5 * h * k2, t + 0.5 * h, params)
    k4 = derivatives(state + h * k3, t + h, params)
    return state + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate(state0, t_max, dt, params):
    t = np.arange(0, t_max + dt / 2, dt)
    states = np.zeros((len(t), 4))
    states[0] = state0
    for i in range(1, len(t)):
        states[i] = rk4_step(states[i - 1], t[i - 1], dt, params)
    return t, states


# ---------- Проверка: предел простого маятника ----------
def check_simple_pendulum(params):
    """Лёгкий второй груз + малый угол -> первый стержень ведёт себя как
    обычный маятник с периодом T = 2*pi*sqrt(L/g)."""
    p = dict(params)
    p["m2"] = 1e-9      # очень лёгкий второй груз
    p["c1"] = 0.0
    p["c2"] = 0.0
    a = np.radians(1.0)  # маленький угол

    t, s = simulate([a, 0.0, a, 0.0], 60.0, 0.005, p)

    # ищем моменты, когда угол проходит через 0 снизу вверх
    th = s[:, 0]
    crossings = []
    for i in range(len(th) - 1):
        if th[i] < 0 and th[i + 1] >= 0:
            # линейная интерполяция момента пересечения нуля
            tc = t[i] - th[i] * (t[i + 1] - t[i]) / (th[i + 1] - th[i])
            crossings.append(tc)

    T_measured = np.mean(np.diff(crossings))
    T_theory = 2 * np.pi * np.sqrt(p["L1"] / p["g"])
    T_theory *= 1 + a**2 / 16   # поправка на конечную амплитуду
    return T_measured, T_theory


