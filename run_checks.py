"""Проверки корректности симуляции. Запуск: python run_checks.py

1. RK4 на уравнении с точным решением (y' = -y) и сравнение с scipy
2. Порядок RK4 по ошибке энергии двойного маятника
3. Уравнения движения против уравнений Эйлера-Лагранжа (численно)
4. Нормальные моды: периоды и формы колебаний
5. Трение: dE/dt равно мощности сил трения
6. Предел простого маятника
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.special import ellipk

from config import load_config
from integrator import integrate, rk4_step
from physics import (derivatives, energy_scale, friction_power, initial_state,
                     linear_matrices, linearized_system, make_rhs, mass_matrix,
                     normal_modes, numeric_jacobian, total_energy)

BASE = "configs/chaos.json"      # базовые параметры системы (m, L, g)
results = []


def report(name, ok, info):
    results.append(ok)
    print(f"[{'OK  ' if ok else 'FAIL'}] {name}: {info}")


def measure_period(t, x):
    """Период по переходам x через ноль снизу вверх (с интерполяцией)."""
    i = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    tc = t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i])
    return np.mean(np.diff(tc))


# ---------------------------------------------------------------- 1
def check_rk4_exact():
    """y' = -y, y(0) = 1: точное решение exp(-t)."""
    f = lambda y, t: -y
    errs = []
    steps = [0.4, 0.2, 0.1, 0.05]
    for h in steps:
        t, y = integrate(f, np.array([1.0]), 2.0, h)
        errs.append(abs(y[-1, 0] - np.exp(-2.0)))
    orders = [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    report("RK4 на y'=-y: порядок", abs(orders[-1] - 4) < 0.2,
           "отношения ошибок при h/2: " + ", ".join(f"{2**o:.1f}" for o in orders)
           + " (ожидается ~16)")

    # гармонический осциллятор y'' = -y: сравнение с точным решением cos(t)
    osc = lambda y, t: np.array([y[1], -y[0]])
    t, y = integrate(osc, np.array([1.0, 0.0]), 20.0, 0.01)
    err = np.max(np.abs(y[:, 0] - np.cos(t)))
    report("RK4 на осцилляторе", err < 1e-8, f"max|y - cos t| = {err:.1e} за 20 с")


def check_against_scipy(cfg):
    """Эталон: другой метод (DOP853 с жёсткими допусками), регулярный режим."""
    p = cfg
    y0 = np.radians([40.0, 0.0, -20.0, 0.0])
    T = 5.0
    t, y = integrate(make_rhs(p), y0, T, 0.001)
    ref = solve_ivp(lambda tt, yy: derivatives(yy, tt, p), (0, T), y0,
                    method="DOP853", rtol=1e-12, atol=1e-12)
    err = np.max(np.abs(y[-1] - ref.y[:, -1]))
    report("RK4 против scipy DOP853", err < 1e-7,
           f"max отличие конечного состояния {err:.1e} (T = {T:g} с)")


def check_state_order(cfg):
    """Ошибка состояния двойного маятника относительно DOP853 (rtol 1e-13)."""
    y0 = np.radians([2.0, 0.0, 3.0, 0.0])
    T = 20.0
    te = np.arange(0, T + 1e-9, 0.04)
    ref = solve_ivp(lambda tt, yy: derivatives(yy, tt, cfg), (0, T), y0, method="DOP853",
                    rtol=1e-13, atol=1e-13, t_eval=te).y.T
    errs = []
    dts = [0.04, 0.02, 0.01, 0.005]
    for dt in dts:
        t, y = integrate(make_rhs(cfg), y0, T, dt)
        k = int(round(0.04 / dt))
        errs.append(np.max(np.abs(y[::k] - ref)))
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    report("Порядок RK4 по состоянию (малые колебания)", abs(ratios[-1] - 16) < 2,
           "max ошибка состояния при dt=" + ", ".join(map(str, dts)) + ": "
           + ", ".join(f"{e:.1e}" for e in errs) + "; отношения при dt/2: "
           + ", ".join(f"{r:.1f}" for r in ratios) + " (ожидается 16)")


# ---------------------------------------------------------------- 2
def energy_error(cfg, y0, dt, T):
    t, y = integrate(make_rhs(cfg), y0, T, dt)
    E = total_energy(y, cfg)
    return np.max(np.abs(E - E[0])) / energy_scale(cfg)


def check_energy_convergence(cfg):
    dts = np.array([0.04, 0.02, 0.01, 0.005, 0.0025])
    cases = {
        "регулярный режим (40°, -20°)": np.radians([40.0, 0.0, -20.0, 0.0]),
        "хаотический режим (120°, -10°)": initial_state(cfg),
    }
    plt.figure(figsize=(6, 4.5))
    for name, y0 in cases.items():
        errs = np.array([energy_error(cfg, y0, dt, 10.0) for dt in dts])
        ratios = errs[:-1] / errs[1:]
        plt.loglog(dts, errs, "o-", label=name)
        info = ", ".join(f"{r:.0f}" for r in ratios)
        # проверяем порядок на мелких шагах, где асимптотика уже наступила
        ok = ratios[-1] > 12
        report(f"Ошибка энергии, {name}", ok,
               f"E_err/E_scale при dt={dts[0]}..{dts[-1]}: "
               + ", ".join(f"{e:.1e}" for e in errs) + f"; отношения при dt/2: {info}")
    plt.loglog(dts, 1e-3 * (dts / dts[0]) ** 4 * 1e-1, "k--", label="наклон 4")
    plt.xlabel("dt, с")
    plt.ylabel("max|E - E0| / E_scale")
    plt.title("Сходимость RK4 по энергии")
    plt.grid(alpha=0.3, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/convergence.png", dpi=120)
    plt.close()


# ---------------------------------------------------------------- 3
def lagrangian_pieces(q, qd, p):
    """T, V, R записаны независимо от physics.py - через декартовы координаты."""
    th1, th2 = q
    w1, w2 = qd
    vx1, vy1 = p["L1"] * np.cos(th1) * w1, p["L1"] * np.sin(th1) * w1
    vx2 = vx1 + p["L2"] * np.cos(th2) * w2
    vy2 = vy1 + p["L2"] * np.sin(th2) * w2
    T = 0.5 * p["m1"] * (vx1 ** 2 + vy1 ** 2) + 0.5 * p["m2"] * (vx2 ** 2 + vy2 ** 2)
    y1 = -p["L1"] * np.cos(th1)
    y2 = y1 - p["L2"] * np.cos(th2)
    V = p["m1"] * p["g"] * y1 + p["m2"] * p["g"] * y2
    R = 0.5 * p["c1"] * w1 ** 2 + 0.5 * p["c2"] * (w2 - w1) ** 2
    return T - V, R


def check_euler_lagrange(cfg):
    """d/dt(dL/dq') - dL/dq + dR/dq' = 0 на случайных состояниях (с трением)."""
    p = dict(cfg, c1=0.3, c2=0.2)
    rng = np.random.default_rng(1)
    eps = 1e-4

    def grad(fun, x):
        g = np.zeros(2)
        for i in range(2):
            d = np.zeros(2)
            d[i] = eps
            g[i] = (fun(x + d) - fun(x - d)) / (2 * eps)
        return g

    def p_mom(s):      # dL/dq' в состоянии s = [th1, w1, th2, w2]
        q, qd = np.array([s[0], s[2]]), np.array([s[1], s[3]])
        return grad(lambda v: lagrangian_pieces(q, v, p)[0], qd)

    worst = 0.0
    for _ in range(30):
        s = rng.uniform(-3, 3, 4)
        q, qd = np.array([s[0], s[2]]), np.array([s[1], s[3]])
        f = derivatives(s, 0.0, p)                      # [w1, a1, w2, a2]
        h = 1e-4
        dp_dt = (p_mom(s + h * f) - p_mom(s - h * f)) / (2 * h)   # d/dt по траектории
        dL_dq = grad(lambda v: lagrangian_pieces(v, qd, p)[0], q)
        dR_dqd = grad(lambda v: lagrangian_pieces(q, v, p)[1], qd)
        resid = dp_dt - dL_dq + dR_dqd
        scale = np.max(np.abs(dL_dq)) + 1.0
        worst = max(worst, np.max(np.abs(resid)) / scale)
    report("Уравнения Эйлера-Лагранжа (численно, с трением)", worst < 1e-5,
           f"max невязка (отн.) = {worst:.1e} на 30 случайных состояниях")


def check_linearization(cfg):
    """Численный якобиан derivatives в нижнем положении против аналитических A, B."""
    p = dict(cfg, m2=0.6, L2=0.8, c1=0.2, c2=0.1)
    A_an, B_an = linearized_system(p)
    A_num, B_num = numeric_jacobian(p)
    err = max(np.max(np.abs(A_an - A_num)), np.max(np.abs(B_an - B_num)))
    report("Линеаризация: численный якобиан против A, B из DERIVATION.md", err < 1e-6,
           f"max отличие = {err:.1e}")


# ---------------------------------------------------------------- 4
def check_normal_modes(cfg):
    sets = {
        "m1=m2, L1=L2": dict(cfg),
        "m1=1, m2=0.5, L1=1, L2=0.7": dict(cfg, m1=1.0, m2=0.5, L1=1.0, L2=0.7),
    }
    amp = 1e-3
    for label, p in sets.items():
        p = dict(p, c1=0.0, c2=0.0)
        omegas, V = normal_modes(p)
        M0, _ = linear_matrices(p)
        for i in range(2):
            v = V[:, i] / np.max(np.abs(V[:, i]))
            y0 = np.array([amp * v[0], 0.0, amp * v[1], 0.0])
            t, y = integrate(make_rhs(p), y0, 60.0, 0.002)
            comp = 0 if abs(v[0]) >= abs(v[1]) else 2
            T_sim = measure_period(t, y[:, comp])
            T_th = 2 * np.pi / omegas[i]
            rel = abs(T_sim - T_th) / T_th

            # чистота моды: вклад другой моды (в M0-ортогональных координатах)
            q = y[:, [0, 2]].T
            eta = V.T @ M0 @ q
            leak = np.max(np.abs(eta[1 - i])) / np.max(np.abs(eta[i]))
            report(f"Нормальная мода {i + 1} ({label})", rel < 1e-5 and leak < 1e-2,
                   f"T_sim = {T_sim:.6f} с, T_теор = {T_th:.6f} с, "
                   f"отн. ошибка {rel:.1e}, примесь другой моды {leak:.1e}, "
                   f"θ2/θ1 = {v[1] / v[0]:.3f}")

    # смесь мод: наблюдается биение, обе частоты присутствуют в спектре
    p = dict(cfg, c1=0.0, c2=0.0)
    omegas, _ = normal_modes(p)
    t, y = integrate(make_rhs(p), np.radians([1.0, 0.0, 0.0, 0.0]), 80.0, 0.005)
    nfft = 16 * len(t)      # дополнение нулями для точного положения пиков
    spec = np.abs(np.fft.rfft(y[:, 0] * np.hanning(len(t)), n=nfft))
    freq = 2 * np.pi * np.fft.rfftfreq(nfft, 0.005)
    top = []
    for _ in range(2):          # два самых больших пика
        k = np.argmax(spec)
        top.append(freq[k])
        spec[max(0, k - 600):k + 600] = 0
    top = sorted(top)
    ok = all(abs(a - b) < 0.02 for a, b in zip(top, omegas))
    report("Спектр свободных колебаний", ok,
           f"пики FFT: {top[0]:.3f}, {top[1]:.3f} рад/с; теория: "
           f"{omegas[0]:.3f}, {omegas[1]:.3f} рад/с")


# ---------------------------------------------------------------- 5
def check_friction(cfg_path):
    cfg = load_config(cfg_path)
    y0 = initial_state(cfg)
    dt = 0.002
    t, y = integrate(make_rhs(cfg), y0, cfg["t_max"], dt)
    E = total_energy(y, cfg)
    P = friction_power(y, cfg)
    scale = energy_scale(cfg)

    # (а) мгновенная скорость изменения энергии вдоль траектории
    worst = 0.0
    for s in y[::50]:
        f = derivatives(s, 0.0, cfg)
        h = 1e-6
        dE = (total_energy(s + h * f, cfg) - total_energy(s - h * f, cfg)) / (2 * h)
        worst = max(worst, abs(dE[0] - friction_power(s, cfg)[0]))
    report("Трение: dE/dt = мощность трения (мгновенно)", worst / np.max(np.abs(P)) < 1e-6,
           f"max|dE/dt - P| / max|P| = {worst / np.max(np.abs(P)):.1e}")

    # (б) интеграл мощности по времени против изменения энергии в симуляции
    work = cumulative_trapezoid(P, t, initial=0.0)
    diff = np.max(np.abs((E - E[0]) - work)) / scale
    report("Трение: E(t) - E(0) = интеграл P dt", diff < 1e-5,
           f"max отличие / E_scale = {diff:.1e}; "
           f"потеряно {(E[0] - E[-1]) / scale * 100:.1f}% E_scale за {cfg['t_max']:g} с")
    report("Трение: энергия не растёт", np.all(np.diff(E) < 1e-9 * scale),
           f"max прирост за шаг = {np.max(np.diff(E)) / scale:.1e}")


# ---------------------------------------------------------------- 6
def check_simple_pendulum(cfg):
    p = dict(cfg, m2=1e-9, c1=0.0, c2=0.0)
    a = np.radians(1.0)
    t, y = integrate(make_rhs(p), np.array([a, 0.0, a, 0.0]), 60.0, 0.005)
    T_sim = measure_period(t, y[:, 0])
    T_exact = 4 * np.sqrt(p["L1"] / p["g"]) * ellipk(np.sin(a / 2) ** 2)
    rel = abs(T_sim - T_exact) / T_exact
    report("Предел простого маятника (m2 -> 0, 1°)", rel < 1e-6,
           f"T_sim = {T_sim:.6f} с, T_точн = {T_exact:.6f} с, "
           f"2π√(L/g) = {2 * np.pi * np.sqrt(p['L1'] / p['g']):.6f} с")


def main():
    cfg = load_config(BASE)
    check_rk4_exact()
    check_against_scipy(cfg)
    check_state_order(cfg)
    check_energy_convergence(cfg)
    check_euler_lagrange(cfg)
    check_linearization(cfg)
    check_normal_modes(cfg)
    check_friction("configs/damped.json")
    check_simple_pendulum(cfg)
    print(f"\nПройдено {sum(results)} из {len(results)}")
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main()
