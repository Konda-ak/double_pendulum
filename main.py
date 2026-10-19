import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from mechanics import config, simulate, energy, check_simple_pendulum


def main():
    # 1. Проверка: предел простого маятника
    T_meas, T_th = check_simple_pendulum(config)
    print("Проверка простого маятника:")
    print(f"  период из симуляции: {T_meas:.6f} с")
    print(f"  период по формуле:   {T_th:.6f} с")
    print(f"  разница: {abs(T_meas - T_th) / T_th * 100:.5f} %")

    # 2. Основная симуляция
    state0 = np.radians([config["theta1"], config["omega1"],
                         config["theta2"], config["omega2"]])
    t, states = simulate(state0, config["t_max"], config["dt"], config)

    th1, th2 = states[:, 0], states[:, 2]
    L1, L2 = config["L1"], config["L2"]
    x1 = L1 * np.sin(th1)
    y1 = -L1 * np.cos(th1)
    x2 = x1 + L2 * np.sin(th2)
    y2 = y1 - L2 * np.cos(th2)

    E = np.array([energy(s, config) for s in states])
    print(f"Максимальное изменение энергии за {config['t_max']:g} с: "
          f"{np.max(np.abs(E - E[0])):.2e} Дж")

    # 3. Графики: маятник слева, энергия справа
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 2, 1)
    ax_e = fig.add_subplot(2, 2, 2)
    ax_d = fig.add_subplot(2, 2, 4)

    R = L1 + L2 + 0.2
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    ax.set_title("Двойной маятник")
    line, = ax.plot([], [], "o-", lw=2, color="black", markersize=8)
    trace, = ax.plot([], [], "-", lw=1, color="tab:red", alpha=0.7)
    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

    ax_e.plot(t, E)
    ax_e.set_title("Полная энергия")
    ax_e.set_ylabel("E, Дж")
    ax_e.grid(alpha=0.3)

    ax_d.plot(t, E - E[0], color="tab:red")
    ax_d.set_title("Отклонение энергии от начальной")
    ax_d.set_xlabel("t, с")
    ax_d.set_ylabel("ΔE, Дж")
    ax_d.grid(alpha=0.3)

    trail = 300  # длина следа
    step = 2     # пропуск кадров

    def update(frame):
        i = frame * step
        line.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
        j = max(0, i - trail)
        trace.set_data(x2[j:i], y2[j:i])
        time_text.set_text(f"t = {t[i]:.1f} с")
        return line, trace, time_text

    ani = FuncAnimation(fig, update, frames=len(t) // step,
                        interval=config["dt"] * step * 1000, blit=True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
