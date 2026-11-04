"""Визуализация: анимация маятника, углы от времени, фазовые портреты, энергия."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from physics import energy_scale, positions, total_energy


def wrap_angle(a):
    """Приведение угла к диапазону [-pi, pi]. Только для графиков:
    интегрируется всегда непрерывный угол."""
    return (a + np.pi) % (2 * np.pi) - np.pi


def wrapped_for_plot(a):
    """Приведённый угол, где на скачках через +-pi вставлен NaN,
    чтобы matplotlib не рисовал вертикальные линии через весь график."""
    w = wrap_angle(a)
    jumps = np.where(np.abs(np.diff(w)) > np.pi)[0]
    w = w.copy()
    w[jumps] = np.nan
    return w


def build_figure(cfg, t, y):
    """Рисует всю фигуру. Возвращает (fig, update, n_frames),
    update(i) переставляет маятник и маркеры в момент t[i]."""
    L1, L2 = cfg["L1"], cfg["L2"]
    x1, y1, x2, y2 = [a for a in positions(y, cfg)]
    th1, w1, th2, w2 = y.T
    E = total_energy(y, cfg)

    fig = plt.figure(figsize=(15, 7))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1, 1])
    ax = fig.add_subplot(gs[:, 0])
    ax_ang = fig.add_subplot(gs[0, 1])
    ax_en = fig.add_subplot(gs[1, 1])
    ax_p1 = fig.add_subplot(gs[0, 2])
    ax_p2 = fig.add_subplot(gs[1, 2])

    # --- маятник ---
    R = L1 + L2 + 0.2
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    ax.set_title("Двойной маятник")
    ax.plot(x2, y2, color="tab:red", lw=0.4, alpha=0.25)      # вся траектория
    rod, = ax.plot([], [], "o-", lw=2, color="black", markersize=8)
    trace, = ax.plot([], [], "-", lw=1.2, color="tab:red", alpha=0.8)
    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

    # --- углы от времени (в [-pi, pi]) ---
    ax_ang.plot(t, wrapped_for_plot(th1), label=r"$\theta_1$")
    ax_ang.plot(t, wrapped_for_plot(th2), label=r"$\theta_2$")
    lim_a = np.nanmax(np.abs([wrapped_for_plot(th1), wrapped_for_plot(th2)]))
    if lim_a > 0.8 * np.pi:          # большие углы: показываем весь диапазон [-pi, pi]
        ax_ang.set_ylim(-np.pi - 0.2, np.pi + 0.2)
        ax_ang.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
        ax_ang.set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    else:                            # малые углы: подгоняем масштаб под данные
        ax_ang.set_ylim(-1.3 * lim_a, 1.3 * lim_a)
    ax_ang.set_title("Углы от времени")
    ax_ang.set_ylabel("угол, рад")
    ax_ang.legend(loc="upper right")
    ax_ang.grid(alpha=0.3)
    vline_a = ax_ang.axvline(0, color="gray", lw=0.8)

    # --- энергия ---
    ax_en.plot(t, (E - E[0]) / energy_scale(cfg), color="black")
    ax_en.set_title("Изменение полной энергии")
    ax_en.set_xlabel("t, с")
    ax_en.set_ylabel(r"$(E - E_0)\,/\,E_{scale}$")
    ax_en.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useOffset=False)
    ax_en.grid(alpha=0.3)
    vline_e = ax_en.axvline(0, color="gray", lw=0.8)

    # --- фазовые портреты ---
    dots = []
    for axp, th, w, name in [(ax_p1, th1, w1, "1"), (ax_p2, th2, w2, "2")]:
        axp.plot(wrapped_for_plot(th), np.where(np.isnan(wrapped_for_plot(th)), np.nan, w),
                 lw=0.6, color="tab:blue")
        lim = np.nanmax(np.abs(wrapped_for_plot(th)))
        if lim > 0.8 * np.pi:
            axp.set_xlim(-np.pi - 0.2, np.pi + 0.2)
            axp.set_xticks([-np.pi, 0, np.pi])
            axp.set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
        else:
            axp.set_xlim(-1.3 * lim, 1.3 * lim)
        axp.set_title(f"Фазовый портрет, стержень {name}")
        axp.set_xlabel(rf"$\theta_{name}$, рад")
        axp.set_ylabel(rf"$\omega_{name}$, рад/с")
        axp.grid(alpha=0.3)
        dot, = axp.plot([], [], "o", color="tab:red")
        dots.append(dot)

    trail = int(round(5.0 / (t[1] - t[0])))        # след за последние 5 секунд

    def update(i):
        rod.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
        j = max(0, i - trail)
        trace.set_data(x2[j:i + 1], y2[j:i + 1])
        time_text.set_text(f"t = {t[i]:.2f} с")
        vline_a.set_xdata([t[i], t[i]])
        vline_e.set_xdata([t[i], t[i]])
        dots[0].set_data([wrap_angle(th1[i])], [w1[i]])
        dots[1].set_data([wrap_angle(th2[i])], [w2[i]])
        return [rod, trace, time_text, vline_a, vline_e] + dots

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig, update, len(t)


def animate(fig, update, n_frames, dt, frame_skip=4):
    """Запускает анимацию (кадры пропускаются, чтобы шла примерно в реальном времени)."""
    frames = range(0, n_frames, frame_skip)
    return FuncAnimation(fig, update, frames=frames, interval=dt * frame_skip * 1000,
                         blit=False, repeat=True)
