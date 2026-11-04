"""Запуск симуляции: python main.py [конфиг.json] [--save картинка.png] [--no-show]

Примеры:
    python main.py                          # configs/chaos.json
    python main.py configs/small.json       # малые колебания
    python main.py configs/damped.json      # затухание с трением
    python main.py configs/chaos.json --save figures/chaos.png --no-show
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np

from config import load_config
from integrator import integrate
from physics import energy_scale, initial_state, make_rhs, total_energy
from visualization import animate, build_figure


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config", nargs="?", default="configs/chaos.json")
    ap.add_argument("--save", help="сохранить итоговую картинку (png)")
    ap.add_argument("--no-show", action="store_true", help="не открывать окно")
    args = ap.parse_args()

    cfg = load_config(args.config)
    f = make_rhs(cfg)                      # правая часть f(y, t) для интегратора
    t, y = integrate(f, initial_state(cfg), cfg["t_max"], cfg["dt"])

    E = total_energy(y, cfg)
    print(f"Конфиг: {args.config}; dt = {cfg['dt']} с, t_max = {cfg['t_max']} с")
    print(f"Энергия: E(0) = {E[0]:.4f} Дж, E({cfg['t_max']:g}) = {E[-1]:.4f} Дж")
    print(f"max|E - E(0)| / E_scale = {np.max(np.abs(E - E[0])) / energy_scale(cfg):.2e}")

    fig, update, n_frames = build_figure(cfg, t, y)
    fig.suptitle(args.config)
    if args.save:
        update(n_frames - 1)
        fig.savefig(args.save, dpi=110)
        print("Сохранено:", args.save)
    if not args.no_show:
        ani = animate(fig, update, n_frames, cfg["dt"])   # ссылка нужна, иначе GC
        plt.show()


if __name__ == "__main__":
    main()
