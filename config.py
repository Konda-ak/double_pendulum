"""Чтение конфигурации. Все параметры лежат в JSON-файле, в коде их нет."""
import json

REQUIRED = ["m1", "m2", "L1", "L2", "g", "c1", "c2",
            "theta1_deg", "theta2_deg", "omega1_deg", "omega2_deg",
            "dt", "t_max"]


def load_config(path):
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    missing = [k for k in REQUIRED if k not in cfg]
    if missing:
        raise KeyError(f"В {path} не хватает параметров: {missing}")
    if cfg["dt"] <= 0 or cfg["t_max"] <= 0:
        raise ValueError("dt и t_max должны быть положительными")
    return cfg
