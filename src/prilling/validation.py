"""Caso de validación: torre industrial de Torkashvand et al. (2019).

Torkashvand, B., Gilassi, S., Moghadam, R. (2019). Analysis of heat and
moisture transfer during drying of urea particles. Chem. Ind. Chem. Eng. Q.
25(2), 171-182. https://doi.org/10.2298/CICEQ180827030T  (Tabla 1)
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from .psd import RosinRammler
from .tower import TowerConfig, simulate_tower

# Tabla 1 del artículo: observado en planta y calculado por los autores [°C]
PLANT = pd.DataFrame(
    [
        ("Temperatura del producto (urea)", None, 60.0, 64.5),
        ("Aire de salida (cabeza)", 79.9, 65.0, 69.4),
        ("Aire a 60 m del fondo (C)", 60.0, 61.0, 63.2),
        ("Aire a 40 m del fondo (B)", 40.0, 55.0, 58.3),
        ("Aire a 20 m del fondo (A)", 20.0, 42.0, 46.1),
    ],
    columns=["variable", "altura_desde_fondo_m", "planta_C", "torkashvand_C"],
)


def n_rr_default() -> float:
    return 6.0


def plant_config(**overrides) -> TowerConfig:
    """Condiciones de la Tabla 1. Supuestos propios marcados en docs/modelo.md §5.

    PSD: el artículo no la informa. Se interpreta "0.85 mm" como RADIO medio
    (las figuras usan r = 0.00085 m), es decir d50 = 1.70 mm, con una
    Rosin-Rammler de n = 6 [SUPUESTO].
    """
    n = n_rr_default()
    d50 = 1.70e-3
    d63 = d50 / np.log(2.0) ** (1.0 / n)
    cfg = TowerConfig(
        H=79.9, D=17.0, m_melt=75e3 / 3600, m_air=775e3 / 3600,
        T_air_in=35 + 273.15, T_melt_in=140 + 273.15, Y_water=0.037,
        psd=RosinRammler(d63=d63, n=n),
    )
    return replace(cfg, **overrides)


def model_values(res) -> list[float]:
    vals = [res.T_product - 273.15]
    for h in PLANT["altura_desde_fondo_m"].iloc[1:]:
        vals.append(res.T_air_at_height_from_bottom(float(h)) - 273.15)
    return vals


def compare(res) -> pd.DataFrame:
    df = PLANT.copy()
    df["modelo_C"] = model_values(res)
    df["error_modelo_C"] = df["modelo_C"] - df["planta_C"]
    df["error_torkashvand_C"] = df["torkashvand_C"] - df["planta_C"]
    return df


def rmse(df: pd.DataFrame, col: str = "error_modelo_C") -> float:
    return float(np.sqrt(np.mean(df[col] ** 2)))


def calibrate_h_factor(base: TowerConfig | None = None, bounds=(0.4, 1.5)) -> tuple[float, pd.DataFrame]:
    """Ajusta h_factor minimizando el RMSE contra las 5 mediciones de planta."""
    base = base or plant_config()

    def obj(f):
        return rmse(compare(simulate_tower(replace(base, h_factor=f))))

    opt = minimize_scalar(obj, bounds=bounds, method="bounded", options={"xatol": 0.01})
    return float(opt.x), compare(simulate_tower(replace(base, h_factor=opt.x)))


def plant_energy_balance() -> dict:
    """¿Cierran los datos de planta? Balance con las propiedades del modelo."""
    from .properties import UREA
    from .tower import _air_enthalpy_funcs

    cfg = plant_config()
    h_air, _ = _air_enthalpy_funcs(cfg.Y_water)
    Q_melt = cfg.m_melt * float(UREA.enthalpy(cfg.T_melt_in) - UREA.enthalpy(60 + 273.15))
    Q_air = cfg.m_air * (h_air(65 + 273.15) - h_air(35 + 273.15))
    return {"Q_cedido_urea_kW": Q_melt / 1e3, "Q_ganado_aire_kW": Q_air / 1e3,
            "descuadre_%": 100 * (Q_melt - Q_air) / Q_melt}
