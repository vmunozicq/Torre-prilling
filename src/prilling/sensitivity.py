"""Análisis de sensibilidad: tornado (uno a la vez) y mapa de altura mínima."""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from .droplet import simulate_droplet
from .psd import RosinRammler
from .tower import TowerConfig, simulate_tower

# Parámetro -> (descripción, función que aplica un factor multiplicativo)
PARAMS = {
    "d63 (tamaño de gota)": lambda c, f: replace(c, psd=RosinRammler(c.psd.d63 * f, c.psd.n)),
    "Flujo de aire": lambda c, f: replace(c, m_air=c.m_air * f),
    "Flujo de fundido": lambda c, f: replace(c, m_melt=c.m_melt * f),
    "Altura de caída": lambda c, f: replace(c, H=c.H * f),
    "T aire entrada (°C)": lambda c, f: replace(c, T_air_in=273.15 + (c.T_air_in - 273.15) * f),
    "k sólido urea": lambda c, f: replace(c, material=c.material.with_(k_s=c.material.k_s * f)),
    "Calor latente": lambda c, f: replace(c, material=c.material.with_(L_fus=c.material.L_fus * f)),
    "Coef. convectivo (h_factor)": lambda c, f: replace(c, h_factor=c.h_factor * f),
}


def tornado(cfg: TowerConfig, delta: float = 0.2, params: dict | None = None, progress=None) -> pd.DataFrame:
    params = params or PARAMS
    base = simulate_tower(cfg)
    rows = []
    for i, (name, apply) in enumerate(params.items()):
        out = {"parametro": name}
        for tag, f in (("bajo", 1 - delta), ("alto", 1 + delta)):
            r = simulate_tower(apply(cfg, f))
            out[f"T_producto_{tag}_C"] = r.T_product - 273.15
            out[f"liquido_residual_{tag}_%"] = 100 * r.offspec_liquid_fraction
            out[f"T_aire_salida_{tag}_C"] = r.T_air_out - 273.15
        rows.append(out)
        if progress:
            progress(i + 1, len(params))
    df = pd.DataFrame(rows)
    df["T_producto_base_C"] = base.T_product - 273.15
    df["rango_T_producto_C"] = (df["T_producto_alto_C"] - df["T_producto_bajo_C"]).abs()
    return df.sort_values("rango_T_producto_C", ascending=False).reset_index(drop=True)


def min_height_map(material, diameters_mm, T_air_C, u_air: float = 0.9, T_melt_in_C: float = 140.0,
                   nusselt: str = "Ranz-Marshall", h_factor: float = 1.0, n_shells: int = 12) -> pd.DataFrame:
    """Altura de caída necesaria para solidificar por completo una gota aislada,
    en aire a temperatura constante con velocidad ascendente u_air."""
    rows = []
    for Ta in T_air_C:
        for d in diameters_mm:
            r = simulate_droplet(d * 1e-3, material, 400.0, T_melt_in_C + 273.15, Ta + 273.15, u_air,
                                 n_shells=n_shells, nusselt=nusselt, h_factor=h_factor, n_out=300, rtol=1e-5)
            rows.append({"d_mm": d, "T_aire_C": Ta, "z_solidificacion_m": r.z_solid, "t_solidificacion_s": r.t_solid})
    return pd.DataFrame(rows)


def design_curve(material, H: float, T_air_C: float, u_air: float = 0.9, **kw) -> float:
    """Diámetro máximo [mm] que solidifica en una altura H (bisección sobre d)."""
    lo, hi = 0.3, 6.0
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        z = min_height_map(material, [mid], [T_air_C], u_air, **kw)["z_solidificacion_m"].iloc[0]
        if np.isfinite(z) and z <= H:
            lo = mid
        else:
            hi = mid
    return lo
