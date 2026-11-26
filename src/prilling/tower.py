"""Torre de prilling en contracorriente: gotas (lagrangiano) + aire (balance 1D).

Acoplamiento por iteración de punto fijo sobre el perfil de temperatura del
aire (docs/modelo.md §4):

  h_aire(T_a(z)) = h_aire(T_in) + (m_fundido/m_aire) * sum_j w_j [H_j(z) - H_j(H)]

que es el balance de energía en el volumen de control entre la altura z y el
fondo, con pared adiabática.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .droplet import DropletResult, simulate_droplet
from .properties import CP_VAPOR, UREA, Material, air_table
from .psd import RosinRammler, solid_diameter


@dataclass
class TowerConfig:
    H: float = 79.9               # m, altura de caída libre
    D: float = 17.0               # m, diámetro de la torre
    m_melt: float = 75e3 / 3600   # kg/s
    m_air: float = 775e3 / 3600   # kg/s (aire húmedo)
    T_air_in: float = 35 + 273.15  # K
    T_melt_in: float = 140 + 273.15  # K
    Y_water: float = 0.037        # kg agua / kg aire seco (solo afecta cp)
    psd: RosinRammler = field(default_factory=lambda: RosinRammler(d63=1.9e-3, n=6.0))
    n_classes: int = 10
    material: Material = UREA
    nusselt: str = "Ranz-Marshall"
    radiation: bool = False
    n_shells: int = 16
    v0: float = 0.0
    h_factor: float = 1.0
    rtol: float = 1e-5
    n_z: int = 81
    max_iter: int = 40
    tol_K: float = 0.05
    relax: float = 0.7

    def __post_init__(self) -> None:
        for name in ("H", "D", "m_melt", "m_air", "T_air_in", "T_melt_in", "h_factor"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} debe ser positivo (recibido {getattr(self, name)})")
        if self.Y_water < 0:
            raise ValueError("Y_water no puede ser negativa")
        if self.T_melt_in <= self.material.T_melt:
            raise ValueError("El fundido debe entrar sobre su punto de fusión")
        if self.T_air_in >= self.material.T_melt:
            raise ValueError("El aire debe entrar bajo el punto de fusión del material")

    @property
    def area(self) -> float:
        return np.pi * self.D**2 / 4.0


@dataclass
class TowerResult:
    config: TowerConfig
    z: np.ndarray                 # m, distancia desde la cabeza
    T_air: np.ndarray             # K
    u_air: np.ndarray             # m/s
    classes: pd.DataFrame
    droplets: list[DropletResult]
    iterations: int
    converged: bool
    energy_balance_error: float   # relativo
    T_product: float              # K, temperatura de mezcla del producto
    T_air_out: float              # K
    offspec_liquid_fraction: float  # fracción másica con líquido residual al fondo (desde d_crit)
    d_crit: float = float("nan")    # m, gota más grande que solidifica en H

    def T_air_at_height_from_bottom(self, h_from_bottom: float) -> float:
        return float(np.interp(self.config.H - h_from_bottom, self.z, self.T_air))

    def summary(self) -> dict:
        return {
            "T_producto_C": self.T_product - 273.15,
            "T_aire_salida_C": self.T_air_out - 273.15,
            "masa_con_liquido_residual_%": 100 * self.offspec_liquid_fraction,
            "d_critico_gota_mm": self.d_crit * 1e3,
            "iteraciones": self.iterations,
            "convergio": self.converged,
            "error_balance_energia_%": 100 * self.energy_balance_error,
        }


def _air_enthalpy_funcs(Y: float):
    """h(T) y T(h) del aire húmedo por kg de mezcla, integrando cp de la tabla."""
    air = air_table()
    cp_mix = (air.cp + Y * CP_VAPOR) / (1.0 + Y)
    h = np.concatenate(([0.0], np.cumsum(0.5 * (cp_mix[1:] + cp_mix[:-1]) * np.diff(air.T))))
    return (lambda T: np.interp(T, air.T, h)), (lambda hh: np.interp(hh, h, air.T))


def simulate_tower(cfg: TowerConfig, progress=None) -> TowerResult:
    mat = cfg.material
    air = air_table()
    h_air, T_of_h = _air_enthalpy_funcs(cfg.Y_water)
    z = np.linspace(0.0, cfg.H, cfg.n_z)
    d_cls, w_cls, _ = cfg.psd.classes(cfg.n_classes)
    G_air = cfg.m_air / cfg.area  # kg/(m2 s)

    # Estimación inicial: producto a T_in + 20 K, perfil lineal del aire
    dh_melt = mat.enthalpy(cfg.T_melt_in) - mat.enthalpy(cfg.T_air_in + 20.0)
    T_out0 = T_of_h(h_air(cfg.T_air_in) + cfg.m_melt * dh_melt / cfg.m_air)
    T_air = np.interp(z, [0.0, cfg.H], [T_out0, cfg.T_air_in])

    converged = False
    for it in range(1, cfg.max_iter + 1):
        T_prof = T_air.copy()
        u_prof = G_air / np.array([air.at(T)[0] for T in T_prof])
        def Ta_f(zz, T_prof=T_prof):
            return float(np.interp(zz, z, T_prof))

        def ua_f(zz, u_prof=u_prof):
            return float(np.interp(zz, z, u_prof))

        drops, Hz = [], np.empty((len(d_cls), len(z)))
        for j, d in enumerate(d_cls):
            r = simulate_droplet(d, mat, cfg.H, cfg.T_melt_in, Ta_f, ua_f,
                                 n_shells=cfg.n_shells, nusselt=cfg.nusselt,
                                 radiation=cfg.radiation, v0=cfg.v0, h_factor=cfg.h_factor, n_out=300, rtol=cfg.rtol)
            drops.append(r)
            Hz[j] = np.interp(z, r.z, r.h_mean)
        H_mix = (w_cls[:, None] * Hz).sum(axis=0)
        T_new = T_of_h(h_air(cfg.T_air_in) + cfg.m_melt / cfg.m_air * (H_mix - H_mix[-1]))
        change = np.max(np.abs(T_new - T_air))
        T_air = T_air + cfg.relax * (T_new - T_air)
        if progress:
            progress(it, change)
        if change < cfg.tol_K:
            converged = True
            T_air = T_new
            break

    # Resultados por clase
    rows = []
    for d, w, r in zip(d_cls, w_cls, drops, strict=True):
        e = r.exit
        rows.append({
            "d_gota_mm": d * 1e3,
            "d_prill_mm": float(solid_diameter(d, mat.rho_l, mat.rho_s)) * 1e3,
            "fraccion_masica": w,
            "t_residencia_s": e["t_residencia_s"],
            "v_salida_m_s": e["v_salida_m_s"],
            "T_media_salida_C": e["T_media_salida_C"],
            "T_superficie_salida_C": e["T_superficie_salida_C"],
            "T_centro_salida_C": e["T_centro_salida_C"],
            "fraccion_liquida_salida": e["fraccion_liquida_salida"],
            "z_solidificacion_m": e["z_solidificacion_m"],
            "solidifica": e["fraccion_liquida_salida"] < 1e-4,
        })
    df = pd.DataFrame(rows)

    h_prod = float((w_cls * Hz[:, -1]).sum())
    T_product = float(mat.temperature(h_prod))
    T_air_out = float(T_air[0])
    Q_melt = cfg.m_melt * (float(mat.enthalpy(cfg.T_melt_in)) - h_prod)
    Q_air = cfg.m_air * (h_air(T_air_out) - h_air(cfg.T_air_in))
    err = abs(Q_air - Q_melt) / Q_melt
    d_crit = critical_diameter(df, cfg.H)
    offspec = float(1.0 - cfg.psd.cdf(d_crit)) if np.isfinite(d_crit) else 0.0
    u_air = G_air / np.array([air.at(T)[0] for T in T_air])
    return TowerResult(cfg, z, T_air, u_air, df, drops, it, converged, err,
                       T_product, T_air_out, offspec, d_crit)


def critical_diameter(df: pd.DataFrame, H: float) -> float:
    """Diámetro de gota [m] que solidifica justo al llegar al fondo.

    Ajuste potencial z_solid = a d^b sobre las (hasta) 3 clases sólidas más
    grandes, evaluado en z_solid = H. Evita que la fracción blanda salte de
    clase en clase con una PSD discretizada en pocas clases. Devuelve inf si no
    hay clases suficientes para extrapolar y todas solidifican.
    """
    ok = df[df["solidifica"] & np.isfinite(df["z_solidificacion_m"])]
    if len(ok) < 2:
        return 0.0 if len(ok) == 0 else float("inf")
    tail = ok.tail(3)
    b, ln_a = np.polyfit(np.log(tail["d_gota_mm"] * 1e-3), np.log(tail["z_solidificacion_m"]), 1)
    return float(np.exp((np.log(H) - ln_a) / b))
