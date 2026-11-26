"""Etapas aguas abajo de la torre: tamizado, reciclo y enfriador de lecho fluidizado.

Ver docs/modelo.md §7. Diagrama:

  fundido fresco ──► (+ reciclo redisuelto) ──► TORRE ──► TAMIZ ──► ENFRIADOR ──► producto
                          ▲                                 │
                          └──── finos + gruesos/blandos ◄───┘
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from .properties import CP_VAPOR, air_table
from .tower import TowerConfig, TowerResult, simulate_tower


# --------------------------------------------------------------------------
# Tamizado
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class ScreenResult:
    on_spec: float        # fracción másica en especificación (sólida y en rango)
    fines: float          # bajo d_min
    oversize: float       # sobre d_max (sólidos)
    soft: float           # con líquido residual (cualquier tamaño en rango o sobre)
    T_on_spec: float      # K, temperatura de mezcla de la fracción en especificación

    @property
    def offspec(self) -> float:
        return self.fines + self.oversize + self.soft


def screen(res: TowerResult, d_min_mm: float, d_max_mm: float) -> ScreenResult:
    """Tamiz ideal (corte nítido) por diámetro del prill sólido, sobre la PSD continua.

    Los prills con líquido residual se cuentan como fuera de especificación
    aunque estén en rango: se deforman y se pegan en la malla (supuesto P2).
    """
    if not 0 < d_min_mm < d_max_mm:
        raise ValueError("Se requiere 0 < d_min < d_max")
    cfg = res.config
    mat = cfg.material
    df = res.classes
    # Cortes del tamiz (prill sólido) llevados a diámetro de gota fundida
    k = (mat.rho_s / mat.rho_l) ** (1.0 / 3.0)
    F = cfg.psd.cdf
    d_lo, d_hi = d_min_mm * 1e-3 * k, d_max_mm * 1e-3 * k
    d_c = min(res.d_crit, 1e9)
    fines = float(F(min(d_lo, d_c)))
    soft = float(1.0 - F(d_c))
    oversize = float(max(F(d_c) - F(max(d_hi, d_lo)), 0.0))
    on_spec = float(max(F(min(d_hi, d_c)) - F(d_lo), 0.0))
    # Temperatura de mezcla de lo que pasa (aproximada con las clases)
    w = df["fraccion_masica"].to_numpy()
    d = df["d_prill_mm"].to_numpy()
    ok = (d >= d_min_mm) & (d <= d_max_mm) & df["solidifica"].to_numpy()
    h_exit = np.array([r.h_mean[-1] for r in res.droplets])
    T_ok = float(mat.temperature((w[ok] * h_exit[ok]).sum() / w[ok].sum())) if ok.any() else float("nan")
    return ScreenResult(on_spec, fines, oversize, soft, T_ok)


# --------------------------------------------------------------------------
# Reciclo
# --------------------------------------------------------------------------
@dataclass
class RecycleResult:
    m_fresh: float          # kg/s, fundido fresco
    m_tower: float          # kg/s, carga de la torre (fresco + reciclo)
    m_recycle: float        # kg/s
    tower: TowerResult
    screen: ScreenResult
    iterations: int
    converged: bool

    @property
    def recycle_ratio(self) -> float:
        return self.m_recycle / self.m_fresh

    def summary(self) -> dict:
        s = self.screen
        return {
            "fundido_fresco_t_h": self.m_fresh * 3.6,
            "carga_torre_t_h": self.m_tower * 3.6,
            "reciclo_t_h": self.m_recycle * 3.6,
            "razon_reciclo_%": 100 * self.recycle_ratio,
            "en_especificacion_%": 100 * s.on_spec,
            "finos_%": 100 * s.fines,
            "gruesos_%": 100 * s.oversize,
            "blandos_%": 100 * s.soft,
            "T_producto_tamizado_C": s.T_on_spec - 273.15,
            "T_aire_salida_C": self.tower.T_air_out - 273.15,
            "iteraciones": self.iterations,
            "convergio": self.converged,
        }


def solve_recycle(cfg: TowerConfig, d_min_mm: float, d_max_mm: float, max_iter: int = 15,
                  tol: float = 0.002, progress=None) -> RecycleResult:
    """Estado estacionario con reciclo total de lo fuera de especificación.

    El reciclo se redisuelve/refunde y vuelve con el fundido (supuesto P1), así
    que la torre procesa m_torre = m_fresco + m_reciclo con la misma PSD de gota
    y el mismo aire. Balance: m_reciclo = f_fuera(m_torre) * m_torre, resuelto
    por sustitución sucesiva; converge si f_fuera crece poco con la carga.
    """
    m_fresh = cfg.m_melt
    m_tower = m_fresh
    res = scr = None
    for it in range(1, max_iter + 1):
        res = simulate_tower(replace(cfg, m_melt=m_tower))
        scr = screen(res, d_min_mm, d_max_mm)
        if scr.on_spec <= 0.05:
            raise ValueError("Menos de 5 % en especificación: el reciclo no tiene estado estacionario")
        m_new = m_fresh / scr.on_spec  # producto en especificación = fundido fresco
        if progress:
            progress(it, m_new * 3.6)
        if abs(m_new - m_tower) / m_tower < tol:
            m_tower = m_new
            return RecycleResult(m_fresh, m_tower, m_tower - m_fresh, res, scr, it, True)
        m_tower = m_new
    return RecycleResult(m_fresh, m_tower, m_tower - m_fresh, res, scr, max_iter, False)


# --------------------------------------------------------------------------
# Enfriador de lecho fluidizado
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class CoolerResult:
    T_in: float
    T_out: float
    m_air: float
    Q_kW: float


def fluid_bed_cooler(m_solid: float, cp_s: float, T_in: float, m_air: float, T_air_in: float,
                     Y_water: float = 0.01) -> CoolerResult:
    """Lecho fluidizado bien mezclado; el aire sale a la temperatura del lecho.

    Supuesto P3: en lechos fluidizados de partículas de ~2 mm el gas alcanza la
    temperatura del lecho a pocos centímetros del distribuidor (Kunii y
    Levenspiel, 1991, cap. 11), de modo que el enfriador es limitado por el
    balance de energía, no por la transferencia:

        m_s cp_s (T_in - T) = m_a cp_a (T - T_a,in)
    """
    if m_solid <= 0 or m_air <= 0:
        raise ValueError("Los flujos deben ser positivos")
    cp_a = (air_table().at(0.5 * (T_in + T_air_in))[3] + Y_water * CP_VAPOR) / (1 + Y_water)
    Cs, Ca = m_solid * cp_s, m_air * cp_a
    T = (Cs * T_in + Ca * T_air_in) / (Cs + Ca)
    return CoolerResult(T_in, T, m_air, Cs * (T_in - T) / 1e3)


def cooler_air_required(m_solid: float, cp_s: float, T_in: float, T_target: float, T_air_in: float,
                        Y_water: float = 0.01) -> float:
    """Flujo de aire [kg/s] para llevar el producto a T_target (mismo supuesto P3)."""
    if not T_air_in < T_target < T_in:
        return 0.0 if T_target >= T_in else float("inf")
    cp_a = (air_table().at(0.5 * (T_target + T_air_in))[3] + Y_water * CP_VAPOR) / (1 + Y_water)
    return m_solid * cp_s * (T_in - T_target) / (cp_a * (T_target - T_air_in))


def process_table(rec: RecycleResult, cooler: CoolerResult) -> pd.DataFrame:
    """Tabla de corrientes del diagrama de bloques [t/h, °C]."""
    s = rec.screen
    rows = [
        ("1 Fundido fresco", rec.m_fresh * 3.6, rec.tower.config.T_melt_in - 273.15),
        ("2 Reciclo (redisuelto)", rec.m_recycle * 3.6, rec.tower.config.T_melt_in - 273.15),
        ("3 Alimentación a la torre", rec.m_tower * 3.6, rec.tower.config.T_melt_in - 273.15),
        ("4 Prills al tamiz", rec.m_tower * 3.6, rec.tower.T_product - 273.15),
        ("5 Fuera de especificación", rec.m_tower * s.offspec * 3.6, float("nan")),
        ("6 En especificación al enfriador", rec.m_tower * s.on_spec * 3.6, cooler.T_in - 273.15),
        ("7 Producto", rec.m_tower * s.on_spec * 3.6, cooler.T_out - 273.15),
    ]
    return pd.DataFrame(rows, columns=["corriente", "flujo_t_h", "T_C"])
