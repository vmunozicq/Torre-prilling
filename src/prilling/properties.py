"""Propiedades del aire y del material fundido/sólido.

Aire: CoolProp (fluido pseudo-puro "Air", 1 atm), tabulado una vez y luego
interpolado (las llamadas a CoolProp dentro del integrador serían lentas).
Si CoolProp no está instalado se usan correlaciones de gas ideal + Sutherland.

Material: valores de literatura, contrastados con `thermo`/`chemicals`.
Ver docs/modelo.md §2 y SOURCES.md: cada valor tiene su fuente y las marcas
[VERIFICAR] indican valores que Vicente debe confirmar en bibliografía.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import lru_cache

import numpy as np

P_ATM = 101_325.0  # Pa
SIGMA_SB = 5.670374419e-8  # W/(m2 K4)
CP_VAPOR = 1880.0  # J/(kg K), vapor de agua, para el cp del aire húmedo


# --------------------------------------------------------------------------
# Aire
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class AirProperties:
    """Tabla de propiedades del aire seco a 1 atm, interpolable en T [K]."""

    T: np.ndarray
    rho: np.ndarray
    mu: np.ndarray
    k: np.ndarray
    cp: np.ndarray
    source: str

    def at(self, T: float):
        """Devuelve (rho, mu, k, cp) a la temperatura T [K]."""
        return (
            float(np.interp(T, self.T, self.rho)),
            float(np.interp(T, self.T, self.mu)),
            float(np.interp(T, self.T, self.k)),
            float(np.interp(T, self.T, self.cp)),
        )


def _air_coolprop(T: np.ndarray) -> AirProperties:
    from CoolProp.CoolProp import PropsSI

    rho = np.array([PropsSI("D", "T", t, "P", P_ATM, "Air") for t in T])
    mu = np.array([PropsSI("V", "T", t, "P", P_ATM, "Air") for t in T])
    k = np.array([PropsSI("L", "T", t, "P", P_ATM, "Air") for t in T])
    cp = np.array([PropsSI("C", "T", t, "P", P_ATM, "Air") for t in T])
    return AirProperties(T, rho, mu, k, cp, "CoolProp (Lemmon et al. 2000/2004)")


def _air_fallback(T: np.ndarray) -> AirProperties:
    # Gas ideal; Sutherland para mu y k (White, Viscous Fluid Flow).
    R = 287.05
    rho = P_ATM / (R * T)
    mu = 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)
    k = 0.0241 * (T / 273.15) ** 1.5 * (273.15 + 194.0) / (T + 194.0)
    cp = 1002.5 + 275e-6 * (T - 200.0) ** 2  # ajuste simple 250-450 K
    return AirProperties(T, rho, mu, k, cp, "gas ideal + Sutherland (respaldo)")


@lru_cache(maxsize=1)
def air_table(T_min: float = 240.0, T_max: float = 480.0, n: int = 121) -> AirProperties:
    T = np.linspace(T_min, T_max, n)
    try:
        return _air_coolprop(T)
    except Exception:  # pragma: no cover - solo si falta CoolProp
        return _air_fallback(T)


# --------------------------------------------------------------------------
# Material de la gota
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Material:
    """Propiedades constantes por fase (supuesto S5 de docs/modelo.md)."""

    name: str
    T_melt: float      # K, punto de fusión / solidificación
    L_fus: float       # J/kg, calor latente de fusión
    rho_l: float       # kg/m3, densidad del fundido
    rho_s: float       # kg/m3, densidad del sólido
    cp_l: float        # J/(kg K)
    cp_s: float        # J/(kg K)
    k_l: float         # W/(m K)
    k_s: float         # W/(m K)
    emissivity: float = 0.9
    sources: dict = field(default_factory=dict, compare=False, hash=False)

    def with_(self, **kw) -> Material:
        return replace(self, **kw)

    # Método de entalpía: h(T) con referencia h = 0 en sólido a T_melt.
    def enthalpy(self, T):
        T = np.asarray(T, dtype=float)
        return np.where(
            T < self.T_melt,
            self.cp_s * (T - self.T_melt),
            self.L_fus + self.cp_l * (T - self.T_melt),
        )

    def temperature(self, h):
        """Inversa de enthalpy(): en la meseta 0 <= h <= L la T es T_melt."""
        h = np.asarray(h, dtype=float)
        T = np.where(h < 0.0, self.T_melt + h / self.cp_s, self.T_melt)
        return np.where(h > self.L_fus, self.T_melt + (h - self.L_fus) / self.cp_l, T)

    def liquid_fraction(self, h):
        return np.clip(np.asarray(h, dtype=float) / self.L_fus, 0.0, 1.0)


def _thermo_check(cas: str) -> dict:
    """Valores de `thermo` para contrastar (no se usan directamente)."""
    try:
        from thermo import Chemical

        c = Chemical(cas)
        Tm = c.Tm
        out = {"Tm": Tm, "Hfus": c.Hfus, "MW": c.MW}
        c.calculate(T=Tm - 5)
        out.update(cp_s=c.Cps, rho_s=c.rhos)
        c.calculate(T=Tm + 5)
        out.update(cp_l=c.Cpl, rho_l=c.rhol, k_l=c.kl)
        return out
    except Exception as e:  # pragma: no cover
        return {"error": str(e)}


def thermo_comparison(material: Material, cas: str) -> list[dict]:
    """Tabla literatura vs. thermo, para docs/modelo.md y la app."""
    t = _thermo_check(cas)
    if "error" in t:
        return []
    rows = [
        ("Punto de fusión [°C]", material.T_melt - 273.15, t["Tm"] - 273.15),
        ("Calor latente [kJ/kg]", material.L_fus / 1e3, t["Hfus"] / 1e3),
        ("cp sólido cerca de Tm [kJ/kg K]", material.cp_s / 1e3, (t["cp_s"] or np.nan) / 1e3),
        ("cp líquido cerca de Tm [kJ/kg K]", material.cp_l / 1e3, (t["cp_l"] or np.nan) / 1e3),
        ("Densidad fundido [kg/m3]", material.rho_l, t["rho_l"] or np.nan),
        ("Densidad sólido [kg/m3]", material.rho_s, t["rho_s"] or np.nan),
        ("k fundido [W/m K]", material.k_l, t["k_l"] or np.nan),
    ]
    return [
        {"propiedad": p, "modelo": round(a, 4), "thermo": round(b, 4),
         "dif_%": round(100 * (b - a) / a, 1) if a else np.nan}
        for p, a, b in rows
    ]


# Urea (CAS 57-13-6). Valores por defecto: ver SOURCES.md.
UREA = Material(
    name="Urea",
    T_melt=132.7 + 273.15,  # Torkashvand et al. 2019, Tabla 2 (thermo: 133.4 °C)
    L_fus=224e3,            # Torkashvand et al. 2019, Tabla 2 (thermo: 231 kJ/kg)
    rho_l=1220.0,           # Torkashvand et al. 2019, Tabla 2
    rho_s=1335.0,           # densidad cristalina de la urea [VERIFICAR fuente]
    cp_l=1920.0,            # thermo/chemicals a Tm+5 K (estimación) [VERIFICAR]
    cp_s=1620.0,            # thermo/chemicals a Tm-5 K [VERIFICAR]
    k_l=0.31,               # thermo/chemicals a Tm+5 K (estimación) [VERIFICAR]
    k_s=0.60,               # [VERIFICAR] supuesto; variable clave del análisis de sensibilidad
    emissivity=0.9,         # [VERIFICAR] supuesto
    sources={
        "T_melt": "Torkashvand et al. (2019) Tabla 2",
        "L_fus": "Torkashvand et al. (2019) Tabla 2",
        "rho_l": "Torkashvand et al. (2019) Tabla 2",
        "rho_s": "VERIFICAR",
        "cp_l": "thermo 0.6 (estimación)",
        "cp_s": "thermo 0.6",
        "k_l": "thermo 0.6 (estimación)",
        "k_s": "VERIFICAR (supuesto)",
    },
)

# Nitrato de amonio (CAS 6484-52-2). NO VALIDADO: no modela las transiciones
# sólido-sólido (IV/III/II/I) ni el agua del fundido del prill poroso.
AMMONIUM_NITRATE = Material(
    name="Nitrato de amonio (no validado)",
    T_melt=169.7 + 273.15,  # thermo/chemicals (Tm = 442.85 K)
    L_fus=73.2e3,           # thermo/chemicals (Hfus = 5.86 kJ/mol)
    rho_l=1430.0,           # [VERIFICAR] thermo da 1029 kg/m3, sospechoso
    rho_s=1725.0,           # thermo/chemicals: 1720 kg/m3
    cp_l=1775.0,            # thermo/chemicals (estimación)
    cp_s=1660.0,            # thermo/chemicals
    k_l=0.18,               # thermo/chemicals (estimación) [VERIFICAR]
    k_s=0.50,               # [VERIFICAR] supuesto
    emissivity=0.9,
    sources={"*": "thermo 0.6 salvo marcas VERIFICAR"},
)

MATERIALS = {"Urea": (UREA, "57-13-6"), "Nitrato de amonio": (AMMONIUM_NITRATE, "6484-52-2")}
