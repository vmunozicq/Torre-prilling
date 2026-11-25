"""Correlaciones de arrastre y transferencia de calor para esferas."""
from __future__ import annotations

import numpy as np


def drag_coefficient(Re: float) -> float:
    """Clift & Gauvin (1970), en Clift, Grace & Weber (1978), Tabla 5.2.

    Válida para Re < 3e5. Para Re -> 0 tiende a Stokes (24/Re).
    """
    Re = max(Re, 1e-12)
    return 24.0 / Re * (1.0 + 0.15 * Re**0.687) + 0.42 / (1.0 + 4.25e4 * Re**-1.16)


def nusselt_ranz_marshall(Re: float, Pr: float) -> float:
    """Ranz & Marshall (1952): Nu = 2 + 0.6 Re^1/2 Pr^1/3."""
    return 2.0 + 0.6 * np.sqrt(max(Re, 0.0)) * Pr ** (1.0 / 3.0)


def nusselt_whitaker(Re: float, Pr: float, mu_ratio: float = 1.0) -> float:
    """Whitaker (1972): Nu = 2 + (0.4 Re^1/2 + 0.06 Re^2/3) Pr^0.4 (mu/mu_s)^1/4.

    Rango declarado: 3.5 < Re < 7.6e4, 0.71 < Pr < 380.
    """
    Re = max(Re, 0.0)
    return 2.0 + (0.4 * Re**0.5 + 0.06 * Re ** (2.0 / 3.0)) * Pr**0.4 * mu_ratio**0.25


NUSSELT = {"Ranz-Marshall": nusselt_ranz_marshall, "Whitaker": nusselt_whitaker}


def terminal_velocity(d: float, rho_p: float, rho_a: float, mu_a: float, g: float = 9.81) -> float:
    """Velocidad terminal en aire quieto: resuelve el balance peso-empuje-arrastre."""
    from scipy.optimize import brentq

    def f(v):
        Re = rho_a * v * d / mu_a
        return (rho_p - rho_a) * g * np.pi * d**3 / 6 - 0.5 * rho_a * v**2 * np.pi * d**2 / 4 * drag_coefficient(Re)

    return brentq(f, 1e-6, 100.0)
