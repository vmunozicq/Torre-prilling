"""Distribución de tamaño de gota/prill (PSD)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

RAYLEIGH_FACTOR = (1.5 * 4.508) ** (1.0 / 3.0)  # = 1.89


def rayleigh_drop_diameter(d_jet: float) -> float:
    """Ruptura de Rayleigh de un chorro laminar: d_gota = 1.89 d_chorro.

    Longitud de onda de máximo crecimiento lambda = 4.508 d_chorro (Rayleigh 1878);
    una gota por longitud de onda: pi d^3/6 = pi d_j^2/4 * lambda.
    """
    return RAYLEIGH_FACTOR * d_jet


@dataclass(frozen=True)
class RosinRammler:
    """Distribución acumulada en masa: F(d) = 1 - exp[-(d/d63)^n]."""

    d63: float  # m, diámetro característico (63.2 % en masa bajo tamaño)
    n: float    # índice de uniformidad (mayor = más estrecha)

    def cdf(self, d):
        d = np.asarray(d, dtype=float)
        return 1.0 - np.exp(-(np.maximum(d, 0.0) / self.d63) ** self.n)

    def quantile(self, F):
        F = np.asarray(F, dtype=float)
        return self.d63 * (-np.log(1.0 - F)) ** (1.0 / self.n)

    @property
    def d50(self) -> float:
        return float(self.quantile(0.5))

    def classes(self, n_classes: int = 12, F_lo: float = 0.005, F_hi: float = 0.995):
        """Discretiza en clases de igual ancho en d entre los cuantiles F_lo y F_hi.

        Devuelve (diámetro representativo, fracción másica, bordes). El
        representativo es la media geométrica de los bordes; las fracciones se
        renormalizan para sumar 1.
        """
        edges = np.linspace(self.quantile(F_lo), self.quantile(F_hi), n_classes + 1)
        w = np.diff(self.cdf(edges))
        w = w / w.sum()
        d_rep = np.sqrt(edges[:-1] * edges[1:])
        return d_rep, w, edges


def solid_diameter(d_melt, rho_l: float, rho_s: float):
    """Contracción por cambio de densidad al solidificar (masa constante, sin poros)."""
    return np.asarray(d_melt) * (rho_l / rho_s) ** (1.0 / 3.0)
