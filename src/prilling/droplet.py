"""Trayectoria y estado térmico de una gota/prill en caída (marco lagrangiano).

Ecuaciones (ver docs/modelo.md §3):
  dz/dt = v                                   (z hacia abajo desde la cabeza de la torre)
  dv/dt = g (1 - rho_a/rho_p) - (3/4) Cd rho_a w|w| / (rho_p d),   w = v + u_aire
  Conducción radial transitoria con cambio de fase (método de entalpía,
  volúmenes finitos en N capas esféricas) y convección (+ radiación opcional)
  en la superficie.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .correlations import NUSSELT, drag_coefficient
from .properties import SIGMA_SB, Material, air_table

G = 9.81


@dataclass
class DropletResult:
    d: float                 # diámetro de la gota [m]
    t: np.ndarray            # tiempo [s]
    z: np.ndarray            # distancia caída desde la cabeza [m]
    v: np.ndarray            # velocidad hacia abajo [m/s]
    T_surface: np.ndarray    # K
    T_center: np.ndarray     # K
    T_mean: np.ndarray       # K, temperatura equivalente de la entalpía media
    h_mean: np.ndarray       # J/kg, entalpía específica media (ref.: sólido a T_melt)
    liquid_fraction: np.ndarray  # fracción másica líquida media
    h_conv: np.ndarray       # W/(m2 K), coeficiente convectivo
    Re: np.ndarray
    z_solid: float           # altura caída a la que la gota queda 100 % sólida (nan si no)
    t_solid: float
    reached_bottom: bool
    shells_r: np.ndarray     # radios de los nodos [m]
    T_profile_exit: np.ndarray  # perfil radial de T a la salida [K]

    @property
    def exit(self) -> dict:
        return {
            "d_mm": self.d * 1e3,
            "t_residencia_s": float(self.t[-1]),
            "v_salida_m_s": float(self.v[-1]),
            "T_media_salida_C": float(self.T_mean[-1] - 273.15),
            "T_superficie_salida_C": float(self.T_surface[-1] - 273.15),
            "T_centro_salida_C": float(self.T_center[-1] - 273.15),
            "fraccion_liquida_salida": float(self.liquid_fraction[-1]),
            "z_solidificacion_m": float(self.z_solid),
        }


def _const(value: float) -> Callable[[float], float]:
    return lambda z: value


def simulate_droplet(
    d: float,
    material: Material,
    H: float,
    T_drop0: float,
    T_air: Callable[[float], float] | float,
    u_air: Callable[[float], float] | float = 0.0,
    *,
    n_shells: int = 20,
    nusselt: str = "Ranz-Marshall",
    radiation: bool = False,
    v0: float = 0.0,
    internal_resistance: bool = True,
    h_factor: float = 1.0,
    h_const: float | None = None,
    n_out: int = 400,
    rtol: float = 1e-6,
    t_max: float = 600.0,
) -> DropletResult:
    """Integra una gota de diámetro d desde z=0 hasta z=H.

    h_factor multiplica el coeficiente convectivo (1.0 = correlación pura); es
    el único parámetro de calibración del modelo (docs/modelo.md §5). h_const fija
    un coeficiente constante y solo se usa en los tests de verificación.
    T_air(z) [K] y u_air(z) [m/s, positivo hacia arriba] pueden ser constantes
    o funciones de la altura caída z. Con internal_resistance=False se usa el
    modelo concentrado (Bi -> 0) en una sola celda.
    """
    if d <= 0 or H <= 0:
        raise ValueError("El diámetro y la altura deben ser positivos")
    if n_shells < 1:
        raise ValueError("n_shells debe ser >= 1")
    T_air_f = T_air if callable(T_air) else _const(float(T_air))
    u_air_f = u_air if callable(u_air) else _const(float(u_air))
    air = air_table()
    nuss = NUSSELT[nusselt]
    m = material

    N = n_shells if internal_resistance else 1
    R = d / 2.0
    r_edges = np.linspace(0.0, R, N + 1)
    r_nodes = 0.5 * (r_edges[:-1] + r_edges[1:])
    rho_p = m.rho_l  # supuesto S6: geometría y masa constantes (densidad del fundido)
    vol = 4.0 / 3.0 * np.pi * (r_edges[1:] ** 3 - r_edges[:-1] ** 3)
    mass = rho_p * vol
    A_s = 4.0 * np.pi * R**2
    m_total = mass.sum()

    def sph_res(a, b, k):
        # Resistencia de conducción de una cáscara esférica entre radios a < b.
        return (1.0 / a - 1.0 / b) / (4.0 * np.pi * k)

    def conductivity(h):
        f = m.liquid_fraction(h)
        return m.k_s + f * (m.k_l - m.k_s)

    def film(Ts, z, v):
        Ta = T_air_f(z)
        rho_a, mu_a, k_a, cp_a = air.at(0.5 * (Ts + Ta))
        w = v + u_air_f(z)
        Re = rho_a * abs(w) * d / mu_a
        Pr = cp_a * mu_a / k_a
        h = h_const if h_const is not None else h_factor * nuss(Re, Pr) * k_a / d
        return Ta, rho_a, mu_a, w, Re, h

    def rhs(t, y):
        z, v = y[0], y[1]
        hs = y[2:]
        T = m.temperature(hs)
        Ts_est = T[-1]
        Ta, rho_a, mu_a, w, Re, hc = film(Ts_est, z, v)
        # Movimiento
        Cd = drag_coefficient(Re)
        dv = G * (1.0 - rho_a / rho_p) - 0.75 * Cd * rho_a * w * abs(w) / (rho_p * d)
        # Superficie: convección + radiación linealizada hacia el entorno a T_aire
        h_tot = hc
        if radiation:
            h_tot += m.emissivity * SIGMA_SB * (Ts_est**2 + Ta**2) * (Ts_est + Ta)
        dh = np.empty(N)
        if internal_resistance:
            k = conductivity(hs)
            # Flujos entre nodos i -> i+1 a través de r_edges[i+1]
            Rin = sph_res(r_nodes[:-1], r_edges[1:-1], k[:-1]) + sph_res(r_edges[1:-1], r_nodes[1:], k[1:])
            Q = (T[:-1] - T[1:]) / Rin
            R_out = sph_res(r_nodes[-1], R, k[-1]) + 1.0 / (h_tot * A_s)
            Q_out = (T[-1] - Ta) / R_out
            dh[:] = 0.0
            dh[:-1] -= Q
            dh[1:] += Q
            dh[-1] -= Q_out
            dh /= mass
        else:
            dh[0] = -h_tot * A_s * (T[0] - Ta) / m_total
        return np.concatenate(([v, dv], dh))

    def hit_bottom(t, y):
        return y[0] - H

    hit_bottom.terminal = True
    hit_bottom.direction = 1

    y0 = np.concatenate(([0.0, v0], np.full(N, float(m.enthalpy(T_drop0)))))
    sol = solve_ivp(rhs, (0.0, t_max), y0, method="LSODA", events=hit_bottom,
                    rtol=rtol, atol=np.concatenate(([1e-6, 1e-6], np.full(N, 1e-2))),
                    dense_output=True)
    t_end = sol.t[-1]
    t = np.linspace(0.0, t_end, n_out)
    Y = sol.sol(t)
    z, v, Hs = Y[0], Y[1], Y[2:]
    Ts_nodes = m.temperature(Hs)
    h_mean = (mass[:, None] * Hs).sum(axis=0) / m_total
    f_mean = (mass[:, None] * m.liquid_fraction(Hs)).sum(axis=0) / m_total
    T_mean = m.temperature(h_mean)

    # Temperatura de superficie y coeficiente convectivo a lo largo de la trayectoria
    T_surf = np.empty_like(t)
    hconv = np.empty_like(t)
    Re_arr = np.empty_like(t)
    for i in range(len(t)):
        Ta, _, _, _, Re, hc = film(Ts_nodes[-1, i], z[i], v[i])
        if internal_resistance:
            k_out = conductivity(Hs[-1, i])
            R_cond = sph_res(r_nodes[-1], R, k_out)
            Tn = Ts_nodes[-1, i]
            h_tot = hc + (m.emissivity * SIGMA_SB * (Tn**2 + Ta**2) * (Tn + Ta) if radiation else 0.0)
            q = (Ts_nodes[-1, i] - Ta) / (R_cond + 1.0 / (h_tot * A_s))
            T_surf[i] = Ts_nodes[-1, i] - q * R_cond
        else:
            T_surf[i] = Ts_nodes[0, i]
        hconv[i], Re_arr[i] = hc, Re

    # Solidificación completa: todas las capas con fracción líquida ~0
    f_max = m.liquid_fraction(Hs).max(axis=0)
    idx = np.where(f_max <= 1e-6)[0]
    if len(idx):
        i0 = idx[0]
        z_solid, t_solid = float(z[i0]), float(t[i0])
    else:
        z_solid = t_solid = float("nan")

    return DropletResult(
        d=d, t=t, z=z, v=v, T_surface=T_surf, T_center=Ts_nodes[0], T_mean=T_mean,
        h_mean=h_mean, liquid_fraction=f_mean, h_conv=hconv, Re=Re_arr,
        z_solid=z_solid, t_solid=t_solid, reached_bottom=bool(sol.status == 1),
        shells_r=r_nodes, T_profile_exit=Ts_nodes[:, -1],
    )


def simulate_droplet_lumped(d, material, H, T_drop0, T_air, u_air=0.0, **kw) -> DropletResult:
    """Modelo concentrado de tres etapas (sensible-latente-sensible), Bi -> 0."""
    kw.pop("n_shells", None)
    return simulate_droplet(d, material, H, T_drop0, T_air, u_air, internal_resistance=False, **kw)


def biot_number(h: float, d: float, k_s: float) -> float:
    """Bi = h (R/3) / k_s, con longitud característica V/A = R/3 de la esfera."""
    return h * (d / 6.0) / k_s
