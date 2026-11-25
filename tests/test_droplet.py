"""Verificación del modelo de gota contra soluciones analíticas."""
import numpy as np
import pytest
from scipy.optimize import brentq

from prilling.correlations import terminal_velocity
from prilling.droplet import simulate_droplet, simulate_droplet_lumped
from prilling.properties import UREA, air_table

# Material "inerte" sin cambio de fase en el rango probado (siempre líquido)
INERT = UREA.with_(T_melt=200.0, k_l=0.5, k_s=0.5, cp_l=1800.0, cp_s=1800.0)


def sphere_series(Bi, Fo, n_terms=60):
    """Esfera con convección (Incropera §5.6): theta_centro y theta_media."""
    zetas = []
    for n in range(n_terms):
        lo, hi = n * np.pi + 1e-9, (n + 1) * np.pi - 1e-9
        zetas.append(brentq(lambda z: 1 - z / np.tan(z) - Bi, lo, hi))
    z = np.array(zetas)
    C = 4 * (np.sin(z) - z * np.cos(z)) / (2 * z - np.sin(2 * z))
    center = np.sum(C * np.exp(-z**2 * Fo))
    mean = np.sum(3 * C * (np.sin(z) - z * np.cos(z)) / z**3 * np.exp(-z**2 * Fo))
    return center, mean


def test_conduction_matches_analytical_series():
    d, h, Ta, T0 = 2e-3, 300.0, 300.0, 400.0
    R = d / 2
    r = simulate_droplet(d, INERT, 5.0, T0, Ta, 0.0, n_shells=40, h_const=h, rtol=1e-8)
    alpha = INERT.k_l / (INERT.rho_l * INERT.cp_l)
    Bi = h * R / INERT.k_l
    for i in (len(r.t) // 4, len(r.t) // 2, len(r.t) - 1):
        Fo = alpha * r.t[i] / R**2
        c, m = sphere_series(Bi, Fo)
        assert (r.T_center[i] - Ta) / (T0 - Ta) == pytest.approx(c, abs=0.01)
        assert (r.T_mean[i] - Ta) / (T0 - Ta) == pytest.approx(m, abs=0.01)


def test_lumped_matches_exponential():
    d, h, Ta, T0 = 2e-3, 200.0, 300.0, 400.0
    r = simulate_droplet_lumped(d, INERT, 5.0, T0, Ta, 0.0, h_const=h, rtol=1e-8)
    tau = INERT.rho_l * INERT.cp_l * d / (6 * h)
    expected = Ta + (T0 - Ta) * np.exp(-r.t / tau)
    np.testing.assert_allclose(r.T_mean, expected, atol=0.05)


def test_solidification_time_quasi_steady():
    """Solidificación de una esfera a T_m con convección, Ste -> 0:
    t_f = rho L R / (3 h dT) * (1 + Bi/2)   (London & Seban 1943)."""
    mat = UREA.with_(cp_s=5.0, cp_l=5.0)  # Stefan ~ 0
    d, h, dT = 2e-3, 250.0, 100.0
    R = d / 2
    Ta = mat.T_melt - dT
    r = simulate_droplet(d, mat, 1e4, mat.T_melt + 0.01, Ta, 0.0, n_shells=60, h_const=h,
                         rtol=1e-7, t_max=60.0, n_out=4000)
    Bi = h * R / mat.k_s
    t_exact = mat.rho_l * mat.L_fus * R / (3 * h * dT) * (1 + Bi / 2)
    assert r.t_solid == pytest.approx(t_exact, rel=0.03)


def test_high_conductivity_converges_to_lumped():
    mat = UREA.with_(k_s=500.0, k_l=500.0)
    a = simulate_droplet(1.5e-3, mat, 60.0, 413.15, 308.15, 1.0)
    b = simulate_droplet_lumped(1.5e-3, mat, 60.0, 413.15, 308.15, 1.0)
    assert a.T_mean[-1] == pytest.approx(b.T_mean[-1], abs=0.3)
    assert a.z_solid == pytest.approx(b.z_solid, rel=0.02)


def test_reaches_terminal_velocity_in_still_air():
    d = 1.5e-3
    r = simulate_droplet(d, UREA, 80.0, 413.15, 308.15, 0.0)
    # T de película al final ~ aire frío (el prill ya se enfrió)
    rho_a, mu_a, _, _ = air_table().at(0.5 * (r.T_surface[-1] + 308.15))
    assert r.v[-1] == pytest.approx(terminal_velocity(d, UREA.rho_l, rho_a, mu_a), rel=0.02)


def test_energy_conservation_with_phase_change():
    """Calor perdido por convección integrado = caída de entalpía de la gota."""
    d = 2.0e-3
    r = simulate_droplet(d, UREA, 70.0, 413.15, 308.15, 1.0, n_out=3000)
    A = np.pi * d**2
    m = UREA.rho_l * np.pi * d**3 / 6
    q = r.h_conv * A * (r.T_surface - 308.15)
    Q = np.trapezoid(q, r.t) if hasattr(np, "trapezoid") else np.trapz(q, r.t)
    assert Q == pytest.approx(m * (r.h_mean[0] - r.h_mean[-1]), rel=0.01)


def test_bigger_drops_solidify_lower():
    zs = [simulate_droplet(d, UREA, 200.0, 413.15, 308.15, 1.0).z_solid for d in (1e-3, 1.5e-3, 2e-3)]
    assert zs[0] < zs[1] < zs[2]
