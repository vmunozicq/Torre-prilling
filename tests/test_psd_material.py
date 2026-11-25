import numpy as np
import pytest

from prilling.properties import UREA
from prilling.psd import RAYLEIGH_FACTOR, RosinRammler, rayleigh_drop_diameter, solid_diameter


def test_rosin_rammler():
    rr = RosinRammler(2e-3, 5.0)
    assert rr.cdf(2e-3) == pytest.approx(1 - np.exp(-1))
    assert rr.cdf(rr.quantile(0.37)) == pytest.approx(0.37)
    d, w, edges = rr.classes(12)
    assert w.sum() == pytest.approx(1.0)
    assert np.all(np.diff(d) > 0) and len(edges) == 13


def test_rayleigh():
    assert RAYLEIGH_FACTOR == pytest.approx(1.89, abs=0.005)
    assert rayleigh_drop_diameter(1e-3) == pytest.approx(1.89e-3, rel=0.003)


def test_solid_diameter_mass_conservation():
    ds = solid_diameter(2e-3, UREA.rho_l, UREA.rho_s)
    assert UREA.rho_s * ds**3 == pytest.approx(UREA.rho_l * (2e-3) ** 3)


@pytest.mark.parametrize("T", [300.0, UREA.T_melt - 1, UREA.T_melt + 1, 420.0])
def test_enthalpy_roundtrip(T):
    assert UREA.temperature(UREA.enthalpy(T)) == pytest.approx(T)


def test_liquid_fraction_plateau():
    assert UREA.liquid_fraction(UREA.L_fus / 2) == pytest.approx(0.5)
    assert UREA.temperature(UREA.L_fus / 2) == pytest.approx(UREA.T_melt)


def test_min_height_grows_with_diameter_and_air_temperature():
    from prilling.sensitivity import min_height_map

    mp = min_height_map(UREA, [1.0, 2.0], [25.0, 45.0], u_air=0.9, n_shells=8)
    z = mp.set_index(["T_aire_C", "d_mm"])["z_solidificacion_m"]
    assert z[(25.0, 1.0)] < z[(25.0, 2.0)]
    assert z[(25.0, 2.0)] < z[(45.0, 2.0)]
