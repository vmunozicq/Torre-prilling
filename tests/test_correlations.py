import pytest

from prilling.correlations import drag_coefficient, nusselt_ranz_marshall, nusselt_whitaker, terminal_velocity
from prilling.properties import UREA, air_table


def test_drag_stokes_limit():
    # Re -> 0: Cd -> 24/Re
    assert drag_coefficient(1e-3) == pytest.approx(24 / 1e-3, rel=5e-3)


@pytest.mark.parametrize("Re, cd_ref", [(1.0, 26.5), (10.0, 4.2), (100.0, 1.09), (1000.0, 0.46)])
def test_drag_standard_curve(Re, cd_ref):
    # Curva estándar de la esfera (Clift, Grace & Weber 1978, Tabla 5.1), +-6 %
    assert drag_coefficient(Re) == pytest.approx(cd_ref, rel=0.06)


def test_nusselt_conduction_limit():
    assert nusselt_ranz_marshall(0.0, 0.7) == pytest.approx(2.0)
    assert nusselt_whitaker(0.0, 0.7) == pytest.approx(2.0)


def test_terminal_velocity_stokes():
    rho_a, mu_a, _, _ = air_table().at(300.0)
    d = 10e-6  # Re ~ 0.002, régimen de Stokes
    v_stokes = (UREA.rho_l - rho_a) * 9.81 * d**2 / (18 * mu_a)
    assert terminal_velocity(d, UREA.rho_l, rho_a, mu_a) == pytest.approx(v_stokes, rel=0.02)


def test_air_properties_coolprop():
    rho, mu, k, cp = air_table().at(300.0)
    # Gas ideal para rho; Incropera, Tabla A.4 (aire, 300 K) para el resto
    assert rho == pytest.approx(101325 / (287.05 * 300.0), rel=0.005)
    assert mu == pytest.approx(184.6e-7, rel=0.02)
    assert k == pytest.approx(26.3e-3, rel=0.02)
    assert cp == pytest.approx(1007, rel=0.01)
