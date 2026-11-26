import numpy as np
import pytest

from prilling.process import cooler_air_required, fluid_bed_cooler, screen
from prilling.tower import simulate_tower
from prilling.validation import plant_config


def test_cooler_energy_balance():
    r = fluid_bed_cooler(20.0, 1777.0, 340.0, 30.0, 300.0)
    assert 300.0 < r.T_out < 340.0
    cp_a = (r.Q_kW * 1e3) / (30.0 * (r.T_out - 300.0))
    assert 1000 < cp_a < 1040  # el calor que pierde el sólido lo gana el aire


def test_cooler_air_required_inverts_cooler():
    m_air = cooler_air_required(20.0, 1777.0, 340.0, 323.15, 300.0)
    assert fluid_bed_cooler(20.0, 1777.0, 340.0, m_air, 300.0).T_out == pytest.approx(323.15, abs=0.05)


@pytest.fixture(scope="module")
def tower_cal():
    return simulate_tower(plant_config(h_factor=0.61, n_classes=8, n_shells=12))


def test_screen_closes_mass(tower_cal):
    s = screen(tower_cal, 1.0, 2.4)
    assert s.on_spec + s.offspec == pytest.approx(1.0, abs=1e-9)
    assert 0 < s.soft < 0.3


def test_critical_diameter_consistent(tower_cal):
    df = tower_cal.classes
    d_c = tower_cal.d_crit * 1e3
    assert df.loc[df["solidifica"], "d_gota_mm"].max() <= d_c + 1e-9
    assert np.all(df.loc[~df["solidifica"], "d_gota_mm"] >= d_c - 0.2)


def test_rejects_bad_screen(tower_cal):
    with pytest.raises(ValueError):
        screen(tower_cal, 2.4, 1.0)
