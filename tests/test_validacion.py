from dataclasses import replace

import numpy as np
import pytest

from prilling.tower import simulate_tower
from prilling.validation import compare, plant_config, plant_energy_balance, rmse


@pytest.fixture(scope="module")
def plant():
    return simulate_tower(plant_config(n_classes=8, n_shells=12))


def test_converges_and_closes_energy(plant):
    assert plant.converged
    assert plant.energy_balance_error < 1e-3


def test_air_profile_physical(plant):
    # El aire se calienta al subir: T decrece con z (z medido desde la cabeza)
    assert np.all(np.diff(plant.T_air) <= 1e-6)
    assert plant.T_air[-1] == pytest.approx(plant.config.T_air_in, abs=1e-6)
    assert plant.config.T_air_in < plant.T_product < plant.config.T_melt_in


def test_more_air_gives_colder_product(plant):
    r2 = simulate_tower(replace(plant.config, m_air=plant.config.m_air * 1.3))
    assert r2.T_product < plant.T_product
    assert r2.T_air_out < plant.T_air_out


def test_plant_data_do_not_close_energy_balance():
    # Hallazgo documentado: los datos publicados descuadran ~9-10 %
    assert 5 < plant_energy_balance()["descuadre_%"] < 15


def test_validation_regression(plant):
    """Regresión: el modelo SIN calibrar sobre la torre de Torkashvand et al.
    Si cambias el modelo y esto falla, actualiza docs/modelo.md §6 a sabiendas."""
    df = compare(plant)
    assert rmse(df) == pytest.approx(11.3, abs=1.5)
    assert df.loc[0, "modelo_C"] == pytest.approx(41.8, abs=2.0)


def test_rejects_unphysical_inputs():
    with pytest.raises(ValueError):
        plant_config(m_air=-1.0)
    with pytest.raises(ValueError):
        plant_config(T_melt_in=120 + 273.15)  # bajo el punto de fusión
