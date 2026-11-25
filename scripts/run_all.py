"""Reproduce todos los resultados de docs/modelo.md y del README.

    python scripts/run_all.py          (≈ 5-8 min en un laptop)

Escribe CSV en data/resultados/ y figuras en docs/img/.
"""
from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from prilling.properties import UREA, thermo_comparison  # noqa: E402
from prilling.sensitivity import design_curve, min_height_map, tornado  # noqa: E402
from prilling.tower import simulate_tower  # noqa: E402
from prilling.validation import calibrate_h_factor, compare, plant_config, plant_energy_balance, rmse  # noqa: E402

OUT = ROOT / "data" / "resultados"
IMG = ROOT / "docs" / "img"
OUT.mkdir(parents=True, exist_ok=True)
IMG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "axes.grid": True, "grid.alpha": 0.3})

summary = {}

# 1. Propiedades: literatura vs thermo
tc = thermo_comparison(UREA, "57-13-6")
import pandas as pd  # noqa: E402

pd.DataFrame(tc).to_csv(OUT / "propiedades_vs_thermo.csv", index=False)
print("1/5 propiedades")

# 2. Validación sin calibrar y calibrada
summary["balance_datos_planta"] = plant_energy_balance()
base = simulate_tower(plant_config())
c0 = compare(base)
f_cal, c1 = calibrate_h_factor()
cal = simulate_tower(plant_config(h_factor=f_cal))
c0.assign(caso="sin calibrar").to_csv(OUT / "validacion_sin_calibrar.csv", index=False)
c1.assign(caso=f"h_factor={f_cal:.2f}").to_csv(OUT / "validacion_calibrada.csv", index=False)
summary.update({
    "rmse_sin_calibrar_C": rmse(c0), "rmse_calibrado_C": rmse(c1),
    "rmse_torkashvand_C": rmse(c0, "error_torkashvand_C"), "h_factor_calibrado": f_cal,
    "sin_calibrar": base.summary(), "calibrado": cal.summary(),
})
cal.classes.to_csv(OUT / "clases_caso_calibrado.csv", index=False)
print("2/5 validación, h_factor =", round(f_cal, 3))

fig, ax = plt.subplots(figsize=(6.5, 4.2))
hb = cal.config.H - cal.z
ax.plot(base.T_air - 273.15, base.config.H - base.z, label="Modelo sin calibrar")
ax.plot(cal.T_air - 273.15, hb, label=f"Modelo calibrado (h_factor = {f_cal:.2f})")
pts = c0.iloc[1:]
ax.scatter(pts["planta_C"], pts["altura_desde_fondo_m"], c="k", zorder=5, label="Planta (Torkashvand et al. 2019)")
ax.scatter(pts["torkashvand_C"], pts["altura_desde_fondo_m"], marker="x", c="tab:red", zorder=5, label="Modelo de Torkashvand et al.")
ax.set_xlabel("Temperatura del aire [°C]"); ax.set_ylabel("Altura desde el fondo [m]")
ax.set_title("Perfil de aire en la torre: modelo vs planta"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(IMG / "validacion_perfil_aire.png"); plt.close(fig)

# 3. Trayectorias por clase (caso calibrado)
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
for r in cal.droplets[::2]:
    axs[0].plot(r.T_mean - 273.15, r.z, label=f"{r.d*1e3:.2f} mm")
    axs[1].plot(r.liquid_fraction, r.z)
axs[0].plot(cal.T_air - 273.15, cal.z, "k--", label="aire")
axs[0].invert_yaxis(); axs[0].set_ylabel("Distancia de caída desde la cabeza [m]")
axs[0].set_xlabel("Temperatura media de la gota [°C]"); axs[0].legend(fontsize=8)
axs[1].set_xlabel("Fracción líquida"); axs[0].set_title("Enfriamiento"); axs[1].set_title("Solidificación")
fig.tight_layout(); fig.savefig(IMG / "trayectorias.png"); plt.close(fig)
print("3/5 trayectorias")

# 4. Mapa de altura mínima y curva de diseño
diam = np.round(np.arange(0.8, 3.21, 0.2), 2)
mp = min_height_map(UREA, diam, [25, 35, 45], u_air=float(cal.u_air.mean()), h_factor=f_cal)
mp.to_csv(OUT / "altura_minima_solidificacion.csv", index=False)
fig, ax = plt.subplots(figsize=(6.5, 4.2))
for Ta, g in mp.groupby("T_aire_C"):
    ax.plot(g["d_mm"], g["z_solidificacion_m"], marker="o", ms=3, label=f"aire a {Ta:.0f} °C")
ax.axhline(79.9, color="k", ls=":", label="H torre caso (79,9 m)")
ax.set_xlabel("Diámetro de gota [mm]"); ax.set_ylabel("Altura de caída para solidificar [m]")
ax.set_title(f"Altura mínima (gota aislada, h_factor = {f_cal:.2f})"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(IMG / "altura_minima.png"); plt.close(fig)
summary["d_max_solidifica_79.9m_35C_mm"] = design_curve(UREA, 79.9, 35.0, u_air=float(cal.u_air.mean()), h_factor=f_cal)
print("4/5 mapa")

# 5. Tornado (caso calibrado, ±20 %)
tor = tornado(replace(plant_config(h_factor=f_cal), n_classes=8, n_shells=12), 0.2,
              progress=lambda i, n: print(f"   tornado {i}/{n}"))
tor.to_csv(OUT / "tornado.csv", index=False)
base_T = tor["T_producto_base_C"].iloc[0]
fig, ax = plt.subplots(figsize=(7, 4.2))
y = np.arange(len(tor))[::-1]
ax.barh(y, tor["T_producto_bajo_C"] - base_T, color="tab:blue", label="parámetro −20 %")
ax.barh(y, tor["T_producto_alto_C"] - base_T, color="tab:orange", label="parámetro +20 %")
ax.set_yticks(y, tor["parametro"]); ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel(f"Cambio en T del producto [°C] (base {base_T:.1f} °C)")
ax.set_title("Sensibilidad de la temperatura del producto"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(IMG / "tornado.png"); plt.close(fig)
print("5/5 tornado")

(OUT / "resumen.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float))
print(json.dumps(summary, indent=2, ensure_ascii=False, default=float))
