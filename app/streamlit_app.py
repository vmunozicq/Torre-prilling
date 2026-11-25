"""App Streamlit: torre de prilling con PSD.

    streamlit run app/streamlit_app.py
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from prilling.droplet import biot_number, simulate_droplet  # noqa: E402
from prilling.properties import MATERIALS, air_table, thermo_comparison  # noqa: E402
from prilling.psd import RosinRammler, rayleigh_drop_diameter  # noqa: E402
from prilling.sensitivity import min_height_map  # noqa: E402
from prilling.tower import TowerConfig, simulate_tower  # noqa: E402

RES = ROOT / "data" / "resultados"
IMG = ROOT / "docs" / "img"

st.set_page_config(page_title="Torre de prilling · PSD", page_icon="🧪", layout="wide")


def calibrated_h_factor() -> float:
    try:
        return round(float(json.loads((RES / "resumen.json").read_text())["h_factor_calibrado"]), 2)
    except Exception:
        return 1.0


# ---------------------------------------------------------------- sidebar
st.sidebar.title("Parámetros")
mat_name = st.sidebar.selectbox("Material", list(MATERIALS), index=0)
material, cas = MATERIALS[mat_name]
if "no validado" in material.name:
    st.sidebar.warning("Nitrato de amonio: propiedades preliminares; no modela las transiciones sólido-sólido ni el agua del fundido.")

st.sidebar.subheader("Torre")
H = st.sidebar.slider("Altura de caída H [m]", 20.0, 120.0, 79.9, 0.1)
D = st.sidebar.slider("Diámetro de torre [m]", 5.0, 25.0, 17.0, 0.5)
m_melt_th = st.sidebar.number_input("Flujo de fundido [t/h]", 5.0, 200.0, 75.0, 1.0)
m_air_th = st.sidebar.number_input("Flujo de aire [t/h]", 50.0, 2000.0, 775.0, 5.0)
T_air_in = st.sidebar.slider("T aire de entrada [°C]", 0.0, 50.0, 35.0, 0.5)
T_melt_in = st.sidebar.slider("T fundido [°C]", float(material.T_melt - 273.15 + 1), float(material.T_melt - 273.15 + 30),
                              float(min(140.0, material.T_melt - 273.15 + 30)), 0.5)

st.sidebar.subheader("Distribución de gotas (Rosin-Rammler)")
psd_mode = st.sidebar.radio("Definir por", ["d50 de gota", "Orificio del balde (Rayleigh)"], horizontal=True)
if psd_mode == "d50 de gota":
    d50_mm = st.sidebar.slider("d50 de gota [mm]", 0.8, 3.5, 1.70, 0.05)
else:
    hole = st.sidebar.slider("Diámetro de orificio [mm]", 0.4, 1.8, 0.9, 0.05)
    d50_mm = rayleigh_drop_diameter(hole * 1e-3) * 1e3
    st.sidebar.caption(f"Rayleigh: d50 ≈ 1,89 × orificio = **{d50_mm:.2f} mm**")
n_rr = st.sidebar.slider("Índice de uniformidad n", 2.0, 12.0, 6.0, 0.5)
spec = st.sidebar.slider("Rango de especificación del prill [mm]", 0.5, 4.0, (1.0, 2.4), 0.1)

st.sidebar.subheader("Transferencia de calor")
h_factor = st.sidebar.slider("h_factor (calibración)", 0.3, 1.5, calibrated_h_factor(), 0.01,
                             help="1,00 = correlación pura. El valor por defecto es el calibrado contra la planta de Torkashvand et al. (2019).")
nusselt = st.sidebar.selectbox("Correlación de Nusselt", ["Ranz-Marshall", "Whitaker"])
radiation = st.sidebar.checkbox("Incluir radiación", False)
k_s = st.sidebar.slider("k del sólido [W/m K]", 0.2, 1.2, float(material.k_s), 0.05, help="[VERIFICAR] valor supuesto, ver docs/modelo.md")
material = material.with_(k_s=k_s)

d63 = d50_mm * 1e-3 / np.log(2.0) ** (1.0 / n_rr)
cfg = TowerConfig(H=H, D=D, m_melt=m_melt_th / 3.6, m_air=m_air_th / 3.6, T_air_in=T_air_in + 273.15,
                  T_melt_in=T_melt_in + 273.15, psd=RosinRammler(d63, n_rr), material=material,
                  nusselt=nusselt, radiation=radiation, h_factor=h_factor, n_classes=8, n_shells=12)


@st.cache_data(show_spinner=False)
def run_tower(cfg_repr: str, _cfg: TowerConfig):
    return simulate_tower(_cfg)


@st.cache_data(show_spinner=False)
def run_drop(d_mm, H, T_melt, T_air, u_air, mat_repr, _mat, nusselt, radiation, h_factor):
    return simulate_droplet(d_mm * 1e-3, _mat, H, T_melt, T_air, u_air, nusselt=nusselt,
                            radiation=radiation, h_factor=h_factor, n_shells=20)


# ---------------------------------------------------------------- header
st.title("Torre de prilling: enfriamiento, solidificación y PSD del producto")
st.caption("Proyecto personal de portafolio · datos públicos (bibliografía) · desarrollado con Claude Code. "
           "Ecuaciones, supuestos y límites en docs/modelo.md.")

tab_drop, tab_tower, tab_val, tab_sens, tab_props = st.tabs(
    ["Gota individual", "Torre completa + PSD", "Validación", "Sensibilidad", "Propiedades"])

# ---------------------------------------------------------------- gota
with tab_drop:
    c1, c2 = st.columns([1, 3])
    with c1:
        d_mm = st.slider("Diámetro de gota [mm]", 0.5, 4.0, 1.8, 0.05)
        rho_a = air_table().at(cfg.T_air_in)[0]
        u_air = cfg.m_air / (rho_a * cfg.area)
        st.caption(f"Aire a T constante = {T_air_in:.1f} °C, subiendo a {u_air:.2f} m/s (sin calentamiento: "
                   "caso conservador para la solidificación; la pestaña Torre acopla el aire).")
    r = run_drop(d_mm, H, cfg.T_melt_in, cfg.T_air_in, u_air, repr(material), material, nusselt, radiation, h_factor)
    e = r.exit
    Bi = biot_number(float(np.median(r.h_conv)), r.d, material.k_s)
    with c1:
        st.metric("Tiempo de residencia", f"{e['t_residencia_s']:.1f} s")
        st.metric("Velocidad final", f"{e['v_salida_m_s']:.2f} m/s")
        zs = e["z_solidificacion_m"]
        st.metric("Solidifica a los", f"{zs:.1f} m" if np.isfinite(zs) else "no solidifica")
        st.metric("T media a la salida", f"{e['T_media_salida_C']:.1f} °C")
        st.metric("Biot (mediana)", f"{Bi:.2f}", help="Bi > 0,1: el gradiente interno importa; el modelo concentrado no sirve.")
    with c2:
        fig = make_subplots(rows=1, cols=3, subplot_titles=("Temperatura", "Fracción líquida", "Perfil radial a la salida"))
        for y, name, dash in ((r.T_surface, "superficie", "dot"), (r.T_mean, "media", None), (r.T_center, "centro", "dash")):
            fig.add_trace(go.Scatter(x=y - 273.15, y=r.z, name=name, line=dict(dash=dash)), 1, 1)
        fig.add_trace(go.Scatter(x=r.liquid_fraction, y=r.z, name="fracción líquida", showlegend=False), 1, 2)
        fig.add_trace(go.Scatter(x=r.shells_r * 1e3, y=r.T_profile_exit - 273.15, mode="lines+markers", showlegend=False), 1, 3)
        fig.update_yaxes(autorange="reversed", title_text="caída desde la cabeza [m]", row=1, col=1)
        fig.update_yaxes(autorange="reversed", row=1, col=2)
        fig.update_xaxes(title_text="°C", row=1, col=1)
        fig.update_xaxes(title_text="radio [mm]", row=1, col=3)
        fig.update_yaxes(title_text="°C", row=1, col=3)
        fig.update_layout(height=460, margin=dict(t=40, b=10))
        st.plotly_chart(fig, width="stretch")

# ---------------------------------------------------------------- torre
with tab_tower:
    st.write("Acopla todas las clases de tamaño con el aire en contracorriente (iteración de punto fijo).")
    if st.button("Simular torre", type="primary") or "tower" in st.session_state:
        with st.spinner("Simulando 8 clases de tamaño acopladas con el aire… (10-30 s)"):
            res = run_tower(repr(cfg), cfg)
        st.session_state["tower"] = True
        s = res.summary()
        df = res.classes.copy()
        in_spec = df["d_prill_mm"].between(*spec)
        ok = in_spec & df["solidifica"]
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("T producto (mezcla)", f"{s['T_producto_C']:.1f} °C")
        m2.metric("T aire de salida", f"{s['T_aire_salida_C']:.1f} °C")
        m3.metric("Masa con líquido residual", f"{s['masa_con_liquido_residual_%']:.1f} %")
        m4.metric("Producto en especificación", f"{100 * df.loc[ok, 'fraccion_masica'].sum():.1f} %")
        if not s["convergio"]:
            st.warning("La iteración no convergió: revisa los parámetros.")

        c1, c2 = st.columns(2)
        with c1:
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=res.T_air - 273.15, y=cfg.H - res.z, name="aire", line=dict(color="black", dash="dash")))
            for rr in res.droplets:
                fig.add_trace(go.Scatter(x=rr.T_mean - 273.15, y=cfg.H - rr.z, name=f"{rr.d*1e3:.2f} mm", opacity=0.8))
            fig.update_layout(title="Temperaturas en la torre", xaxis_title="°C", yaxis_title="altura desde el fondo [m]", height=420)
            st.plotly_chart(fig, width="stretch")
        with c2:
            colors = np.where(~df["solidifica"], "#d62728", np.where(in_spec, "#2ca02c", "#ff7f0e"))
            fig = go.Figure(go.Bar(x=df["d_prill_mm"], y=100 * df["fraccion_masica"], marker_color=colors,
                                   customdata=df[["T_media_salida_C", "fraccion_liquida_salida"]],
                                   hovertemplate="d = %{x:.2f} mm<br>%{y:.1f} % masa<br>T = %{customdata[0]:.1f} °C<br>líquido = %{customdata[1]:.2f}<extra></extra>"))
            fig.add_vrect(x0=spec[0], x1=spec[1], fillcolor="green", opacity=0.07, line_width=0)
            fig.update_layout(title="PSD del producto (verde: OK, naranja: fuera de rango, rojo: no solidifica)",
                              xaxis_title="diámetro del prill [mm]", yaxis_title="% en masa", height=420)
            st.plotly_chart(fig, width="stretch")
        st.dataframe(df.round(3), width="stretch", hide_index=True)
        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as xw:
            pd.DataFrame([s]).to_excel(xw, sheet_name="resumen", index=False)
            df.to_excel(xw, sheet_name="clases", index=False)
            pd.DataFrame({"z_desde_cabeza_m": res.z, "T_aire_C": res.T_air - 273.15, "u_aire_m_s": res.u_air}).to_excel(xw, sheet_name="perfil_aire", index=False)
        st.download_button("Descargar Excel", buf.getvalue(), "torre_prilling.xlsx")

# ---------------------------------------------------------------- validación
with tab_val:
    st.subheader("Caso: torre industrial de urea (Torkashvand et al., 2019, Tabla 1)")
    st.caption("H = 79,9 m, D = 17 m, 75 t/h de urea a 140 °C, 775 t/h de aire a 35 °C. DOI 10.2298/CICEQ180827030T")
    try:
        summ = json.loads((RES / "resumen.json").read_text())
        a, b, c = st.columns(3)
        a.metric("RMSE sin calibrar", f"{summ['rmse_sin_calibrar_C']:.1f} °C")
        b.metric(f"RMSE con h_factor = {summ['h_factor_calibrado']:.2f}", f"{summ['rmse_calibrado_C']:.1f} °C")
        c.metric("RMSE del modelo de los autores", f"{summ['rmse_torkashvand_C']:.1f} °C")
        st.dataframe(pd.concat([pd.read_csv(RES / "validacion_sin_calibrar.csv"), pd.read_csv(RES / "validacion_calibrada.csv")]).round(1),
                     width="stretch", hide_index=True)
        bal = summ["balance_datos_planta"]
        st.info(f"Los datos de planta no cierran el balance de energía: la urea cede {bal['Q_cedido_urea_kW']:.0f} kW y el aire gana "
                f"{bal['Q_ganado_aire_kW']:.0f} kW ({bal['descuadre_%']:.1f} % de descuadre). Ningún modelo adiabático puede ajustar los 5 puntos a la vez.")
        st.image(str(IMG / "validacion_perfil_aire.png"))
    except FileNotFoundError:
        st.warning("Ejecuta `python scripts/run_all.py` para generar los resultados de validación.")

# ---------------------------------------------------------------- sensibilidad
with tab_sens:
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Tornado (±20 %, caso calibrado)")
        if (IMG / "tornado.png").exists():
            st.image(str(IMG / "tornado.png"))
            st.dataframe(pd.read_csv(RES / "tornado.csv").round(2), width="stretch", hide_index=True)
    with c2:
        st.subheader("Altura mínima para solidificar")
        if st.button("Calcular mapa con los parámetros actuales"):
            with st.spinner("Calculando…"):
                rho_a = air_table().at(cfg.T_air_in)[0]
                mp = min_height_map(material, np.round(np.arange(0.8, 3.61, 0.2), 2), [T_air_in - 10, T_air_in, T_air_in + 10],
                                    u_air=cfg.m_air / (rho_a * cfg.area), T_melt_in_C=T_melt_in, nusselt=nusselt, h_factor=h_factor)
            fig = go.Figure()
            for Ta, g in mp.groupby("T_aire_C"):
                fig.add_trace(go.Scatter(x=g["d_mm"], y=g["z_solidificacion_m"], name=f"aire {Ta:.0f} °C", mode="lines+markers"))
            fig.add_hline(y=H, line_dash="dot", annotation_text=f"H = {H:.1f} m")
            fig.update_layout(xaxis_title="diámetro de gota [mm]", yaxis_title="altura de caída [m]", height=420)
            st.plotly_chart(fig, width="stretch")
        elif (IMG / "altura_minima.png").exists():
            st.image(str(IMG / "altura_minima.png"))

# ---------------------------------------------------------------- propiedades
with tab_props:
    st.subheader(f"{material.name}: valores del modelo vs. `thermo`")
    st.dataframe(pd.DataFrame(thermo_comparison(material, cas)), width="stretch", hide_index=True)
    st.caption("`thermo` estima la densidad de la urea fundida ~35 % por debajo del valor medido: por eso el modelo usa valores de "
               "literatura y deja `thermo` solo como contraste. Fuentes en SOURCES.md.")
    st.write(f"Propiedades del aire: {air_table().source}.")
