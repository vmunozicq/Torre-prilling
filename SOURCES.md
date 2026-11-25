# Fuentes

Toda fila de `data/raw/` tiene una columna `source_id` que apunta a esta tabla. No se redistribuyen PDFs ni figuras: solo valores extraídos.

| source_id | Fuente | DOI o URL | Licencia | Consulta | Uso en el repo |
| --- | --- | --- | --- | --- | --- |
| S01 | Torkashvand, B., Gilassi, S., Moghadam, R. (2019). Analysis of heat and moisture transfer during drying of urea particles. *Chem. Ind. Chem. Eng. Q.* 25(2), 171-182 | https://doi.org/10.2298/CICEQ180827030T | Acceso abierto (verificar licencia exacta) | 2026-11-25 | Caso de validación (Tabla 1) y propiedades de la urea (Tabla 2): `data/raw/torkashvand2019_*.csv` |
| S02 | Clift, R., Grace, J. R., Weber, M. E. (1978). *Bubbles, Drops, and Particles*. Academic Press, cap. 5 | ISBN 0-12-176950-X | Libro | 2026-11-25 | Correlación de arrastre de Clift-Gauvin y curva estándar de la esfera (tests) |
| S03 | Ranz, W. E., Marshall, W. R. (1952). Evaporation from drops. *Chem. Eng. Prog.* 48(3), 141-146 y 48(4), 173-180 | — | Artículo | 2026-11-25 | Correlación de Nusselt por defecto |
| S04 | Whitaker, S. (1972). Forced convection heat transfer correlations for flow in pipes, past flat plates, single cylinders, single spheres… *AIChE J.* 18(2), 361-371 | https://doi.org/10.1002/aic.690180219 | Artículo | 2026-11-25 | Correlación de Nusselt alternativa |
| S05 | Bell, I. H., Wronski, J., Quoilin, S., Lemort, V. (2014). Pure and pseudo-pure fluid thermophysical property evaluation and the open-source thermophysical property library CoolProp. *Ind. Eng. Chem. Res.* 53(6), 2498-2508 | https://doi.org/10.1021/ie4033999 | Software MIT | 2026-11-25 | Propiedades del aire (ecuación de Lemmon et al. 2000 y transporte de Lemmon & Jacobsen 2004) |
| S06 | Bell, C. et al. `thermo` 0.6.1 y `chemicals` (Python) | https://github.com/CalebBell/thermo | MIT | 2026-11-25 | cp y k de la urea; contraste de Tm, ΔHfus y densidades |
| S07 | Lord Rayleigh (1878). On the instability of jets. *Proc. London Math. Soc.* 10, 4-13 | https://doi.org/10.1112/plms/s1-10.1.4 | Dominio público | 2026-11-25 | d_gota = 1,89 d_chorro |
| S08 | London, A. L., Seban, R. A. (1943). Rate of ice formation. *Trans. ASME* 65, 771-778 | — | Artículo | 2026-11-25 | Solución cuasi-estacionaria de solidificación de una esfera (test) |
| S09 | Incropera, F. P., DeWitt, D. P., Bergman, T. L., Lavine, A. S. *Fundamentals of Heat and Mass Transfer*, 6.ª ed., §5.6 y Tabla A.4 | — | Libro | 2026-11-25 | Serie analítica de la esfera con convección y propiedades del aire (tests) |
| S10 | Rosin, P., Rammler, E. (1933). The laws governing the fineness of powdered coal. *J. Inst. Fuel* 7, 29-36 | — | Artículo | 2026-11-25 | Forma de la PSD |

Las referencias S02-S04 y S07-S10 son clásicas; Vicente debe confirmar páginas y ecuaciones contra su copia antes de publicar (marcado `[VERIFICAR]` en `docs/modelo.md`).
