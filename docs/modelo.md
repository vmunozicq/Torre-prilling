# Modelo: torre de prilling con PSD

> **Para Vicente.** Este documento lo redactó Claude Code a partir de la bibliografía, para que lo revises, lo corrijas y lo hagas tuyo antes de publicar. Todo lo marcado **[VERIFICAR]** es un valor o una cita que tienes que confirmar contra tu fuente; todo lo marcado **[DECIDIR]** es un supuesto de modelado que tienes que poder defender en una entrevista. Al final hay preguntas de entrevista que este modelo te va a provocar.

## 1. Pregunta de ingeniería

¿Qué diámetro de gota alcanza a solidificar completamente en una torre de altura H, con un flujo y una temperatura de aire dados, y a qué temperatura y con qué distribución de tamaños sale el producto?

Es la decisión que se toma al dimensionar una torre nueva (altura, flujo de aire) o al operar una existente (orificios del balde, carga de fundido, aire en verano). Una gota que llega al fondo con el núcleo líquido se deforma, se aglomera en la correa y genera finos y terrones.

## 2. Entradas y salidas

| Variable | Tipo | Unidad | Rango probado | Caso base |
| --- | --- | --- | --- | --- |
| H, altura de caída | entrada | m | 20-120 | 79,9 |
| D, diámetro de la torre | entrada | m | 5-25 | 17 |
| ṁ_f, flujo de fundido | entrada | kg/s | 1-55 | 20,8 (75 t/h) |
| ṁ_a, flujo de aire | entrada | kg/s | 15-550 | 215 (775 t/h) |
| T_a,in, aire de entrada | entrada | °C | 0-50 | 35 |
| T_f,in, fundido | entrada | °C | Tm+1 a Tm+30 | 140 |
| d63, n (Rosin-Rammler) | entrada | mm, – | 0,8-3,5; 2-12 | d50 = 1,70; n = 6 |
| h_factor | calibración | – | 0,3-1,5 | 1,00 / 0,61 |
| T del producto (mezcla) | salida | °C | | |
| T del aire de salida, perfil T_a(z) | salida | °C | | |
| Masa con líquido residual al fondo | salida | % | | |
| PSD del prill y fracción en especificación | salida | % masa | | |

## 3. Ecuaciones

### 3.1 Movimiento de la gota (marco lagrangiano, 1D vertical)

z se mide hacia abajo desde la cabeza; v > 0 hacia abajo; el aire sube con u_a = ṁ_a / (ρ_a A).

$$\frac{dz}{dt} = v \qquad \frac{dv}{dt} = g\left(1-\frac{\rho_a}{\rho_p}\right) - \frac{3}{4}\,\frac{C_D\,\rho_a\,w|w|}{\rho_p\,d}, \qquad w = v + u_a$$

Arrastre de Clift-Gauvin (S02), válido hasta Re ≈ 3·10⁵:

$$C_D = \frac{24}{Re}\left(1+0{,}15\,Re^{0{,}687}\right) + \frac{0{,}42}{1+4{,}25\cdot10^{4}\,Re^{-1{,}16}}, \qquad Re = \frac{\rho_a |w| d}{\mu_a}$$

### 3.2 Transferencia de calor en la superficie

$$Nu = 2 + 0{,}6\,Re^{1/2}Pr^{1/3} \;\;\text{(Ranz-Marshall, S03)} \qquad h = h_\text{factor}\cdot\frac{Nu\,k_a}{d}$$

Alternativa: Whitaker (S04). Propiedades del aire a la temperatura de película (T_s + T_a)/2 con CoolProp (S05). Radiación opcional, linealizada hacia un entorno a T_a: h_rad = εσ(T_s² + T_a²)(T_s + T_a). Por defecto va apagada porque así lo trata la referencia de validación (aporta ≈ 4 % del calor; ver tornado).

### 3.3 Conducción con cambio de fase dentro de la gota (método de entalpía)

La gota se divide en N capas esféricas de igual espesor (N = 12-20). Cada capa guarda su entalpía específica h_i (referencia: sólido a T_m):

$$h(T) = \begin{cases} c_{p,s}(T-T_m) & T<T_m \\ L + c_{p,l}(T-T_m) & T>T_m \end{cases} \qquad f_l = \mathrm{clip}(h/L,\,0,\,1)$$

En la meseta 0 ≤ h ≤ L la capa está a T_m y su fracción líquida es h/L. El balance de cada capa es

$$m_i\frac{dh_i}{dt} = Q_{i-1\to i} - Q_{i\to i+1}, \qquad Q_{i\to i+1} = \frac{T_i - T_{i+1}}{R_{i,i+1}}, \qquad R = \frac{1}{4\pi k}\left(\frac{1}{r_a}-\frac{1}{r_b}\right)$$

con k = k_s + f_l (k_l − k_s). La última capa entrega calor al aire a través de media capa de conducción en serie con 1/(hA). El sistema [z, v, h_1…h_N] se integra con LSODA (scipy).

**Por qué no un modelo concentrado de tres etapas.** El número de Biot del prill de urea es Bi = h(R/3)/k_s ≈ 0,1-0,2 con k_s = 0,6 W/m K, y Bi = hR/k_s ≈ 0,3-0,6: el núcleo queda bastante más caliente que la superficie. El modelo concentrado (incluido como `simulate_droplet_lumped` para comparar) solidifica la gota de 2,4 mm unos 20 m más arriba que el modelo radial. Esa diferencia es una buena pregunta de entrevista.

### 3.4 Torre en contracorriente (acoplamiento con el aire)

Volumen de control entre la altura z y el fondo, pared adiabática, aire que entra por abajo:

$$\dot m_a\,[h_a(T_a(z)) - h_a(T_{a,in})] = \dot m_f \sum_j w_j\,[\bar H_j(z) - \bar H_j(H)]$$

donde w_j es la fracción másica de la clase j y H̄_j(z) su entalpía media en la altura z. Se itera por punto fijo (relajación 0,7) sobre el perfil T_a(z) hasta que cambie menos de 0,05 K; converge en 5 a 8 iteraciones. La entalpía del aire incluye la humedad solo en el cp: c_p = (c_p,seco + Y c_p,vapor)/(1 + Y).

### 3.5 PSD

Distribución acumulada en masa de Rosin-Rammler (S10): F(d) = 1 − exp[−(d/d63)ⁿ], discretizada en 8-12 clases entre los cuantiles 0,5 % y 99,5 %. Opcionalmente, d50 a partir del orificio del balde por ruptura de Rayleigh (S07): d = 1,89 d_orificio. El diámetro del prill sólido es d_s = d (ρ_l/ρ_s)^{1/3} (masa constante, sin porosidad).

## 4. Propiedades

| Propiedad | Urea, valor usado | Fuente | `thermo` 0.6 | Comentario |
| --- | --- | --- | --- | --- |
| T_m | 132,7 °C | S01 Tabla 2 | 133,4 °C | concuerdan |
| L | 224 kJ/kg | S01 Tabla 2 | 231 kJ/kg | concuerdan (3 %) |
| ρ_l | 1220 kg/m³ | S01 Tabla 2 | 799 kg/m³ | `thermo` subestima ~35 %: **no usar** |
| ρ_s | 1335 kg/m³ | **[VERIFICAR]** densidad cristalina | 900 kg/m³ | ídem |
| c_p,l | 1,92 kJ/kg K | `thermo` (estimación) **[VERIFICAR]** | | S01 usa 1,334 para ambas fases |
| c_p,s | 1,62 kJ/kg K | `thermo` **[VERIFICAR]** | | |
| k_l | 0,31 W/m K | `thermo` (estimación) **[VERIFICAR]** | | |
| k_s | 0,60 W/m K | **[VERIFICAR] supuesto** | sin dato | S01 imprime un valor con unidades inconsistentes; es el parámetro de material más incierto |
| ε | 0,9 | **[VERIFICAR] supuesto** | | solo si se activa radiación |

Aire: CoolProp, aire seco pseudo-puro a 1 atm, tabulado 240-480 K e interpolado. Hallazgo útil para la entrevista: `thermo` da la densidad de la urea fundida y del nitrato de amonio fundido con errores de 30-35 %, por eso el modelo toma esos valores de literatura y deja `thermo` solo como contraste (pestaña Propiedades de la app).

El nitrato de amonio está disponible en la app como **no validado**: no modela las transiciones sólido-sólido (IV↔III a 32 °C, II↔I a 125 °C), que en un prill de nitrato son justamente lo que produce la degradación por ciclos térmicos, ni el agua del fundido de un prill poroso.

## 5. Supuestos

| # | Supuesto | Por qué es razonable | Cuándo deja de serlo |
| --- | --- | --- | --- |
| S1 | Estado estacionario de la torre | operación continua | arranques, cambios de carga |
| S2 | Flujo de aire uniforme en la sección (tapón, 1D) | estándar en modelos 1D de prilling | torres reales tienen canalización y zonas muertas: **[DECIDIR]** es la principal sospecha detrás del h_factor < 1 |
| S3 | Gotas esféricas y sin interacción (sin coalescencia ni estela) | gotas pequeñas, tensión superficial alta | cerca del balde la nube es densa |
| S4 | Solo movimiento vertical; la gota parte con v₀ = 0 | el balde lanza la gota en horizontal, que se frena en el primer metro | **[DECIDIR]** subestima levemente el recorrido y el tiempo de residencia |
| S5 | Propiedades constantes por fase; solidificación a T_m, sin subenfriamiento | simplifica sin perder el orden de magnitud | la urea puede subenfriarse antes de nuclear; eso retrasaría la solidificación |
| S6 | Masa y diámetro constantes durante la caída (densidad del fundido); sin evaporación de agua ni formación de biuret | la humedad del fundido es 0,3 % | no predice porosidad ni la contracción real |
| S7 | Pared adiabática; radiación apagada por defecto | igual que la referencia de validación | — |
| S8 | Humedad del aire solo en el cp | 3,7 % en peso | no calcula secado (lo que sí hace S01) |
| S9 | PSD de Rosin-Rammler con d50 = 1,70 mm y n = 6 en el caso de planta | S01 no informa la PSD; su "0,85 mm" lo interpreto como **radio**, porque sus figuras usan r = 0,00085 m | **[DECIDIR]** si es diámetro, el caso cambia mucho |
| S10 | Un único parámetro de calibración, h_factor, que multiplica el h de la correlación | junta todo lo que el modelo 1D no ve (S2, S3, S4) en un número interpretable | si hace falta h_factor muy lejos de 1 para ajustar, el problema es de estructura, no de parámetro |

## 6. Validación

### 6.1 Verificación numérica (que el código resuelve bien sus ecuaciones)

`tests/test_droplet.py` compara el integrador con soluciones exactas:

| Caso | Referencia | Tolerancia | Resultado |
| --- | --- | --- | --- |
| Conducción en esfera con convección, sin cambio de fase | serie analítica (S09 §5.6): T del centro y T media | 0,01 en θ | pasa |
| Modelo concentrado | exponencial exacta | 0,05 K | pasa |
| Solidificación de una esfera, Ste → 0 | t_f = ρLR/(3hΔT)·(1 + Bi/2) (S08) | 3 % | pasa |
| Conservación de energía con cambio de fase | ∫q dt = m Δh | 1 % | pasa |
| k → ∞ | modelo radial → concentrado | 0,3 K | pasa |
| Velocidad terminal en aire quieto | balance de fuerzas y Stokes | 2 % | pasa |
| Torre | balance global de energía | 0,1 % | pasa (error ≈ 10⁻¹³) |

### 6.2 Validación contra planta (que las ecuaciones representan la torre)

Caso: torre industrial de urea de Torkashvand et al. (2019), Tabla 1 (S01): H = 79,9 m, D = 17 m, 75 t/h de urea a 140 °C, 775 t/h de aire a 35 °C y 3,7 % de humedad.

**Antes de comparar: los datos de planta no cierran el balance de energía.** Con las propiedades del modelo, la urea que entra a 140 °C y sale a 60 °C cede 7412 kW, mientras que el aire que se calienta de 35 a 65 °C gana 6708 kW: un descuadre de 9,5 %. Puede venir de pérdidas por la pared, de evaporación o de dónde se midió cada temperatura. Ningún modelo adiabático puede acertar a la vez la temperatura del producto y la del aire de salida, así que la meta realista es un error del orden del descuadre.

| Variable | Planta | Modelo de los autores | Este modelo, h_factor = 1 | Este modelo, h_factor = 0,61 |
| --- | --- | --- | --- | --- |
| Temperatura del producto (urea) [°C] | 60,0 | 64,5 | 41,8 | 61,5 |
| Aire de salida (cabeza) [°C] | 65,0 | 69,4 | 70,9 | 67,9 |
| Aire a 60 m del fondo (C) [°C] | 61,0 | 63,2 | 52,7 | 55,6 |
| Aire a 40 m del fondo (B) [°C] | 55,0 | 58,3 | 41,7 | 45,8 |
| Aire a 20 m del fondo (A) [°C] | 42,0 | 46,1 | 36,9 | 39,0 |
| **RMSE [°C]** | | **3,8** | **11,3** | **5,2** |

Lectura:

- **Sin calibrar (h_factor = 1)** el modelo enfría demasiado: el producto sale a 41,8 °C contra 60 °C medidos, RMSE 11,3 °C.
- **Calibrando un solo parámetro** (h_factor = 0,61, ajustado minimizando el RMSE de los 5 puntos) el error baja a 5,2 °C. El modelo de los propios autores, con sus parámetros, tiene RMSE 3,8 °C.
- **Ojo con el caso calibrado:** para que el producto salga a 61 °C, el modelo necesita que las dos clases más grandes (≈ 11 % de la masa, d ≥ 2,1 mm) lleguen al fondo con el núcleo todavía a T_m. En una planta que produce prills sanos eso es poco creíble. Es decir, el h_factor ajusta la temperatura correcta por una razón probablemente equivocada. Lo dejo visible en vez de esconderlo con más parámetros; es la limitación que más conviene que domines.
- El perfil de aire medido sube rápido en la parte baja de la torre (35 → 42 °C en los primeros 20 m y → 55 °C a los 40 m): la planta libera más calor abajo de lo que predice cualquier variante del modelo. Las explicaciones candidatas son gotas más grandes que el d50 supuesto (S9), subenfriamiento (S5) y aire canalizado (S2). Es la limitación principal y conviene decirla tal cual.
- Que h_factor resulte ≈ 0,6 es coherente con literatura que reporta coeficientes efectivos menores en nubes densas de gotas, pero eso **[VERIFICAR]**: no tengo una cita que lo respalde con números.

`tests/test_validacion.py` fija estos resultados como regresión.

### 6.3 Pendiente: tu informe de urea

El plan del portafolio pide validar también contra tu informe de urea. Cuando lo tengas, agrega sus datos a `data/raw/` con su `source_id` y un test en `tests/test_validacion.py`. Si es material de un curso o de una empresa, publica solo valores que puedas citar.

## 7. Resultados de sensibilidad (caso calibrado)

- **Tornado ±20 %** (`docs/img/tornado.png`): el tamaño de gota domina: subir d63 un 20 % lleva la T del producto de 61 a 92 °C y deja 41 % de la masa con líquido residual. Le siguen la altura y el h_factor (unos 25 °C de rango cada uno), el calor latente y la T del aire de entrada (17-18 °C). El flujo de aire pesa menos (8,5 °C) y k_s casi nada en la T de mezcla (1,6 °C), aunque sea la propiedad más incierta: es una buena noticia para el modelo y una buena respuesta de entrevista. Los saltos del líquido residual (0 / 4 / 18 %) son gruesos porque la PSD va en 8 clases.
- **Altura mínima para solidificar** (`docs/img/altura_minima.png`): con aire a 35 °C y h_factor = 0,61, la gota más grande que solidifica en 79,9 m es de 2,27 mm (gota aislada en aire a temperatura constante, por lo que es optimista respecto de la torre acoplada, donde el aire abajo ya viene más caliente).
- La altura necesaria crece aproximadamente con d² a d^2,5: la PSD manda. Ver la columna `z_solidificacion_m` en `data/resultados/clases_caso_calibrado.csv`.

## 8. Fuera de alcance

Secado y biuret (S01 los modela; aquí no), aglomeración, formación de polvo y emisiones de la cabeza de la torre, mecánica del balde (solo se usa Rayleigh para d50), CFD de la distribución de aire, y el nitrato de amonio poroso tipo Prillex como caso validado.

## 9. Preguntas que te harán en una entrevista

1. ¿Por qué un modelo radial y no uno concentrado? (Biot, §3.3.)
2. ¿Qué significa físicamente h_factor ≈ 0,61 y por qué no calibraste tres parámetros? (S2, S10: un parámetro interpretable antes que un ajuste perfecto.)
3. Los datos de la planta no cierran el balance: ¿cómo lo detectaste y qué haces con eso? (§6.2.)
4. ¿Qué cambiarías para el nitrato de amonio? (transiciones sólido-sólido, agua, porosidad.)
5. Si la planta empieza a sacar prills blandos en verano, ¿qué te dice el modelo? (tornado: T del aire, flujo de aire, d63; orificios del balde.)
6. ¿Cómo verificaste que el código resuelve bien sus ecuaciones, independiente de si el modelo es bueno? (§6.1: verificación ≠ validación.)
