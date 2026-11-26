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
| h_factor | calibración | – | 0,3-1,5 | 1,00 / 0,67 |
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

**Por qué no un modelo concentrado de tres etapas.** El número de Biot del prill de urea es Bi = h(R/3)/k_s ≈ 0,1-0,2 con k_s = 0,6 W/m K, y Bi = hR/k_s ≈ 0,3-0,6: el núcleo queda bastante más caliente que la superficie. El modelo concentrado (incluido como `simulate_droplet_lumped` para comparar) solidifica la gota de 2,4 mm unos 25 m más arriba que el modelo radial (39 m contra 63 m). Esa diferencia es una buena pregunta de entrevista.

### 3.4 Torre en contracorriente (acoplamiento con el aire)

Volumen de control entre la altura z y el fondo, pared adiabática, aire que entra por abajo:

$$\dot m_a\,[h_a(T_a(z)) - h_a(T_{a,in})] = \dot m_f \sum_j w_j\,[\bar H_j(z) - \bar H_j(H)]$$

donde w_j es la fracción másica de la clase j y H̄_j(z) su entalpía media en la altura z. Se itera por punto fijo (relajación 0,7) sobre el perfil T_a(z) hasta que cambie menos de 0,05 K; converge en 5 a 8 iteraciones. La entalpía del aire incluye la humedad solo en el cp: c_p = (c_p,seco + Y c_p,vapor)/(1 + Y).

### 3.5 PSD

Distribución acumulada en masa de Rosin-Rammler (S10): F(d) = 1 − exp[−(d/d63)ⁿ], discretizada en 8-12 clases entre los cuantiles 0,5 % y 99,5 %. Opcionalmente, d50 a partir del orificio del balde por ruptura de Rayleigh (S07): d = 1,89 d_orificio. El diámetro del prill sólido es d_s = d (ρ_l/ρ_s)^{1/3} (masa constante, sin porosidad).

## 4. Propiedades

| Propiedad | Urea, valor usado | Fuente | Contraste | Comentario |
| --- | --- | --- | --- | --- |
| T_m | 132,7 °C | S01 Tabla 2 | `thermo` 133,4; S11 132,6 | concuerdan |
| L | 224 kJ/kg | S01 Tabla 2; S12 Tabla II (de Ullmann) | `thermo` 231; S11 13,61 kJ/mol = 227 | concuerdan (3 %). S11 imprime también 251 J/g en una tabla, inconsistente con su propio texto: descartado |
| ρ_l | 1245 kg/m³ | S11: ρ = 1638,5 − 0,96 T [K], válida 133-150 °C | S01 y S12: 1220; `thermo`: 799 | `thermo` subestima 35 %: **no usar** |
| ρ_s | 1335 kg/m³ (20 °C) | S11 | S13: 1335; `thermo`: 900 | ídem |
| c_p,l | 2,25 kJ/kg K | S11 (135 °C) | S13: 2,25; `thermo`: 1,92 | `thermo` subestima 15 % |
| c_p,s | 1,78 kJ/kg K | media de la correlación de S11 (240-400 K) entre 35 y 132,7 °C | a 298 K la correlación da 93,1 J/mol K; NIST: 92,8 | S01 y S12 usan 1,334 para ambas fases, sin distinguir |
| k_l | 0,31 W/m K | `thermo` (estimación) **[VERIFICAR]** | S13 imprime 0,013: imposible (menor que el aire) | no encontré una medición accesible |
| k_s | 0,60 W/m K | **[VERIFICAR] supuesto** | ver abajo | pesa poco: ±20 % mueve la T del producto 1,6 °C |
| ε | 0,9 | **[VERIFICAR] supuesto** | | solo si se activa radiación |

**La conductividad del sólido merece un párrafo, porque te la van a preguntar.** La literatura de prilling (S12 Tabla II, citando a Ullmann; S01 Tabla 2; S13 Tabla 1) repite el valor 2,651·10⁻⁵ kW/m K = 0,0265 W/m K. Ese número no puede ser la conductividad de un cristal molecular: es la del aire. Lo comprobé con el modelo: con k_s = 0,0265 una gota de 1,6 mm necesita 107 m de caída para solidificar en aire a 35 °C (24 m con k_s = 0,6), y las torres industriales de 50-80 m sí producen prills de 1,6 mm sólidos (S12, S13). La única medición directa que encontré es Andersson y Ross (1994, S14), de pago. Mientras no la tengas, el modelo usa 0,6 W/m K, que está en el rango de los sólidos orgánicos cristalinos. Si puedes conseguir S14 por la universidad, ese es el dato que falta.

Aire: CoolProp, aire seco pseudo-puro a 1 atm, tabulado 240-480 K e interpolado. Hallazgo útil para la entrevista: `thermo` da la densidad de la urea fundida y del nitrato de amonio fundido con errores de 30-35 % y el cp del fundido 15 % bajo; por eso el modelo toma esos valores de literatura y deja `thermo` solo como contraste (pestaña Propiedades de la app).

### 4.1 Nitrato de amonio (opción "no validado")

La app permite elegir nitrato de amonio para conectar con productos como el prill poroso de nitrato que se usa en explosivos para minería, pero el modelo no está hecho para él y lo dice. Le faltan tres cosas:

1. **Transiciones sólido-sólido.** El nitrato de amonio cambia de estructura cristalina varias veces al enfriarse (fases I, II, III y IV; una cerca de 125 °C y otra cerca de 32 °C) y cada cambio libera calor y cambia el volumen. Al pasar repetidamente por 32 °C en acopio, el prill se agrieta y se vuelve polvo. Modelarlo exige agregar esas entalpías como "mesetas" adicionales en h(T), igual que la de fusión.
2. **Agua en el fundido.** El prill poroso se hace con una solución concentrada (del orden de 95-97 %), no con fundido puro. La porosidad, que es la gracia del producto porque absorbe el combustible, viene del agua que se evapora en la torre y en el secado posterior. El modelo no tiene transferencia de masa.
3. **Propiedades.** Varias vienen de `thermo` y al menos una está claramente mal (densidad del fundido); están marcadas [VERIFICAR] en el código.

Recomendación: no lo muestres como resultado en el video; menciónalo como extensión natural ("el siguiente paso sería el nitrato de amonio, que agrega transiciones de fase y secado").

## 5. Supuestos

| # | Supuesto | Por qué es razonable | Cuándo deja de serlo |
| --- | --- | --- | --- |
| S1 | Estado estacionario de la torre | operación continua | arranques, cambios de carga |
| S2 | Flujo de aire uniforme en la sección (tapón, 1D) | ver §5.1 | ver §5.1 |
| S3 | Gotas esféricas y sin interacción (sin coalescencia ni estela) | la fracción volumétrica de gotas en la torre es ≈ 0,1 % (S12, supuesto 7), así que cada gota "ve" aire libre | cerca del balde la nube es más densa |
| S4 | Solo movimiento vertical; la gota parte con v₀ = 0 | ver §5.1 | ver §5.1 |
| S5 | Propiedades constantes por fase; solidificación a T_m, sin subenfriamiento | simplifica sin perder el orden de magnitud | la urea puede subenfriarse antes de nuclear; eso retrasaría la solidificación |
| S6 | Masa y diámetro constantes durante la caída (densidad del fundido); sin evaporación de agua ni formación de biuret | la humedad del fundido es 0,3 % | no predice porosidad ni la contracción real |
| S7 | Pared adiabática; radiación apagada por defecto | igual que la referencia de validación; S12 estima la radiación en 0,6 % del calor | — |
| S8 | Humedad del aire solo en el cp | 3,7 % en peso | no calcula secado (lo que sí hace S01) |
| S9 | PSD de Rosin-Rammler con d50 = 1,70 mm y n = 6 en el caso de planta | S01 no informa la PSD. Su "0,85 mm" es un **radio**: el texto lo escribe como "R_p = 0,85 mm" y las figuras 2 y 3 usan r = 0,00085 m (S01, p. 176-177). El n = 6 es supuesto | si la PSD real es más ancha, aumenta la fracción blanda |
| S10 | Un único parámetro de calibración, h_factor, que multiplica el h de la correlación | junta lo que el modelo 1D no ve (S2, S4, S5) en un número interpretable | si hace falta h_factor muy lejos de 1 para ajustar, el problema es de estructura, no de parámetro |

### 5.1 Por qué aire 1D uniforme y caída solo vertical

**Aire en flujo pistón uniforme (S2).** Es la elección estándar de los modelos de torre que sirven para diseñar y operar, por tres razones:

- **Es lo que se mide en planta.** La planta conoce el flujo total de aire y, con suerte, temperaturas a algunas alturas (S01 Tabla 1; S12 informa una velocidad media de 0,63 m/s medida por la empresa). No hay datos para alimentar un perfil radial, así que un modelo 2D tendría más parámetros que datos.
- **El acoplamiento que importa es axial.** La pregunta de diseño es cuánta altura hace falta, y eso depende de cómo sube la temperatura del aire con la altura, que es justo lo que resuelve el balance 1D de §3.4.
- **El costo es bajo y conocido.** Una corrida de la torre tarda segundos, lo que permite el análisis de sensibilidad y la app. Un CFD euleriano-lagrangiano (como S16) tarda horas por caso y igual usa las mismas correlaciones de Ranz-Marshall.

Lo que se pierde: en una torre real el aire entra por persianas en la base y no se reparte parejo; parte sube canalizado por las paredes y las gotas del centro ven aire más caliente. Ese efecto es la principal sospecha detrás de h_factor < 1 (§6.2), y por eso no se esconde en otro parámetro.

**Caída solo vertical, partiendo de v₀ = 0 (S4).** El balde rotatorio lanza la gota casi en horizontal, a unos 5-10 m/s, y el arrastre la frena rápido. Lo calculé con el mismo arrastre del modelo en 2D para una gota de 1,7 mm: lanzada a 5-10 m/s, su velocidad horizontal baja a menos del 10 % en unos 2 s, mientras cae unos 10 m y se desplaza 4 a 7 m hacia la pared. Desde ahí cae vertical a su velocidad terminal, que es lo que el modelo supone desde el inicio. Se ignora porque:

- simplifica el modelo a una dimensión, consistente con S2, y deja la trayectoria en función de un solo dato (la altura);
- el error va hacia el lado conservador: en esos primeros 10 m la velocidad relativa real es mayor, el h también y la gota está más caliente, así que el modelo subestima algo el enfriamiento y sobreestima la altura necesaria;
- el desplazamiento horizontal (4-7 m) es menor que el radio de la torre (8,5 m en el caso de planta), así que las gotas no chocan con la pared, que es la condición para que el supuesto valga.

Lo que se pierde: la gota sí recorre ese tramo inicial a alta velocidad relativa, con un h mayor justo cuando está más caliente. Si quisieras afinarlo, el cambio es pasar a velocidad 2D [v_x, v_z] con el mismo arrastre, sin cambiar nada más.

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

Caso: torre industrial de urea de Torkashvand et al. (2019), Tabla 1 (S01): H = 79,9 m, D = 17 m, 75 t/h de urea a 140 °C, 775 t/h de aire a 35 °C y 3,7 % de humedad. El PDF está en acceso abierto: https://doi.org/10.2298/CICEQ180827030T (Tabla 1 en la p. 176).

**Antes de comparar: los datos de planta no cierran el balance de energía.** Con las propiedades del modelo, la urea que entra a 140 °C y sale a 60 °C cede 7700 kW, mientras que el aire que se calienta de 35 a 65 °C gana 6708 kW: un descuadre de 12,9 %. Puede venir de pérdidas por la pared, de evaporación o de dónde se midió cada temperatura. Ningún modelo adiabático puede acertar a la vez la temperatura del producto y la del aire de salida, así que la meta realista es un error del orden del descuadre.

| Variable | Planta | Modelo de los autores | Este modelo, h_factor = 1 | Este modelo, h_factor = 0,67 |
| --- | --- | --- | --- | --- |
| Temperatura del producto (urea) [°C] | 60,0 | 64,5 | 44,4 | 61,5 |
| Aire de salida (cabeza) [°C] | 65,0 | 69,4 | 72,0 | 69,2 |
| Aire a 60 m del fondo (C) [°C] | 61,0 | 63,2 | 54,4 | 56,4 |
| Aire a 40 m del fondo (B) [°C] | 55,0 | 58,3 | 43,1 | 46,2 |
| Aire a 20 m del fondo (A) [°C] | 42,0 | 46,1 | 37,4 | 39,2 |
| **RMSE [°C]** | | **3,8** | **10,0** | **5,0** |

Lectura:

- **Sin calibrar (h_factor = 1)** el modelo enfría demasiado: el producto sale a 44,4 °C contra 60 °C medidos, RMSE 10,0 °C.
- **Calibrando un solo parámetro** (h_factor = 0,67, ajustado minimizando el RMSE de los 5 puntos) el error baja a 5,0 °C. El modelo de los propios autores, con sus parámetros, tiene RMSE 3,8 °C.
- **Ojo con el caso calibrado:** para que el producto salga a 61 °C, el modelo necesita que las gotas sobre 2,14 mm (6,3 % de la masa) lleguen al fondo con el núcleo todavía a T_m. En una planta que produce prills sanos eso es dudoso, aunque la literatura reciente sí describe solidificación parcial de las gotas grandes en torres reales (S13). El h_factor puede estar ajustando la temperatura correcta por una razón equivocada. Lo dejo visible en vez de esconderlo con más parámetros; es la limitación que más conviene que domines. (En la primera versión este número era 11 %: la PSD en 10 clases lo exageraba. Ahora la fracción blanda sale del diámetro crítico sobre la PSD continua, ver §7.)
- El perfil de aire medido sube rápido en la parte baja de la torre (35 → 42 °C en los primeros 20 m y → 55 °C a los 40 m): la planta libera más calor abajo de lo que predice cualquier variante del modelo.
- Que h_factor resulte < 1 **no** se puede atribuir a la interacción entre gotas: en estas torres la fracción volumétrica de gotas es ≈ 0,1 % y la literatura la desprecia explícitamente (S12, supuesto 7). Busqué una cita que respaldara un h efectivo menor por nube de gotas y no la encontré para estas condiciones, así que la explicación queda entre la mala distribución del aire (S2), una PSD más gruesa que la supuesta (S9) y el subenfriamiento (S5).

`tests/test_validacion.py` fija estos resultados como regresión.

### 6.3 Pendiente: tu informe de urea

El plan del portafolio pide validar también contra tu informe de urea. Cuando lo tengas, agrega sus datos a `data/raw/` con su `source_id` y un test en `tests/test_validacion.py`. Si es material de un curso o de una empresa, publica solo valores que puedas citar. Alternativa si no aparece: Mehrez et al. (2012, S12) publican otra torre industrial (57 m, 18 m de diámetro, prill medio de 1,6 mm por tamizado), aunque sin la temperatura del aire de entrada en su tabla.

## 7. Planta completa: tamiz, reciclo y enfriador

La torre no entrega el producto final. Aguas abajo hay un tamiz y, si el producto sale caliente, un enfriador. El módulo `src/prilling/process.py` resuelve el conjunto en estado estacionario:

```
fundido fresco ─► (+ reciclo redisuelto) ─► TORRE ─► TAMIZ ─► ENFRIADOR ─► producto
                        ▲                              │
                        └── finos + gruesos + blandos ◄┘
```

**Diámetro crítico.** La gota que solidifica justo al llegar al fondo, d_crit, se obtiene ajustando z_solid = a·d^b a las clases sólidas más grandes y resolviendo z_solid = H. Todo lo que está sobre d_crit en la PSD continua es "blando". Esto evita que el resultado salte de clase en clase.

**Tamiz (P2).** Corte ideal por diámetro del prill sólido en [1,0; 2,4] mm **[DECIDIR: rango de especificación]**. Los prills blandos se cuentan como fuera de especificación aunque estén en rango, porque se deforman y se pegan en la malla.

**Reciclo (P1).** Lo fuera de especificación no puede volver como sólido a la boquilla: el balde rotatorio solo pulveriza líquido. En las plantas de urea los finos y gruesos se **redisuelven** (en solución de urea, aguas arriba de la evaporación) y vuelven como fundido. Por eso el modelo lo suma a la alimentación de la torre: m_torre = m_fresco / f_en_especificación. ¿Y una boquilla aparte para el reciclo? Solo tiene sentido si ese reciclo ya volvió a ser fundido, y entonces conviene que esa boquilla haga gotas más chicas que la principal. Pero es un segundo sistema de pulverización para un 10 % del flujo, y refundir urea agrega tiempo a alta temperatura, que es justo lo que forma biuret (S01 lo limita a < 1 %). Redisolver en la corriente principal es la práctica habitual; la boquilla aparte queda como idea, no como diseño.

**Enfriador de lecho fluidizado (P3).** Lecho bien mezclado donde el aire sale a la temperatura del lecho, porque con partículas de ~2 mm el gas se equilibra a pocos centímetros del distribuidor (S15). Así el enfriador lo limita el balance de energía: m_s c_p,s (T_entrada − T) = m_a c_p,a (T − T_a,entrada). Objetivo 50 °C **[DECIDIR: temperatura de acopio o ensacado]**. El enfriador resuelve producto **caliente pero sólido**; no arregla un prill blando, que ya se deformó en la torre.

**Resultados en el caso calibrado (h_factor = 0,67):**

| | Sin reciclo | Con reciclo | Alternativa: gota más chica, sin reciclo |
| --- | --- | --- | --- |
| d50 de gota | 1,70 mm | 1,70 mm | 1,57 mm (orificio ≈ 0,83 mm por Rayleigh) |
| Carga de la torre | 75 t/h | 84,0 t/h (12 % de reciclo) | 75 t/h |
| Prills blandos | 6,3 % | 7,4 % | 1,0 % |
| Finos (< 1 mm) | 3,2 % | 3,2 % | 5,2 % |
| En especificación | | 89,3 % | 93,8 % |
| T a la salida del tamiz | | 58,4 °C | 52,1 °C |
| Aire del enfriador para llegar a 50 °C (aire a 35 °C) | | 74 t/h | 23 t/h |

Lectura, que es la parte que vale en la entrevista:

- **El reciclo empeora la torre.** La torre recibe 12 % más de fundido con el mismo aire, el aire se calienta más y la fracción blanda sube de 6,3 a 7,4 %. El reciclo no resuelve la causa; solo la recircula.
- **La palanca de fondo es el tamaño de gota.** Bajar el d50 de 1,70 a 1,57 mm (orificios de ~0,9 a ~0,83 mm) deja menos de 1 % blando, sube lo que está en especificación y reduce a un tercio el aire del enfriador, a costa de más finos (de 3,2 a 5,2 %). Los finos sí son sólidos y se pueden redisolver sin dañar la torre.
- **¿Y si el producto llega a 92 °C** porque se agrandó la gota? En el tornado, subir d63 un 20 % lleva el producto a 90 °C, pero también deja 40 % de la masa blanda. Ahí el enfriador no sirve: el problema es que la gota no alcanzó a solidificar. Un enfriador se justifica cuando el producto sale sólido pero caliente, por ejemplo en verano con el aire de entrada a 40-45 °C.

## 8. Resultados de sensibilidad (caso calibrado)

- **Tornado ±20 %** (`docs/img/tornado.png`): el tamaño de gota domina. Subir d63 un 20 % lleva la T del producto de 62 a 90 °C y deja 40 % de la masa con líquido residual. Le siguen la altura y el h_factor (unos 24 °C de rango cada uno), el calor latente y la T del aire de entrada (17 °C). El flujo de aire pesa menos (8 °C). k_s casi no mueve la T de mezcla (1,6 °C), pero sí la fracción blanda (5,0 a 8,9 %): por eso vale la pena conseguir el dato de S14.
- **Altura mínima para solidificar** (`docs/img/altura_minima.png`): con aire a 35 °C y h_factor = 0,67, la gota más grande que solidifica en 79,9 m es de 2,28 mm (gota aislada en aire a temperatura constante, optimista respecto de la torre acoplada, donde el aire abajo ya viene más caliente: ahí el diámetro crítico es 2,14 mm).
- La altura necesaria crece aproximadamente con d² a d^2,5: la PSD manda.

## 9. Fuera de alcance

Secado y biuret (S01 los modela; aquí no), aglomeración, formación de polvo y emisiones de la cabeza de la torre, mecánica del balde (solo se usa Rayleigh para d50), CFD de la distribución de aire, dimensionamiento mecánico del tamiz y del enfriador, y el nitrato de amonio poroso como caso validado (§4.1).

## 10. Preguntas que te harán en una entrevista

1. ¿Por qué un modelo radial y no uno concentrado? (Biot, §3.3.)
2. ¿Qué significa físicamente h_factor ≈ 0,67 y por qué no calibraste tres parámetros? (S2, S10: un parámetro interpretable antes que un ajuste perfecto; no es interacción entre gotas, S3.)
3. Los datos de la planta no cierran el balance: ¿cómo lo detectaste y qué haces con eso? (§6.2.)
4. ¿De dónde sacaste la conductividad del sólido? (§4: el valor de la literatura de prilling es el del aire, y lo demostraste con el modelo.)
5. Si la planta empieza a sacar prills blandos en verano, ¿qué haces? (tornado y §7: primero tamaño de gota o aire de la torre; el reciclo y el enfriador no lo arreglan.)
6. ¿Por qué no reciclas directo a la boquilla? (§7: el balde solo pulveriza líquido; redisolver y el biuret.)
7. ¿Qué cambiarías para el nitrato de amonio? (§4.1.)
8. ¿Cómo verificaste que el código resuelve bien sus ecuaciones, independiente de si el modelo es bueno? (§6.1: verificación ≠ validación.)
