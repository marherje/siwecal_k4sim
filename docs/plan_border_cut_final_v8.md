# Plan: `barxy_min_border_distance`, campaña final_v8 y paneles con corte de contención lateral

Estado a 2026-10-28. Retomado con aprobación del usuario de toda la hoja A–E (test_v8 → final_v8 → data → plots →
artefacto; un rescue DAG como máximo por DAG; sin condor_rm ni borrados en EOS).
- 22:15 test_v8 OK (Simulation/tests/border_cut/test_v8): beam z −65500 → TimeOffsetNs 209.65 ns, hits/ev digi 223
  (v7 235), C++ = Python 5.6e-6 mm.
- 22:2x **final_v8 enviada**: cluster 21355387, DAG analysis/condor_final_v8 (2350 chunks, 94 grupos).
- Data: prueba de un run (run 20) OK: 84 ramas escalares + 12 por hit idénticas a Reconstruction, 3.2 GB, 17 min.
- **Corrección del usuario:** la data reprocesada sustituye a Reconstruction (solo la simulación se versiona). El DAG
  `border_data` (solo valtree) se paró con condor_rm; enviado **`reco_border`** (cluster 21357359,
  analysis/condor_reco_border, 103 runs, `--format both --validation --hit-mip-cut 0.5`, 8 GB) con staging en
  Simulation/tests/border_cut/data/<run>/. Cuando termine y esté verificado, el usuario ejecuta
  `analysis/replace_reconstruction_pid.sh --doit` (sobrescribe edm4hep + valtree en Reconstruction/<run>/; decidió
  sin respaldo). A partir de ahí los plots usan FINAL_RECO por defecto (Reconstruction) y FINAL_SIM=final_v8.

## Objetivo
Rehacer, sobre una campaña nueva (**final_v8**: línea H2 + sag k·Q, 94 puntos × 50k) y con
`barxy_min_border_distance > 50 mm` (más el corte de hit de 0.5 MIP de siempre), los perfiles longitudinales en X0,
los paneles resumidos, los paneles por grupo y la resolución + linealidad del artefacto
*SiW-ECAL TB2026 Final Sample* (https://claude.ai/artifact/MdY3nMKMUyDs4cDWk6geG5), en una sección nueva.
final_v7 se conserva en `Processed/final_v7` para comparar.

## Fondos π⁻/μ⁻ con la configuración de v8 (2026-10-29)
`final_v8_bkg` (cluster 21368919, analysis/condor_final_v8_bkg): los 4 puntos de bkg_points.txt (π⁻ y μ⁻ en E1 y E2,
th230, energía plana 5–100 GeV), 25 × 2000 cada uno, seed offset 40000000, H2 + TimeOffsetNs + k·Q por defecto.
final_v7_bkg se conserva. Mismo uso que final_v7_bkg (fondos generales para limpiar los runs de positrones), ahora
coherentes con la configuración de v8.

## Corte de borde: 50 mm (decidido 2026-10-29, sustituye a 40 mm en todo)
Pérdida lateral medida (energía W-ponderada más allá de d del barycentro, a un lado; th230 E4, final_v7 / data):
d = 40 mm 1.1–1.4 %, **d = 50 mm 0.6–0.8 %**, d = 30 mm 2.2–2.5 %; casi independiente de E, se divide por 2 cada
~10 mm. Supervivencia en data (shower-like, bar_x/bar_y de Reconstruction, scratchpad border_survival.py): con 50 mm
quedan < 1000 sucesos en th210 34/52/74 GeV y 99 GeV (22, −48) (618/277/238/67); el scan th210 se queda en 2–27 %
(barycentro a y ≈ −50, solo pasa la cola del haz), th230 E1 en 2–14 %; E4/E5 conservan 50–100 %. La recta MIP→GeV se
ajusta sin corte, así que la calibración no se ve afectada; la resolución/linealidad de th210 tendrá poca estadística.

## Sag medido: mapa B(Q_chip, Q_slab) por umbral → campaña final_v9 (2026-10-29, noche)
- Datos: analysis/cm_sag_collect.py (Q_chip, Q_slab, Q_det, carga previa de la adquisición) → EOS
  sim_tuning/occupancy/sag3/cm_sag_<run>.npz (8 runs, 1.8 M registros). Ajustes: sag_model_fit.py, sag_acq_fit.py.
  Ninguna fórmula aguanta con Q_slab alto (satura a −10…−14 ADC); el término de slab es esencial; por chip, Q_det y la
  carga previa no aportan. **Mapa medido** 16 × 11 por umbral: analysis/sag_table.py → mappings/occupancy_sag_<th>.txt
  (residuo por run ≤ ±0.2 ADC; χ²/ndf 3–8 frente a 2000–7000 de k·Q).
- AdcDigitizer: propiedad OccupancySagTable (bilineal en ln(1+Q), fijado en los bordes; Q_slab = otros chips del
  slab); exclusiva con OccupancySlope. digi_calibration.yml: occupancy_slope 0.0 + occupancy_sag_table por umbral
  (y se limpiaron las claves occupancy_slope duplicadas). job3_digitize.py: OCCUPANCY_SAG_TABLE la sobreescribe.
  k4sim recompilado 19:33.
- final_v8 y final_v8_bkg (k·Q) cancelados con condor_rm (decisión del usuario: relanzar con la corrección). Sus
  carpetas Processed/final_v8{,_bkg} quedan a medio llenar con el modelo viejo: el usuario decide si borrarlas.
- Prueba test_v9 (cluster 21377739, 4 puntos × 4000): a 20 GeV hits/ev sim/data 0.988 (th230 E4) y 1.032 (th210),
  como sim_v7 (el sag pesa poco a 20 GeV; el +3 % de th210 no es sag). 74 GeV pendiente.
- **final_v9 lanzada** 2026-10-29 21:50 (cluster 21378160, analysis/condor_final_v9, 2350 chunks + 94 grupos) y
  **final_v9_bkg** (cluster 21378161). El usuario pidió borrar todo v8: comandos dados (rm de Processed/final_v8{,_bkg},
  tests/border_cut/test_v8 y los condor_*v8 locales; condor_val_v8 NO es la campaña v8).
- Datos (reco_border2) terminado 20:12; 4 runs (107, 111, 116, 134) tenían sucesos vacíos que run_pid_batch descarta
  → rehechos con la receta de Reconstruction (analysis/border_data_node_reco.sh, cluster 21378634) para que el
  valtree vaya fila a fila con ecal_<run>.root. Después: replace_reconstruction_pid.sh (comprobación) → el usuario
  ejecuta --doit.
- Plots con FINAL_SIM=…/Processed/final_v9.

## Estado 2026-10-31 13:00
- final_v9 TERMINADA (31 Oct 04:01, 0 fallos, 188 pid), final_v9_bkg y test_v9 terminadas. test_v9 74 GeV: hits/ev
  sim/data th230 E4 0.996 (v7 1.077, k·Q 1.031), th210 1.072 (v7 1.127); E/hit +3.5 % th230, +1.4 % th210.
- 10 runs de datos con sucesos vacíos (12, 56, 72, 78, 93, 94, 107, 111, 116, 134) rehechos con
  analysis/border_data_node_reco.sh; los 103 pasan replace_reconstruction_pid.sh (comprobación). Pendiente del
  usuario: --doit, borrar staging, borrar v8, commits (lista en la conversación).
- final_profiles.py: el cruce valtree ↔ ecal tree ahora tolera sucesos ausentes en el valtree (run_pid_batch descarta
  los vacíos; pasa en sim nodigi), siguen en orden y cada fila del valtree debe emparejarse.
- validate_v8.py: --cache-tag (caché por campaña; "new" tenía la de la mini-validación) y --label.
- **Job de plots** analysis/final_v9_plots_job.sh (cluster 21410521) → EOS tests/border_cut/figs_v9/{nocut, b50,
  control_v7, validation}, datos desde el staging tests/border_cut/data. Probado en local con 3 grupos.
- Siguiente: revisar figuras; control_v7 debe dar a = 81.9/87.2/88.7; sección nueva del artefacto final
  (MdY3nMKMUyDs4cDWk6geG5, leerlo entero antes) con b50 + supervivencia + resolución/linealidad; sección 3 de la nota
  de upstream (FZqHXKJ…) con la validación de v9.

## Forma del sag de línea de base (estudio 2026-10-29, scratchpad sag_*.py sobre occupancy/cm/cm_trigger_e*.npz)
Ajuste de B (canales sin señal) frente a Q del chip-SCA, runs th210 e291/e287/e286 (20/52/74 GeV):
- B = kQ por el origen (el digitizador) es lo peor (χ²/ndf 2000–3300). La "medida de pedestales" −0.95e-3 es
  B = −0.853 − 0.953e-3·Q ajustada a los tres runs juntos. Cuadrática: no mejora el combinado (700 vs 706) y gira
  dentro del rango medido (Q ≈ 4–7k): descartada. Por run, a + b·ln(1+Q/q0) o a + b·Q^p ajustan (χ²/ndf 15–80);
  forma común con offset libre por run: 270–400 → cambia la forma entre runs.
- **Por slab:** el offset a varía de −0.2 a −3.5 ADC, máximo en el máximo de la cascada (slabs 4–5) y crece con la
  energía en el mismo slab (slab 5: −1.5/−3.2/−3.5); la pendiente local k ≈ −0.5…−1.0e-3 por slab, −0.86e-3 mediana
  por chip (±0.45e-3); slab 12 (COB) −1.0…−1.7e-3. → dos términos: local (Q del chip) + de capa (carga del slab).
  k = −1.5e-3 por el origen absorbe el término de capa; por eso un único k no vale entre energías/umbrales.
- Siguiente: medir B frente a Q_chip y Q_slab juntos, por chip, y poner B = k_chip·Q_chip + f(Q_slab) en el
  AdcDigitizer. El k por umbral es la solución provisional.

## Corrección de ocupación por umbral (decidido 2026-10-29)
El scan de k (nota de upstream, figura kscan) solo sostiene k = −1.5e-3 para **th230** (4 energías); th210 tenía un
único punto (20 GeV) que prefiere otro valor y th220 ninguno. El usuario quiere resultados consistentes: **un k por
umbral** (`occupancy_slope` por set en digi_calibration.yml ya lo admite). final_v8 corre con −1.5e-3 en los tres:
los grupos th230 (32) valen; los th210 (11) y th220 (51) habrá que rehacerlos con su k cuando esté ajustado.
- Figura: `beamline_occupancy_figs.fig_kscan` escribe ya `kscan_<set>.png`; en la nota solo th230 por ahora.
- Pendiente tras final_v8: scan de k para th210 (10, 20, 52, 74 GeV) y th220 (20 GeV y varias posiciones a 74 GeV).
  Las muestras de final_v8 llevan k = −1.5e-3 dentro: el scan post-hoc solo llega a k más negativos (aprox.);
  para cubrir ambos lados hacen falta muestras con OCCUPANCY_SLOPE=0 (DAG pequeño). Medir también la pendiente del sag
  en datos de th230/th220 (cm_trigger_test.py), hasta ahora solo en runs de th210.

## Decisiones tomadas
- Borde = área activa, L = 89.09 mm (centro del pad extremo 86.325 + medio pad 5.53/2).
  `barxy_min_border_distance = min(L − |bar_x|, L − |bar_y|)` [mm]; negativo fuera, NaN sin barycentro.
- **Calibración frente a selección:** la recta MIP→GeV se ajusta **sin** corte de borde en el scan th210 y se aplica
  a todo; el corte de borde solo decide qué sucesos entran en paneles, resolución y linealidad.
- Simulación: final_v8. Artefacto: sección nueva, las actuales no se tocan.

## Hecho
### siwecal-tb2026 (rama `feat/border-distance`, sin commit)
- C++: `EcalShowerVars.h` (campo, nombre al final de `scalarNames()`, cálculo tras `bar_r`, L en
  `ShowerThresholds`), `EcalShowerVariables.cpp` (`appendScalars`, propiedades `ActiveHalfSizeX/Y`).
- Python: `metrics.min_border_distance` (+ copia en `event_viewer/_metrics.py`), `event_data.py` (los PID antiguos
  leen NaN), `SCALAR_NAMES` (`siwecal_common/edm4hep_pid.py` y copia), corte `border_min/max` en `selection.py`,
  `CACHE_FORMAT_VERSION` 5 → 6, test `test_min_border_distance`, READMEs.
- Compilado `gaudi_source/build` (2026-10-28 21:16, `make -j2`).
- Verificado en 2001 sucesos de `final_v7/trees/ecal_th230_E4_e20.root`:
  C++ = Python desde `bar_x/bar_y` (máx. 7.5e-6 mm); las 48 ramas antiguas idénticas bit a bit a final_v7;
  mediana 54.6 mm, 88 % > 40 mm. `pytest siwecal_validation/tests`: 28/28.

### siwecal_k4sim (`test/new-studies`, sin commit)
- Pendiente de commit 1 (nota de upstream/ocupación): geometría H2 (`SND_compact_beamline.xml`, + `_v7.xml` para
  final_v5–v7), `materials.xml`, `RealDigitizer` TimeOffsetNs, `AdcDigitizer` k·Q, `job3_digitize.py`,
  `digi_calibration.yml` (occupancy_slope −1.5e-3), `make_sim_dag.py` (`--compact`, `--beam-z` −65500 y ahora
  `--campaign` obligatorio), lanzadores, `docs/gaudi_pipeline.md`.
- Commit 2: copias de la lista (`event_viewer/*`, `analysis/tests/test_k4siwecalreco_output.py`,
  `analysis/run_sim_chain.sh`) y en `final_panels.py`, `final_panels_select.py`, `final_profiles.py`:
  `--border-cut MM` (tras `showerlike()`; el p90 se sigue calculando sobre toda la muestra), `--cal-file`
  (obligatorio con el corte), `mip_gev_line.txt`, `border_survival.txt`, caché `cache_<campaign>[_b<MM>]/`.
  `final_profiles.py` une ecal tree y valtree fila a fila comprobando (run, event, spill, bcid) en cada fila (en sim
  la clave no es única: los chunks repiten numeración); `FINAL_RECO_TREES` para los ecal trees de data.
- Tests k4sim (con `source init_siwecal_soft.sh`): todo pasa salvo `test_scalar_variable_present[fractal_dimension]`
  y `[barxy_min_border_distance]`: el fichero de referencia `gaudi_jobs/1_mu_beam_pipeline/ecal_sim.edm4hep.root`
  es del 2026-09-07, anterior a ambas variables (ya fallaba con fractal_dimension). Regenerarlo con la build nueva.

### Commits (los hace el usuario; ver memoria `no-commits`)
```
cd ~/public/siwecal-tb2026
git add README.md gaudi_source/README.md gaudi_source/include/k4SiWEcalReco/EcalShowerVars.h \
  gaudi_source/src/components/EcalShowerVariables.cpp siwecal_common/edm4hep_pid.py \
  siwecal_validation/{event_data,metrics,selection,vars_cache}.py siwecal_validation/tests/test_metrics.py \
  event_viewer/_edm4hep_pid.py event_viewer/_metrics.py
git commit -m "barxy_min_border_distance: barycentre distance to the active edge as a PID scalar"

cd ~/public/siwecal_k4sim
git add analysis/make_sim_dag.py docs/gaudi_pipeline.md gaudi_jobs/pid2026_common/job3_digitize.py \
  gaudi_source/AdcDigitizer.cpp gaudi_source/RealDigitizer.cpp mappings/digi_calibration.yml \
  simulation/geometry/SND_compact_beamline.xml simulation/geometry/SND_compact_beamline_v7.xml \
  simulation/geometry/materials.xml simulation/run_script/generic_condor_beam.sh \
  simulation/run_script/launch_beam_final.sh
git commit -m "H2-like upstream line, TimeOffsetNs and k*Q occupancy sag as simulation defaults; --campaign required"

git add event_viewer/_edm4hep_pid.py event_viewer/_metrics.py analysis/tests/test_k4siwecalreco_output.py \
  analysis/run_sim_chain.sh analysis/final_panels.py analysis/final_panels_select.py analysis/final_profiles.py \
  docs/plan_border_cut_final_v8.md
git commit -m "barxy_min_border_distance mirrors; --border-cut and MIP->GeV line file in the final plot scripts"
```
**No recompilar `gaudi_source/build` mientras corra la campaña** (los jobs cargan esa `.so` desde AFS).

## Pendiente
1. **Prueba pequeña** (confirmación del usuario antes de enviar):
   `python3 analysis/make_sim_dag.py --campaign test_v8 --only th230_E4_e20 --nevents 20 --chunk-events 10`.
   El `--submit` crea carpetas en EOS (`eos_mkdir`) y el sandbox no puede escribir en EOS: generar el DAG sin
   `--submit`, que el usuario cree `Processed/test_v8/{,chunks/,chunks/th230_E4_e20/}` y lanzar
   `condor_submit_dag -batch-name test_v8 test_v8.dag` desde `analysis/condor_test_v8/` (comando suelto, sin
   `cd … &&`, para que `excludedCommands` lo saque del sandbox). Comprobar: offset de 209.65 ns en el log de
   digitización, hay hits, el valtree tiene `barxy_min_border_distance`, valor = cálculo desde `bar_x/bar_y`.
   (Mejor salida para test: `Simulation/tests/border_cut/test_v8` vía `--out`.)
2. **final_v8**: `make_sim_dag.py --campaign final_v8` (94 puntos × 50k, ~2 días, ~400 GB; cuota siw-ecal ~52 % →
   ~64 %). Seguimiento: `analysis/sim_dag_status.sh`, `condor_q` y los logs de `analysis/condor_final_v8/`.
3. **Data, en paralelo**: `run_pid_batch.py --format valtree --hit-mip-cut 0.5` sobre
   `Reconstruction/<run>/ecal_<run>.root` (solo lectura) de los ~51 runs de `beam_points.txt`, salida
   `Simulation/tests/border_cut/data/<run>/ecal_<run>.valtree.root` (~60–80 GB), DAG en
   `analysis/condor_border_data/` (script por escribir).
4. **Plots** (`FINAL_SIM=…/Processed/final_v8`, `FINAL_RECO=…/tests/border_cut/data`,
   `FINAL_RECO_TREES=…/Reconstruction`), salida `Simulation/tests/border_cut/figs_v8/` (en un HTCondor job o local,
   ligero: solo valtrees):
   - sin corte: `final_panels.py <out>/nocut` → `mip_gev_line.txt`, panels y resolución de referencia;
   - `final_panels.py <out>/b50 --border-cut 50 --cal-file <out>/nocut/mip_gev_line.txt` → paneles por grupo (46),
     `resolution_vs_energy.png` (resolución + linealidad), `border_survival.txt`;
   - `final_panels_select.py` (20 y 74 GeV) y `final_profiles.py` (th210 20 y 74 GeV) con `--border-cut 50 --cal-file`.
   - Control: `FINAL_SIM=final_v7` sin corte debe reproducir los números del artefacto v3.7
     (a = 81.9 / 87.2 / 88.7 MIP/GeV).
   - Se esperan pérdidas grandes en el scan th210 (y ≈ −50, mediana ~37–39 mm) y en 99 GeV en (−40, 52) y
     (28, −30); marcar grupos con < 1000 sucesos de data.
5. **Artefacto**: leer entero (580 líneas; leídas 1–403), añadir la sección
   "Lateral containment: barxy_min_border_distance > 50 mm (final_v8)" (definición, final_v8, recta MIP→GeV nueva,
   tabla de supervivencia, las cuatro familias de figuras) y una entrada en el historial de versiones, en inglés;
   republicar en la misma URL.

## Notas de entorno (sesión del 2026-10-28)
- `condor_*` y `cvmfs_config probe` corren fuera del sandbox (`sandbox.excludedCommands`) solo como comando suelto.
- CVMFS/EOS/DUST se mantienen montados con `~/claude-config/hooks/keepalive-mounts.sh` (hook `SessionStart`); si
  `/cvmfs/sw.hsf.org` da "Operation not permitted", `cvmfs_config probe sw.hsf.org`.
- El clasificador de auto mode frena `make -jN` con N > 2 en nodos de login: usar `make -j2`.
- Los chunks de cada campaña se borran al cerrar el grupo (`sim_reco_group.sh`); para pruebas cortas, extraer
  sucesos de `trees/` con `rooteventselector -l N`.
