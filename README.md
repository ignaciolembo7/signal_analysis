# Pipeline commands — signal_analysis

Pasos disponibles en este repo (subconjunto de la pipeline original de
`nogse_pipeline` — los pasos de fitting/contraste/tc **no** están acá):

| Paso | Script | Nombre para `run_dataset.sh` |
|---|---|---|
| 00 | `bash_template/steps/00_filter_master_points.sh` | `filter_master_points` |
| 01 | `bash_template/steps/01_ingest_results.sh` | `ingest` |
| 02 | `bash_template/steps/02_rotate_signals.sh` | `rotate` |
| 03 | `bash_template/steps/03_plot_signals.sh` | `plot_signal` |
| 05 | `bash_template/steps/05_make_grad_correction.sh` | `grad_correction` |
| 06 | `bash_template/steps/06_alpha_macro.sh` | `alpha` |
| 06b | `bash_template/steps/06b_plot_monoexp_D_vs_time.sh` | `plot_monoexp_d` |
| 99 | `bash_template/steps/99_export_master_xlsx.sh` | `export_master_xlsx` |

Todos se corren a través de `bash_template/run_dataset.sh <type_subj> <type_seq> <step...>`,
que resuelve el/los paso(s), aplica los defaults de dataset y setea las
variables del entorno común (`master_table_common.sh`). Requiere que
`signal_extraction` ya haya generado las tablas de señal bajo
`Data-BIDS/derivatives/signal_extraction/`.

Correr los comandos desde `PROJECT_ROOT` (el directorio que contiene
`Data-BIDS/`, `Data-signals/`, `analysis/`, etc., y del cual `signal_analysis/` cuelga
como carpeta hermana).

Por default, `ingest` busca resultados bajo:

```text
Data-BIDS/derivatives/signal_extraction/$DWI_LEVEL/$ROI_VARIANT
```

con `DWI_LEVEL=den_gr-topup` y `ROI_VARIANT=plain` si no se sobreescriben.
Para usar otro namespace de `signal_extraction`, pasá `DWI_LEVEL` y
`ROI_VARIANT`, o usá `--results-root` apuntando a una carpeta `Results`
específica, por ejemplo:

```bash
DWI_LEVEL=den_gr-topup ROI_VARIANT=sket1 \
  bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest

bash signal_analysis/bash_template/run_dataset.sh brain ogse \
  --results-root Data-BIDS/derivatives/signal_extraction/den_gr-topup/sket1/sub-MBBL-3/ses-T0/Results \
  ingest
```

---

## Cómo funcionan las variables entre pasos

Cuando pasás múltiples pasos en un mismo comando (ej.
`run_dataset.sh brain ogse ingest rotate ...`), el runner los ejecuta en
secuencia dentro del mismo proceso. Cada paso corre en un **subshell
aislado** (`bash step_script.sh`), así que variables que un paso setea
internamente NO se propagan al siguiente.

Lo que SÍ se propaga:
- Cualquier variable que seteás **antes** de llamar a `run_dataset.sh` (en
  el entorno externo) la heredan TODOS los pasos.
- `MASTER_PARQUET` y `MASTER_LAST_POINTS_APPLIED` se exportan en el proceso
  padre automáticamente cuando usás `MASTER_LAST_POINTS_BY_TD`, así que el
  filtrado se aplica a todos los pasos en secuencia sin repetir el filtro.

Conflictos entre pasos: no hay, porque cada paso corre en subshell aislado.

---

## Sin filtrar — brains (OGSE)

```
mkdir -p logs
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest > logs/01_ingest.log 2>&1 &
ROTATE_EXTRA_ARGS="--solver solve" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse rotate > logs/02_rotate.log 2>&1 &
PLOT_SIGNAL_YCOL=value nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
ALPHA_EXTRA_ARGS="--bvalmax 10 --roi-bvalmax Syringe=7 --roi-bvalmax Right-Lateral-Ventricle=5 --roi-bvalmax Left-Lateral-Ventricle=5" DPROJ_DIRS="long tra" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/06_alpha.log 2>&1 &
```

### Todo en uno (cadena `ingest → rotate → grad_correction`)

```
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest rotate grad_correction > logs/brain_core.log 2>&1 &
```

`plot_signal`, `alpha` y `plot_monoexp_d` normalmente se corren aparte (cada
uno con sus propios filtros/args), como en el bloque de arriba.

---

## Sin filtrar — phantoms (OGSE)

```
PARAMS_XLSX=Data-signals/sequence_parameters_phantoms.xlsx nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest > logs/01_ingest.log 2>&1 &
ROTATE_EXTRA_ARGS="--solver solve" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse rotate > logs/02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_signal > logs/03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction > logs/05_grad_correction.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/06_alpha.log 2>&1 &
```

### Todo en uno

```
PARAMS_XLSX=Data-signals/sequence_parameters_phantoms.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest rotate grad_correction > logs/phantom_core.log 2>&1 &
```

---

## Con filtrar últimos N puntos por td y N — brains (OGSE)

Reemplazá `120:8=6,120:4=4,210=8` con tus valores. Formato: `td=puntos`
(todos los N de ese td) o `td:N=puntos` (td y N específicos), separados por
coma.

```
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest > logs/01_ingest.log 2>&1 &
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse filter_master_points > logs/00_filter.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse rotate > logs/02_rotate.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/03_plot_signal.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/06_alpha.log 2>&1 &
```

### Todo en uno

```
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest filter_master_points rotate grad_correction > logs/brain_filtered.log 2>&1 &
```

> En el modo "todo en uno" con filtrado, `MASTER_LAST_POINTS_BY_TD` dispara
> el filtrado automáticamente antes de cada paso (el pipeline lo detecta y
> aplica la primera vez, luego no lo repite).

---

## Con filtrar últimos N puntos por td y N — phantoms (OGSE)

```
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest > logs/01_ingest.log 2>&1 &
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse filter_master_points > logs/00_filter.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse rotate > logs/02_rotate.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_signal > logs/03_plot_signal.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction > logs/05_grad_correction.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/06_alpha.log 2>&1 &
```

### Todo en uno

```
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest filter_master_points rotate grad_correction > logs/phantom_filtered.log 2>&1 &
```

---

## NOGSE — brains

Diferencias respecto a OGSE:
- **Sin `grad_correction`**: no hay manifest `grad_correction.csv` para
  `brains_nogse`/`phantoms_nogse` en este repo — el paso fallaría por
  manifest faltante. Si lo necesitás, copiá y completá
  `bash_template/manifests/phantoms_ogse/grad_correction.csv` como
  plantilla dentro de una carpeta `brains_nogse/`/`phantoms_nogse/` nueva.
- Eje X de gradiente: `g` (NOGSE) en vez de `bvalue_thorsten`/`g_thorsten` (OGSE).

```
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest > logs/brain_nogse_01_ingest.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse rotate > logs/brain_nogse_02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse plot_signal > logs/brain_nogse_03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse alpha > logs/brain_nogse_06_alpha.log 2>&1 &
```

### Todo en uno

```
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest rotate > logs/brain_nogse_core.log 2>&1 &
```

### Con filtrar últimos N puntos por td y N

```
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest filter_master_points rotate > logs/brain_nogse_filtered.log 2>&1 &
```

---

## NOGSE — phantoms

```
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest > logs/phantom_nogse_01_ingest.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse rotate > logs/phantom_nogse_02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse plot_signal > logs/phantom_nogse_03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse alpha > logs/phantom_nogse_06_alpha.log 2>&1 &
```

### Todo en uno

```
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest rotate > logs/phantom_nogse_core.log 2>&1 &
```

### Con filtrar últimos N puntos por td y N

```
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest filter_master_points rotate > logs/phantom_nogse_filtered.log 2>&1 &
```

---

## Exportar tabla master como `.xlsx` (visualización)

```
MASTER_PARQUET=analysis/brains/ogse_experiments/master.long.parquet MASTER_XLSX=analysis/brains/ogse_experiments/master.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse export_master_xlsx > logs/export_master_xlsx.log 2>&1 &

MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.long.parquet MASTER_XLSX=analysis/phantoms/ogse_experiments/master.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse export_master_xlsx > logs/export_master_xlsx.log 2>&1 &
```

---

## Referencia de pasos adicionales

### `plot_signal` — graficar curvas de señal

Genera plots de señal vs gradiente directamente del `master.long.parquet`.

| Variable | Default | Descripción |
|---|---|---|
| `PLOT_OUT_ROOT` | `$ANALYSIS_ROOT/plots-master/signal` | Directorio de salida |
| `PLOT_ROW_KIND` | `signal_rotated` | `signal_rotated` o `signal` |
| `PLOT_SIGNAL_YCOL` | `value_norm` | `value` o `value_norm` |
| `PLOT_SIGNAL_XCOL` | (según tipo_seq) | Columna del eje X (ej. `g_thorsten`, `g`, `bvalue_thorsten`) |
| `PLOT_STAT` | `avg` | `avg` o `std` |
| `PLOT_SUBJ` | (todos) | Filtro de sujeto |
| `PLOT_ROI` | (todos) | Filtro de ROI |
| `PLOT_DIRECTION` | (todas) | Filtro de dirección (ej. `long`, `tra`) |
| `PLOT_TD_MS` | (todos) | Filtro de td_ms |
| `PLOT_N` | (todos) | Filtro de N |
| `PLOT_SIGNAL_EXTRA_ARGS` | — | Args extra para `plot_*_signal_vs_g.py` |

```
# Plot señal para un ROI y dirección específica
PLOT_ROI=Left-Lateral-Ventricle PLOT_DIRECTION=long \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/plot_signal.log 2>&1 &

# Plot normalizado con g_thorsten como eje X
PLOT_SUBJ=20220622_BRAIN PLOT_DIRECTION=long PLOT_SIGNAL_XCOL=g_thorsten \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/plot_signal.log 2>&1 &
```

---

### `plot_monoexp_d` — graficar D monoexp vs tiempo

Requiere haber corrido un fit de señal monoexponencial previamente (fuera
de este repo) y que sus parquets de fits estén accesibles bajo
`SIGNAL_FITS_ROOT`.

| Variable | Default | Descripción |
|---|---|---|
| `SIGNAL_FITS_ROOT` | `$ANALYSIS_ROOT/fits/<master>/<experiment>_<model>` | Root escaneado para fits de señal |
| `MONOEXP_D_OUT_DIR` | `$ANALYSIS_ROOT/plots-master/monoexp_D_vs_time` | Directorio de salida |
| `PLOT_MONOEXP_D_EXTRA_ARGS` | — | Args extra para `plot_monoexp_D_vs_time.py` |

```
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/plot_monoexp_d.log 2>&1 &

# Con root explícito
SIGNAL_FITS_ROOT=analysis/brains/ogse_experiments/fits/master/ogse_value_norm_vs_bvaluethorsten_monoexp \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/plot_monoexp_d.log 2>&1 &
```

---

### `grad_correction` — calcular factores de corrección de gradiente

Ajusta D0_nogse y D0_monoexp sobre las curvas de señal de la jeringilla (u
otro phantom de referencia) listadas en el manifest, y calcula
`correction_factor = sqrt(D0_nogse / D0_monoexp_avg)`.

**Cómo se promedia D0_monoexp:**

- D0_monoexp **siempre** se promedia entre todas las direcciones del
  manifest que comparten el mismo (subj, sheet, roi, td_ms, N). El D0_nogse
  se mantiene por fila (direction- y N-específico) para capturar el error
  de gradiente por dirección.
- Con `--avg-N` (sin valores): además se promedia entre todos los N → un
  solo D0_monoexp por (subj, sheet, roi, td_ms).
- Con `--avg-N 4 8`: igual, pero el promedio usa solo las filas con N=4 y
  N=8 → el valor se aplica a todos los N.

El factor resultante se escribe en el `master.long.parquet` para todas las
rois que compartan los mismos parámetros de adquisición.

| Variable | Default | Descripción |
|---|---|---|
| `GRAD_CORR_SCRIPT` | `$REPO_ROOT/scripts/data/make_grad_correction_table.py` | Script Python |
| `GRAD_CORR_MANIFEST` | `$MANIFEST_DIR/grad_correction.csv` | CSV con columnas subj, sheet, roi, direction, td_ms, N |
| `GRAD_CORR_OUT_DIR` | `$ANALYSIS_ROOT/fits/grad_correction` | Directorio de salida (.xlsx, .csv) |
| `GRAD_CORR_PLOT_DIR` | `$GRAD_CORR_OUT_DIR/plots` | Directorio de plots comparativos (antes/después de corrección) |
| `GRAD_CORR_ROI` | `Syringe` (brains) / `Water1` (phantoms) | ROI de referencia a matchear en el master table para todas las filas del manifest |
| `GRAD_CORR_EXTRA_ARGS` | — | Args extra para `make_grad_correction_table.py` |

Args útiles en `GRAD_CORR_EXTRA_ARGS`: `--avg-N`, `--avg-N 4 8`,
`--no-fill-missing`, `--row-kind signal_rotated`, `--stat avg`, `--free-M0`.

```bash
# Corrección estándar (D0_monoexp promediado entre dirs, un factor por dir×td×N)
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

# Promediar también entre todos los N → un factor por dir×td
GRAD_CORR_EXTRA_ARGS="--avg-N" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

# Promediar entre dirs y entre N=4 y N=8 solamente
GRAD_CORR_EXTRA_ARGS="--avg-N 1 4 8" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &

# Sin rellenar factores faltantes con promedio entre sujetos
GRAD_CORR_EXTRA_ARGS="--no-fill-missing" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

# Usar otro ROI de referencia (por default: Syringe en brains, Water1 en phantoms)
GRAD_CORR_ROI=Water2 \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &
```

---

### `alpha` — calcular alpha macro (D vs Delta)

Calcula alpha macroscópico a partir de fits de señal generados fuera de
este repo. También genera plots de D vs Delta_app.

| Variable | Default | Descripción |
|---|---|---|
| `ALPHA_N` | `1` | Selector de N para `make_alpha_macro_summary.py` |
| `ALPHA_OUT_DIR` | `$ANALYSIS_ROOT/alpha_macro/master` | Directorio de salida |
| `ALPHA_EXTRA_ARGS` | — | Args extra para `make_alpha_macro_summary.py` |
| `DPROJ_N` | (igual a ALPHA_N) | Selector de N para plots D vs Delta (si difiere) |
| `DPROJ_DIRS` | — | Filtro de direcciones para D-projection |
| `DPROJ_ROIS` | — | Filtro de ROIs para D-projection |
| `PLOT_D0_EXTRA_ARGS` | — | Args extra solo para `plot_D0_vs_Delta.py` |

Args útiles en `ALPHA_EXTRA_ARGS`: `--bvalmax 5`, `--roi-bvalmax AntCC=7`,
`--dirs long tra`.

Salidas en `$ALPHA_OUT_DIR/`:
- `summary_alpha_values.xlsx` — alpha_macro por sujeto/ROI/dirección
- `D_vs_delta_app.combined.xlsx` — D vs Delta_app agregado
- `alpha_macro_vs_roi.png` — plot resumen por ROI
- `<subj>/<roi>/<dir>_*.png` — curvas individuales

```
# Alpha estándar para brains (N=1)
ALPHA_N=1 ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Left-Lateral-Ventricle=5 --roi-bvalmax Right-Lateral-Ventricle=5 --roi-bvalmax Syringe=7 --dirs long tra" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

# Alpha para phantoms
ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Syringe=7" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/alpha.log 2>&1 &
```

`summary_alpha_values.xlsx` es un insumo típico para pasos de fitting de tc
(fuera de este repo, ej. `TC_METHOD=pseudohuber_fixed_macro`).

---

## Variables globales de dataset

Todas heredadas de `master_table_common.sh`, seteables antes de cualquier
comando:

| Variable | Descripción |
|---|---|
| `SIGNALS_ROOT` | Raíz con Excels de parámetros de secuencia (default `$PROJECT_ROOT/Data-signals`) |
| `DWI_LEVEL` | Nivel de DWI de `signal_extraction` para `ingest` (default `den_gr-topup`) |
| `ROI_VARIANT` | Variante de ROI de `signal_extraction` para `ingest` (default `plain`) |
| `RESULTS_ROOT` / `--results-root` (repetible) | Carpeta(s) de resultados a ingestar |
| `PARAMS_XLSX` | Excel de parámetros de secuencia (obligatorio para `ingest` en phantoms) |
| `ANALYSIS_ROOT` | Raíz de salida de análisis (default `$PROJECT_ROOT/analysis/<brains\|phantoms>/<experiment>`) |
| `MASTER_PARQUET` | Parquet master a usar/generar |
| `MANIFEST_DIR` | Carpeta de manifests para el dataset actual (`bash_template/manifests/<tipo_subj>_<tipo_seq>/`) |
