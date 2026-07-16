# Pipeline commands - signal_analysis

Available steps in this repository (a subset of the original `nogse_pipeline`;
contrast fitting and tc fitting steps are not included here):

| Step | Script | Name for `run_dataset.sh` |
|---|---|---|
| 00 | `bash_template/steps/00_filter_master_points.sh` | `filter_master_points` |
| 01 | `bash_template/steps/01_ingest_results.sh` | `ingest` |
| 02 | `bash_template/steps/02_rotate_signals.sh` | `rotate` |
| 03 | `bash_template/steps/03_plot_signals.sh` | `plot_signal` |
| 05 | `bash_template/steps/05_make_grad_correction.sh` | `grad_correction` |
| 06 | `bash_template/steps/06_alpha_macro.sh` | `alpha` |
| 06b | `bash_template/steps/06b_plot_monoexp_D_vs_time.sh` | `plot_monoexp_d` |
| 99 | `bash_template/steps/99_export_master_xlsx.sh` | `export_master_xlsx` |

All steps are run through `bash_template/run_dataset.sh <type_subj> <type_seq> <step...>`,
which resolves the requested step or steps, applies dataset defaults, and sets
the shared environment variables from `master_table_common.sh`. The pipeline
requires `signal_extraction` to have already generated signal tables under
`Data-BIDS/derivatives/signal_extraction/`.

Run commands from `PROJECT_ROOT`, the directory that contains `Data-BIDS/`,
`Data-signals/`, `analysis/`, and the sibling `signal_analysis/` directory.

By default, `ingest` looks for results under:

```text
Data-BIDS/derivatives/signal_extraction/$DWI_LEVEL/$ROI_VARIANT
```

with `DWI_LEVEL=den_gr-topup` and `ROI_VARIANT=plain` unless overridden. To use
a different `signal_extraction` namespace, pass `DWI_LEVEL` and `ROI_VARIANT`,
or use `--results-root` pointing to a specific `Results` folder, for example:

```bash
DWI_LEVEL=den_gr-topup ROI_VARIANT=sket1 \
  bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest

bash signal_analysis/bash_template/run_dataset.sh brain ogse \
  --results-root Data-BIDS/derivatives/signal_extraction/den_gr-topup/sket1/sub-MBBL-3/ses-T0/Results \
  ingest
```

---

## How Variables Behave Across Steps

When multiple steps are passed in one command, for example
`run_dataset.sh brain ogse ingest rotate ...`, the runner executes them in
sequence within the same parent process. Each step runs in an isolated
subshell (`bash step_script.sh`), so variables set internally by one step do
not propagate to the next step.

Variables that do propagate:

- Any variable set before calling `run_dataset.sh` in the external environment
  is inherited by all steps.
- `MASTER_PARQUET` and `MASTER_LAST_POINTS_APPLIED` are exported in the parent
  process automatically when `MASTER_LAST_POINTS_BY_TD` is used, so the
  filtering is applied across the full step sequence without repeating it.

There are no variable conflicts between steps because each step runs in an
isolated subshell.

---

## Unfiltered - Brains (OGSE)

```bash
mkdir -p logs
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest > logs/01_ingest.log 2>&1 &
ROTATE_EXTRA_ARGS="--solver solve" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse rotate > logs/02_rotate.log 2>&1 &
PLOT_SIGNAL_YCOL=value nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
ALPHA_EXTRA_ARGS="--bvalmax 10 --roi-bvalmax Syringe=7 --roi-bvalmax Right-Lateral-Ventricle=5 --roi-bvalmax Left-Lateral-Ventricle=5" DPROJ_DIRS="long tra" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/06_alpha.log 2>&1 &
```

### All-in-one Core Chain (`ingest -> rotate -> grad_correction`)

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest rotate grad_correction > logs/brain_core.log 2>&1 &
```

`plot_signal`, `alpha`, and `plot_monoexp_d` are usually run separately because
each one often needs its own filters and extra arguments.

---

## Unfiltered - Phantoms (OGSE)

```bash
PARAMS_XLSX=Data-signals/sequence_parameters_phantoms.xlsx nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest > logs/01_ingest.log 2>&1 &
ROTATE_EXTRA_ARGS="--solver solve" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse rotate > logs/02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_signal > logs/03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction > logs/05_grad_correction.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/06_alpha.log 2>&1 &
```

### All-in-one Core Chain

```bash
PARAMS_XLSX=Data-signals/sequence_parameters_phantoms.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest rotate grad_correction > logs/phantom_core.log 2>&1 &
```

---

## Filter the Last N Points by td and N - Brains (OGSE)

Replace `120:8=6,120:4=4,210=8` with the values needed for the dataset. The
format is `td=points` for all N values at that td, or `td:N=points` for a
specific td and N, separated by commas.

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest > logs/01_ingest.log 2>&1 &
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse filter_master_points > logs/00_filter.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse rotate > logs/02_rotate.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/03_plot_signal.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/06_alpha.log 2>&1 &
```

### All-in-one Filtered Core Chain

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest filter_master_points rotate grad_correction > logs/brain_filtered.log 2>&1 &
```

> In all-in-one filtered mode, `MASTER_LAST_POINTS_BY_TD` triggers filtering
> automatically before each step. The pipeline detects when it has already been
> applied and does not repeat the filter.

---

## Filter the Last N Points by td and N - Phantoms (OGSE)

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest > logs/01_ingest.log 2>&1 &
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse filter_master_points > logs/00_filter.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse rotate > logs/02_rotate.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_signal > logs/03_plot_signal.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction > logs/05_grad_correction.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/06_alpha.log 2>&1 &
```

### All-in-one Filtered Core Chain

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest filter_master_points rotate grad_correction > logs/phantom_filtered.log 2>&1 &
```

---

## NOGSE - Brains

Differences from OGSE:

- No `grad_correction` step is provided for NOGSE in this repository. There is
  no `grad_correction.csv` manifest for `brains_nogse` or `phantoms_nogse`, so
  the step would fail because the manifest is missing. If needed, copy and
  complete `bash_template/manifests/phantoms_ogse/grad_correction.csv` as a
  template under a new `brains_nogse/` or `phantoms_nogse/` manifest folder.
- Gradient x-axis: `g` for NOGSE instead of `bvalue_thorsten` or `g_thorsten`
  for OGSE.

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest > logs/brain_nogse_01_ingest.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse rotate > logs/brain_nogse_02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse plot_signal > logs/brain_nogse_03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse alpha > logs/brain_nogse_06_alpha.log 2>&1 &
```

### All-in-one Core Chain

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest rotate > logs/brain_nogse_core.log 2>&1 &
```

### Filter the Last N Points by td and N

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest filter_master_points rotate > logs/brain_nogse_filtered.log 2>&1 &
```

---

## NOGSE - Phantoms

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest > logs/phantom_nogse_01_ingest.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse rotate > logs/phantom_nogse_02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse plot_signal > logs/phantom_nogse_03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse alpha > logs/phantom_nogse_06_alpha.log 2>&1 &
```

### All-in-one Core Chain

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest rotate > logs/phantom_nogse_core.log 2>&1 &
```

### Filter the Last N Points by td and N

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest filter_master_points rotate > logs/phantom_nogse_filtered.log 2>&1 &
```

---

## Export the Master Table as `.xlsx` for Inspection

```bash
MASTER_PARQUET=analysis/brains/ogse_experiments/master.long.parquet MASTER_XLSX=analysis/brains/ogse_experiments/master.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse export_master_xlsx > logs/export_master_xlsx.log 2>&1 &

MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.long.parquet MASTER_XLSX=analysis/phantoms/ogse_experiments/master.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse export_master_xlsx > logs/export_master_xlsx.log 2>&1 &
```

---

## Additional Step Reference

### `plot_signal` - Plot Signal Curves

Generates signal-vs-gradient plots directly from `master.long.parquet`.

| Variable | Default | Description |
|---|---|---|
| `PLOT_OUT_ROOT` | `$ANALYSIS_ROOT/plots-master/signal` | Output directory |
| `PLOT_ROW_KIND` | `signal_rotated` | `signal_rotated` or `signal` |
| `PLOT_SIGNAL_YCOL` | `value_norm` | `value` or `value_norm` |
| `PLOT_SIGNAL_XCOL` | depends on `type_seq` | X-axis column, e.g. `g_thorsten`, `g`, or `bvalue_thorsten` |
| `PLOT_STAT` | `avg` | `avg` or `std` |
| `PLOT_SUBJ` | all | Subject filter |
| `PLOT_ROI` | all | ROI filter |
| `PLOT_DIRECTION` | all | Direction filter, e.g. `long` or `tra` |
| `PLOT_TD_MS` | all | `td_ms` filter |
| `PLOT_N` | all | N filter |
| `PLOT_SIGNAL_EXTRA_ARGS` | none | Extra arguments for `plot_*_signal_vs_g.py` |

```bash
# Plot signal for a specific ROI and direction
PLOT_ROI=Left-Lateral-Ventricle PLOT_DIRECTION=long \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/plot_signal.log 2>&1 &

# Plot normalized signal using g_thorsten as the x-axis
PLOT_SUBJ=20220622_BRAIN PLOT_DIRECTION=long PLOT_SIGNAL_XCOL=g_thorsten \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/plot_signal.log 2>&1 &
```

---

### `plot_monoexp_d` - Plot Monoexponential D vs Time

Requires monoexponential signal fits generated outside this repository, with
their fit parquets available under `SIGNAL_FITS_ROOT`.

| Variable | Default | Description |
|---|---|---|
| `SIGNAL_FITS_ROOT` | `$ANALYSIS_ROOT/fits/<master>/<experiment>_<model>` | Root scanned for signal fits |
| `MONOEXP_D_OUT_DIR` | `$ANALYSIS_ROOT/plots-master/monoexp_D_vs_time` | Output directory |
| `PLOT_MONOEXP_D_EXTRA_ARGS` | none | Extra arguments for `plot_monoexp_D_vs_time.py` |

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/plot_monoexp_d.log 2>&1 &

# Explicit fit root
SIGNAL_FITS_ROOT=analysis/brains/ogse_experiments/fits/master/ogse_value_norm_vs_bvaluethorsten_monoexp \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/plot_monoexp_d.log 2>&1 &
```

---

### `grad_correction` - Compute Gradient Correction Factors

Fits `D0_nogse` and `D0_monoexp` on signal curves from the syringe or another
reference phantom ROI listed in the manifest, then computes:

```text
correction_factor = sqrt(D0_nogse / D0_monoexp_avg)
```

How `D0_monoexp` is averaged:

- `D0_monoexp` is always averaged across all directions in the manifest that
  share the same `(subj, sheet, roi, td_ms, N)`. `D0_nogse` remains row-specific
  by direction and N so direction-specific gradient errors are preserved.
- With `--avg-N` and no values, the fit is additionally averaged across all N
  values, yielding one `D0_monoexp` per `(subj, sheet, roi, td_ms)`.
- With `--avg-N 4 8`, the same averaging is applied using only rows with
  `N=4` and `N=8`, and the result is applied to all N values.

The resulting factor is written back into `master.long.parquet` for all ROIs
that share the same acquisition parameters.

| Variable | Default | Description |
|---|---|---|
| `GRAD_CORR_SCRIPT` | `$REPO_ROOT/scripts/data/make_grad_correction_table.py` | Python script |
| `GRAD_CORR_MANIFEST` | `$MANIFEST_DIR/grad_correction.csv` | CSV with columns `subj`, `sheet`, `roi`, `direction`, `td_ms`, `N` |
| `GRAD_CORR_OUT_DIR` | `$ANALYSIS_ROOT/fits/grad_correction` | Output directory for `.xlsx` and `.csv` files |
| `GRAD_CORR_PLOT_DIR` | `$GRAD_CORR_OUT_DIR/plots` | Output directory for before/after correction plots |
| `GRAD_CORR_ROI` | `Syringe` for brains, `Water1` for phantoms | Reference ROI matched in the master table for all manifest rows |
| `GRAD_CORR_EXTRA_ARGS` | none | Extra arguments for `make_grad_correction_table.py` |

Useful `GRAD_CORR_EXTRA_ARGS`: `--avg-N`, `--avg-N 4 8`,
`--no-fill-missing`, `--row-kind signal_rotated`, `--stat avg`, `--free-M0`.

```bash
# Standard correction: D0_monoexp averaged across directions, one factor per direction x td x N
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

# Also average across all N values: one factor per direction x td
GRAD_CORR_EXTRA_ARGS="--avg-N" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

# Average across directions and only N=1, N=4, and N=8
GRAD_CORR_EXTRA_ARGS="--avg-N 1 4 8" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &

# Do not fill missing factors with averages from other subjects
GRAD_CORR_EXTRA_ARGS="--no-fill-missing" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

# Use another reference ROI. Defaults are Syringe for brains and Water1 for phantoms.
GRAD_CORR_ROI=Water2 \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &
```

---

### `alpha` - Compute Macroscopic Alpha (D vs Delta)

Computes macroscopic alpha from signal fits generated outside this repository.
Also generates D-vs-Delta_app plots.

| Variable | Default | Description |
|---|---|---|
| `ALPHA_N` | `1` | N selector for `make_alpha_macro_summary.py` |
| `ALPHA_OUT_DIR` | `$ANALYSIS_ROOT/alpha_macro/master` | Output directory |
| `ALPHA_EXTRA_ARGS` | none | Extra arguments for `make_alpha_macro_summary.py` |
| `DPROJ_N` | same as `ALPHA_N` | N selector for D-vs-Delta plots when it should differ from `ALPHA_N` |
| `DPROJ_DIRS` | none | Direction filter for D-projection |
| `DPROJ_ROIS` | none | ROI filter for D-projection |
| `PLOT_D0_EXTRA_ARGS` | none | Extra arguments only for `plot_D0_vs_Delta.py` |

Useful `ALPHA_EXTRA_ARGS`: `--bvalmax 5`, `--roi-bvalmax AntCC=7`,
`--dirs long tra`.

Outputs under `$ALPHA_OUT_DIR/`:

- `summary_alpha_values.xlsx` - alpha_macro by subject, ROI, and direction
- `D_vs_delta_app.combined.xlsx` - aggregated D vs Delta_app
- `alpha_macro_vs_roi.png` - summary plot by ROI
- `<subj>/<roi>/<dir>_*.png` - individual curves

```bash
# Standard alpha for brains (N=1)
ALPHA_N=1 ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Left-Lateral-Ventricle=5 --roi-bvalmax Right-Lateral-Ventricle=5 --roi-bvalmax Syringe=7 --dirs long tra" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

# Alpha for phantoms
ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Syringe=7" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/alpha.log 2>&1 &
```

`summary_alpha_values.xlsx` is a typical input for downstream tc fitting steps
outside this repository, for example `TC_METHOD=pseudohuber_fixed_macro`.

---

## Global Dataset Variables

All variables are inherited from `master_table_common.sh` and can be set before
any command:

| Variable | Description |
|---|---|
| `SIGNALS_ROOT` | Root containing sequence-parameter Excel files; default `$PROJECT_ROOT/Data-signals` |
| `DWI_LEVEL` | `signal_extraction` DWI level used by `ingest`; default `den_gr-topup` |
| `ROI_VARIANT` | `signal_extraction` ROI variant used by `ingest`; default `plain` |
| `RESULTS_ROOT` / `--results-root` repeatable | Results folder or folders to ingest |
| `PARAMS_XLSX` | Sequence-parameter Excel file; required for phantom `ingest` |
| `ANALYSIS_ROOT` | Analysis output root; default `$PROJECT_ROOT/analysis/<brains\|phantoms>/<experiment>` |
| `MASTER_PARQUET` | Master parquet to use or generate |
| `MANIFEST_DIR` | Manifest folder for the current dataset: `bash_template/manifests/<type_subj>_<type_seq>/` |
