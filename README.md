# signal_analysis

Master-table ingestion, signal rotation, plotting, gradient correction, and
summary analysis for `signal_extraction` outputs.

This repository consumes canonical pipeline outputs from:

- `Data-BIDS/derivatives/signal_extraction/<DWI_LEVEL>/<ROI_VARIANT>/sub-*/ses-*/Results/`
- sequence-parameter workbooks under `Data-signals/`
- step manifests under `signal_analysis/bash_template/manifests/`

It writes analysis products under:

```text
analysis/<brains|phantoms>/<ogse_experiments|nogse_experiments>/
```

The canonical table is:

```text
analysis/<brains|phantoms>/<experiment>/master.long.parquet
```

Excel files are exported only for inspection. The pipeline appends, filters,
rotates, and updates the parquet master table.

## Repository Layout

```text
signal_analysis/
├── README.md
├── requirements.txt
├── environment.yml
├── pyproject.toml
├── assets/
│   └── dirs/                         rotation direction tables
├── bash_template/
│   ├── helpers/master_table_common.sh
│   ├── manifests/
│   │   ├── brains_ogse/
│   │   │   └── grad_correction.csv
│   │   └── phantoms_ogse/
│   │       ├── grad_correction.csv
│   │       └── grad_correction-PHANTOM3.csv
│   ├── run_dataset.sh
│   └── steps/
│       ├── 00_filter_master_points.sh
│       ├── 01_ingest_results.sh
│       ├── 02_rotate_signals.sh
│       ├── 03_plot_signals.sh
│       ├── 05_make_grad_correction.sh
│       ├── 06_alpha_macro.sh
│       ├── 06b_plot_monoexp_D_vs_time.sh
│       └── 99_export_master_xlsx.sh
├── scripts/
│   ├── data/
│   ├── plotting/
│   └── summary/
├── src/
│   ├── data_processing/
│   ├── fitting/
│   ├── monoexp_fitting/
│   ├── ogse_fitting/
│   ├── pipeline/
│   ├── plotting/
│   ├── signal_rotation/
│   ├── tc_fittings/
│   └── tools/
└── tests/
```

Run commands from `PROJECT_ROOT`, the directory that contains `Data-BIDS/`,
`Data-signals/`, `analysis/`, and the sibling pipeline repositories:

```text
PROJECT_ROOT/
├── Data-BIDS/
├── Data-signals/
├── analysis/
├── dicom_to_bids/
├── img_preprocessing/
├── img_segmentation/
├── signal_extraction/
└── signal_analysis/
```

## Installation

Clone the repository inside the same `PROJECT_ROOT` that contains `Data-BIDS/`:

```bash
cd /path/to/PROJECT_ROOT
git clone git@github.com:YOUR_USERNAME/signal_analysis.git
cd signal_analysis
```

If you use HTTPS instead of SSH:

```bash
cd /path/to/PROJECT_ROOT
git clone https://github.com/YOUR_USERNAME/signal_analysis.git
cd signal_analysis
```

Create a Python environment and install the repository:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
python -m pip install -r requirements.txt
```

If you prefer conda:

```bash
conda env create -f environment.yml
conda activate signal_analysis_env
python -m pip install -e .
```

Prepare the local log directory:

```bash
cd /path/to/PROJECT_ROOT
mkdir -p logs
```

Make sure scripts are executable:

```bash
chmod +x signal_analysis/bash_template/run_dataset.sh \
  signal_analysis/bash_template/helpers/*.sh \
  signal_analysis/bash_template/steps/*.sh
```

## Pipeline Inputs

### Signal Extraction Results

`ingest` reads result workbooks from the canonical `signal_extraction` layout:

```text
Data-BIDS/derivatives/signal_extraction/<DWI_LEVEL>/<ROI_VARIANT>/sub-<SUBJECT_ID>/ses-<SESSION_ID>/Results/*_results.xlsx
```

By default:

```text
DWI_LEVEL=den_gr-topup
ROI_VARIANT=plain
RESULTS_ROOT=Data-BIDS/derivatives/signal_extraction/$DWI_LEVEL/$ROI_VARIANT
```

Examples:

```bash
DWI_LEVEL=den_gr-topup ROI_VARIANT=plain \
  bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest

DWI_LEVEL=den_gr-topup ROI_VARIANT=sket1 \
  bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest

DWI_LEVEL=den_gr-topup ROI_VARIANT=sket2 \
  bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest
```

To ingest one explicit `Results` folder:

```bash
bash signal_analysis/bash_template/run_dataset.sh brain ogse \
  --results-root Data-BIDS/derivatives/signal_extraction/den_gr-topup/sket1/sub-MBBL-3/ses-T0/Results \
  ingest
```

The old nested layout `Results/<subject>` is not part of the pipeline contract.
Each `Results` folder must contain `*_results.xlsx` directly.

### Sequence Parameters

The ingest step matches result rows against sequence-parameter Excel workbooks:

| Dataset | Default workbook |
|---|---|
| `brain` | `Data-signals/sequence_parameters_brains.xlsx` |
| `phantom` | `Data-signals/sequence_parameters_phantoms.xlsx` |

Override with:

```bash
PARAMS_XLSX=/path/to/sequence_parameters.xlsx \
  bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest
```

### Manifests

Some downstream steps need explicit manifest CSV files. The runner selects the
manifest directory from:

```text
signal_analysis/bash_template/manifests/<brains|phantoms>_<ogse|nogse>/
```

Current tracked manifests:

```text
bash_template/manifests/brains_ogse/grad_correction.csv
bash_template/manifests/phantoms_ogse/grad_correction.csv
bash_template/manifests/phantoms_ogse/grad_correction-PHANTOM3.csv
```

`grad_correction.csv` columns:

| Column | Meaning |
|---|---|
| `subj` | Subject/group key used to match master-table rows |
| `sheet` | Sequence-parameter workbook sheet |
| `roi` | Reference ROI used to estimate the correction factor |
| `direction` | Direction selector, e.g. `long` or `tra` |
| `td_ms` | Diffusion time in ms |
| `N` | Oscillation count |
| `Hz` | Frequency label used by the acquisition |
| `model` | Fit model, usually `monoexp` |

NOGSE has no tracked `grad_correction.csv` manifest in this repository. Running
`grad_correction` for NOGSE will fail unless you create the corresponding
manifest folder and CSV.

Override the manifest directory with:

```bash
MANIFEST_DIR=signal_analysis/bash_template/manifests/phantoms_ogse \
  bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction
```

## Outputs

Default output roots:

| Dataset | Sequence | Default `ANALYSIS_ROOT` |
|---|---|---|
| `brain` | `ogse` | `analysis/brains/ogse_experiments` |
| `brain` | `nogse` | `analysis/brains/nogse_experiments` |
| `phantom` | `ogse` | `analysis/phantoms/ogse_experiments` |
| `phantom` | `nogse` | `analysis/phantoms/nogse_experiments` |

Typical output layout:

```text
analysis/brains/ogse_experiments/
├── master.long.parquet
├── master.last_points.long.parquet
├── master.xlsx
├── data/
│   └── tables/
├── rotated/
├── plots-master/
│   ├── signal/
│   └── monoexp_D_vs_time/
├── fits/
│   └── grad_correction/
└── alpha_macro/
```

Important outputs:

| Output | Created by | Meaning |
|---|---|---|
| `master.long.parquet` | `ingest` | Canonical long-format master table |
| `master.last_points.long.parquet` | `filter_master_points` | Filtered master table |
| rotated signal rows | `rotate` | Tensor-rotated signal directions appended to the master table |
| `fits/grad_correction/*` | `grad_correction` | Correction tables, plots, and updated master factors |
| `plots-master/signal/*` | `plot_signal` | Signal-vs-gradient figures |
| `plots-master/monoexp_D_vs_time/*` | `plot_monoexp_d` | Monoexponential D-vs-time figures |
| `alpha_macro/*` | `alpha` | Macroscopic alpha tables and plots |
| `master.xlsx` | `export_master_xlsx` | Inspection-only Excel export |

## Runner Interface

All steps are run through:

```bash
bash signal_analysis/bash_template/run_dataset.sh <type_subj> <type_seq> <step...>
```

Valid subject types:

```text
brain, brains, phantom, phantoms
```

Valid sequence types:

```text
ogse, nogse
```

Equivalent option form:

```bash
bash signal_analysis/bash_template/run_dataset.sh --type-subj brain --type-seq ogse ingest
```

Get help:

```bash
bash signal_analysis/bash_template/run_dataset.sh --help
bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest --help
bash signal_analysis/bash_template/run_dataset.sh brain ogse rotate --help
```

Step-specific settings are environment variables placed before the command:

```bash
VAR=value OTHER_VAR=value \
  bash signal_analysis/bash_template/run_dataset.sh brain ogse <step>
```

Extra Python flags go through the corresponding `*_EXTRA_ARGS` variable:

```bash
GRAD_CORR_EXTRA_ARGS="--avg-N --no-fill-missing" \
  bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction
```

## Available Steps

| Step name | Purpose |
|---|---|
| `ingest` | Read `Results/*_results.xlsx` into `master.long.parquet` |
| `filter_master_points` | Create a filtered master table using last-points rules |
| `rotate` | Rotate signal tensor directions |
| `plot_signal` | Plot signal curves from the master table |
| `grad_correction` | Build and embed gradient-correction factors |
| `alpha` | Build alpha_macro summaries and D-vs-Delta plots |
| `plot_monoexp_d` | Plot monoexponential D vs `td_ms` or `Delta_app_ms` |
| `export_master_xlsx` | Export the selected master parquet to Excel |

The runner also exposes some fitting and contrast steps from the shared code
base. Use `--help` to see the full current list.

## How Variables Behave Across Steps

When multiple steps are passed in one command, for example:

```bash
bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest rotate grad_correction
```

the runner executes them in sequence within the same parent process. Each step
runs in an isolated subshell (`bash step_script.sh`), so variables set
internally by one step do not propagate to the next step.

Variables that do propagate:

- Any variable set before calling `run_dataset.sh` in the external environment
  is inherited by all steps.
- `MASTER_PARQUET` and `MASTER_LAST_POINTS_APPLIED` are exported in the parent
  process automatically when `MASTER_LAST_POINTS_BY_TD` is used.

There are no variable conflicts between steps because each step runs in an
isolated subshell.

## Quick Start

Brain OGSE, plain CC ROI, topup-corrected DWI:

```bash
mkdir -p logs
DWI_LEVEL=den_gr-topup ROI_VARIANT=plain \
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest rotate grad_correction \
  > logs/brain_core.log 2>&1 &
```

Brain OGSE, skeleton ROI:

```bash
DWI_LEVEL=den_gr-topup ROI_VARIANT=sket2 \
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest rotate grad_correction \
  > logs/brain_core_sket2.log 2>&1 &
```

Phantom OGSE:

```bash
PARAMS_XLSX=Data-signals/sequence_parameters_phantoms.xlsx \
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest rotate grad_correction \
  > logs/phantom_core.log 2>&1 &
```

## Brain OGSE Commands

Unfiltered run:

```bash
mkdir -p logs
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest > logs/01_ingest.log 2>&1 &
ROTATE_EXTRA_ARGS="--solver solve" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse rotate > logs/02_rotate.log 2>&1 &
PLOT_SIGNAL_YCOL=value nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
ALPHA_EXTRA_ARGS="--bvalmax 10 --roi-bvalmax Syringe=7 --roi-bvalmax Right-Lateral-Ventricle=5 --roi-bvalmax Left-Lateral-Ventricle=5" DPROJ_DIRS="long tra" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/06_alpha.log 2>&1 &
```

All-in-one core chain:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest rotate grad_correction > logs/brain_core.log 2>&1 &
```

`plot_signal`, `alpha`, and `plot_monoexp_d` are usually run separately because
each one often needs its own filters and extra arguments.

## Phantom OGSE Commands

Unfiltered run:

```bash
PARAMS_XLSX=Data-signals/sequence_parameters_phantoms.xlsx nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest > logs/01_ingest.log 2>&1 &
ROTATE_EXTRA_ARGS="--solver solve" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse rotate > logs/02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_signal > logs/03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction > logs/05_grad_correction.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/06_alpha.log 2>&1 &
```

All-in-one core chain:

```bash
PARAMS_XLSX=Data-signals/sequence_parameters_phantoms.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest rotate grad_correction > logs/phantom_core.log 2>&1 &
```

## Filtered Runs

`MASTER_LAST_POINTS_BY_TD` keeps the last `POINTS` b-step values within each
`td_ms` or `(td_ms, N)` group. The format is:

```text
td=points
td:N=points
```

Multiple rules are comma-separated:

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8"
```

Brain OGSE:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest > logs/01_ingest.log 2>&1 &
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse filter_master_points > logs/00_filter.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse rotate > logs/02_rotate.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/03_plot_signal.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/06_alpha.log 2>&1 &
```

Brain OGSE all-in-one filtered core chain:

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse ingest filter_master_points rotate grad_correction > logs/brain_filtered.log 2>&1 &
```

Phantom OGSE:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest > logs/01_ingest.log 2>&1 &
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse filter_master_points > logs/00_filter.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse rotate > logs/02_rotate.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_signal > logs/03_plot_signal.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction > logs/05_grad_correction.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/06_alpha.log 2>&1 &
```

Phantom OGSE all-in-one filtered core chain:

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse ingest filter_master_points rotate grad_correction > logs/phantom_filtered.log 2>&1 &
```

In all-in-one filtered mode, `MASTER_LAST_POINTS_BY_TD` triggers filtering
automatically before downstream steps. The pipeline detects when it has already
been applied and does not repeat the filter.

## NOGSE Commands

NOGSE differs from OGSE in two important ways:

- No `grad_correction` manifest is tracked for NOGSE in this repository.
- The gradient x-axis is `g`, not `bvalue_thorsten` or `g_thorsten`.

Brain NOGSE:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest > logs/brain_nogse_01_ingest.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse rotate > logs/brain_nogse_02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse plot_signal > logs/brain_nogse_03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse alpha > logs/brain_nogse_06_alpha.log 2>&1 &
```

Brain NOGSE all-in-one core chain:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest rotate > logs/brain_nogse_core.log 2>&1 &
```

Brain NOGSE filtered core chain:

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain nogse ingest filter_master_points rotate > logs/brain_nogse_filtered.log 2>&1 &
```

Phantom NOGSE:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest > logs/phantom_nogse_01_ingest.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse rotate > logs/phantom_nogse_02_rotate.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse plot_signal > logs/phantom_nogse_03_plot_signal.log 2>&1 &
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse alpha > logs/phantom_nogse_06_alpha.log 2>&1 &
```

Phantom NOGSE all-in-one core chain:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest rotate > logs/phantom_nogse_core.log 2>&1 &
```

Phantom NOGSE filtered core chain:

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom nogse ingest filter_master_points rotate > logs/phantom_nogse_filtered.log 2>&1 &
```

## Step Reference

### `plot_signal`

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

Examples:

```bash
PLOT_ROI=Left-Lateral-Ventricle PLOT_DIRECTION=long \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/plot_signal.log 2>&1 &

PLOT_SUBJ=20220622_BRAIN PLOT_DIRECTION=long PLOT_SIGNAL_XCOL=g_thorsten \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_signal > logs/plot_signal.log 2>&1 &
```

### `plot_monoexp_d`

Requires monoexponential signal fits generated outside this repository, with
their fit parquets available under `SIGNAL_FITS_ROOT`.

| Variable | Default | Description |
|---|---|---|
| `SIGNAL_FITS_ROOT` | `$ANALYSIS_ROOT/fits/<master>/<experiment>_<model>` | Root scanned for signal fits |
| `MONOEXP_D_OUT_DIR` | `$ANALYSIS_ROOT/plots-master/monoexp_D_vs_time` | Output directory |
| `PLOT_MONOEXP_D_EXTRA_ARGS` | none | Extra arguments for `plot_monoexp_D_vs_time.py` |

Examples:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/plot_monoexp_d.log 2>&1 &

SIGNAL_FITS_ROOT=analysis/brains/ogse_experiments/fits/master/ogse_value_norm_vs_bvaluethorsten_monoexp \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse plot_monoexp_d > logs/plot_monoexp_d.log 2>&1 &
```

### `grad_correction`

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
| `GRAD_CORR_MANIFEST` | `$MANIFEST_DIR/grad_correction.csv` | Manifest CSV |
| `GRAD_CORR_OUT_DIR` | `$ANALYSIS_ROOT/fits/grad_correction` | Output directory for `.xlsx` and `.csv` files |
| `GRAD_CORR_PLOT_DIR` | `$GRAD_CORR_OUT_DIR/plots` | Output directory for before/after correction plots |
| `GRAD_CORR_ROI` | `Syringe` for brains, `Water1` for phantoms | Reference ROI matched in the master table |
| `GRAD_CORR_EXTRA_ARGS` | none | Extra arguments for `make_grad_correction_table.py` |

Useful `GRAD_CORR_EXTRA_ARGS`: `--avg-N`, `--avg-N 4 8`,
`--no-fill-missing`, `--row-kind signal_rotated`, `--stat avg`, `--free-M0`.

Examples:

```bash
nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_EXTRA_ARGS="--avg-N" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_EXTRA_ARGS="--avg-N 1 4 8" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_EXTRA_ARGS="--no-fill-missing" \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_ROI=Water2 \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &
```

### `alpha`

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

Examples:

```bash
ALPHA_N=1 ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Left-Lateral-Ventricle=5 --roi-bvalmax Right-Lateral-Ventricle=5 --roi-bvalmax Syringe=7 --dirs long tra" nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Syringe=7" nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse alpha > logs/alpha.log 2>&1 &
```

`summary_alpha_values.xlsx` is a typical input for downstream tc fitting steps
outside this repository, for example `TC_METHOD=pseudohuber_fixed_macro`.

## Export Master Table to Excel

```bash
MASTER_PARQUET=analysis/brains/ogse_experiments/master.long.parquet MASTER_XLSX=analysis/brains/ogse_experiments/master.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh brain ogse export_master_xlsx > logs/export_master_xlsx.log 2>&1 &

MASTER_PARQUET=analysis/phantoms/ogse_experiments/master.long.parquet MASTER_XLSX=analysis/phantoms/ogse_experiments/master.xlsx \
  nohup bash signal_analysis/bash_template/run_dataset.sh phantom ogse export_master_xlsx > logs/export_master_xlsx.log 2>&1 &
```

## Global Dataset Variables

All variables are inherited from `master_table_common.sh` and can be set before
any command:

| Variable | Description |
|---|---|
| `PY` | Python interpreter |
| `SIGNALS_ROOT` | Root containing sequence-parameter Excel files; default `$PROJECT_ROOT/Data-signals` |
| `DWI_LEVEL` | `signal_extraction` DWI level used by `ingest`; default `den_gr-topup` |
| `ROI_VARIANT` | `signal_extraction` ROI variant used by `ingest`; default `plain` |
| `RESULTS_ROOT` / `--results-root` repeatable | Results folder or folders to ingest |
| `PARAMS_XLSX` | Sequence-parameter Excel file; required for phantom `ingest` |
| `ANALYSIS_ROOT` | Analysis output root; default `$PROJECT_ROOT/analysis/<brains\|phantoms>/<experiment>` |
| `MASTER_PARQUET` | Master parquet to use or generate |
| `MANIFEST_DIR` | Manifest folder for the current dataset |
| `MASTER_LAST_POINTS_BY_TD` | Optional td/N last-points filtering rules |
| `MASTER_LAST_POINTS_PARQUET` | Optional output path for the filtered master table |

## Validation

Run the test suite from `PROJECT_ROOT`:

```bash
python -m unittest discover -s signal_analysis/tests
```

Syntax-check the shell entry points:

```bash
bash -n signal_analysis/bash_template/run_dataset.sh
bash -n signal_analysis/bash_template/steps/01_ingest_results.sh
bash -n signal_analysis/bash_template/steps/02_rotate_signals.sh
bash -n signal_analysis/bash_template/steps/05_make_grad_correction.sh
```
