# signal_analysis

Master-table ingestion, signal rotation, plotting, gradient correction, and
summary analysis for `signal_extraction` outputs.

This repository consumes canonical pipeline outputs from:

- `Data-BIDS/derivatives/signal_extraction/<SIGNAL_EXTRACTION_TAG>/sub-*/ses-*/Results/`
- sequence-parameter workbooks under `Data-BIDS/`
- step manifests under `signal_analysis/manifests/`

It writes analysis products under:

```text
analysis/<SIGNAL_EXTRACTION_TAG>--<brains|phantoms>/<ogse_experiments|nogse_experiments>/
```

The canonical table is:

```text
analysis/<SIGNAL_EXTRACTION_TAG>--<brains|phantoms>/<experiment>/master.long.parquet
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
├── helpers/
│   └── master_table_common.sh
├── manifests/
│   ├── brains_ogse/
│   │   └── grad_correction.csv
│   └── phantoms_ogse/
│       ├── grad_correction.csv
│       └── grad_correction-PHANTOM3.csv
├── run_dataset.sh
├── steps/
│   ├── 00_filter_master_points.sh
│   ├── 01_ingest_results.sh
│   ├── 02_rotate_signals.sh
│   ├── 03_plot_signals.sh
│   ├── 04_fit_signals.sh
│   ├── 05_make_grad_correction.sh
│   ├── 06_alpha_macro.sh
│   ├── 06b_plot_monoexp_D_vs_time.sh
│   └── 99_export_master_xlsx.sh
├── scripts/
│   ├── data/
│   ├── fitting/
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
`analysis/`, and the sibling pipeline repositories:

```text
PROJECT_ROOT/
├── Data-BIDS/
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
chmod +x signal_analysis/run_dataset.sh \
  signal_analysis/helpers/*.sh \
  signal_analysis/steps/*.sh
```

## Pipeline Inputs

### Signal Extraction Results

`ingest` reads result workbooks from the canonical flat `signal_extraction`
layout:

```text
Data-BIDS/derivatives/signal_extraction/<SIGNAL_EXTRACTION_TAG>/sub-<SUBJECT_ID>/ses-<SESSION_ID>/Results/*_results.xlsx
```

By default, `SIGNAL_EXTRACTION_TAG` is built as
`<DWI_LEVEL>--<ROI_VARIANT>`:

```text
brains:   DWI_LEVEL=den_gr-topup, ROI_VARIANT=plain  -> den_gr-topup--plain
phantoms: DWI_LEVEL=raw,          ROI_VARIANT=manual -> raw--manual
RESULTS_ROOT=Data-BIDS/derivatives/signal_extraction/$SIGNAL_EXTRACTION_TAG
```

Examples:

```bash
DWI_LEVEL=den_gr-topup ROI_VARIANT=plain \
  bash signal_analysis/run_dataset.sh brain ogse ingest

DWI_LEVEL=den_gr-topup ROI_VARIANT=sket1 \
  bash signal_analysis/run_dataset.sh brain ogse ingest

DWI_LEVEL=den_gr-topup ROI_VARIANT=sket2 \
  bash signal_analysis/run_dataset.sh brain ogse ingest
```

To ingest one explicit `Results` folder:

```bash
bash signal_analysis/run_dataset.sh brain ogse \
  --results-root Data-BIDS/derivatives/signal_extraction/den_gr-topup--sket1/sub-MBBL-3/ses-T0/Results \
  ingest
```

To ingest old nested results, set `SIGNAL_EXTRACTION_LAYOUT=nested`, which uses
`Data-BIDS/derivatives/signal_extraction/$DWI_LEVEL/$ROI_VARIANT`. The older
`Results/<subject>` layout is not part of the current pipeline contract. Each
`Results` folder must contain `*_results.xlsx` directly.

### Sequence Parameters

The ingest step matches result rows against sequence-parameter Excel workbooks:

| Dataset | Default workbook |
|---|---|
| `brain` | `Data-BIDS/sequence_parameters_brains.xlsx` |
| `phantom` | `Data-BIDS/sequence_parameters_phantoms.xlsx` |

Override with:

```bash
PARAMS_XLSX=/path/to/sequence_parameters.xlsx \
  bash signal_analysis/run_dataset.sh phantom ogse ingest
```

The ingest step uses the matched `sequence_parameters` row to attach acquisition
metadata and, when needed, to resolve the result-table layout. Layout resolution
is tried in this order:

1. Legacy filename tokens such as `10bval_6dir`.
2. Optional `sequence_parameters` columns: `ndirs` and `nbvals`.
3. The canonical Balseiro signal-extraction table shape: 2 b0 rows followed by
   60 data rows, interpreted as 6 directions by 10 b-values.

For new result filenames such as
`sub-BRAIN-1_ses-T0_acq-hz000d40b2000s5_results.xlsx`, adding `ndirs` and
`nbvals` to the matching sequence-parameter row is the explicit, preferred
metadata path. Existing files without those columns are still ingested through
the table-shape fallback.

### Manifests

Some downstream steps need explicit manifest CSV files. The runner selects the
manifest directory from:

```text
signal_analysis/manifests/<brains|phantoms>_<ogse|nogse>/
```

Current tracked manifests:

```text
manifests/brains_ogse/grad_correction.csv
manifests/phantoms_ogse/grad_correction.csv
manifests/phantoms_ogse/grad_correction-PHANTOM3.csv
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
MANIFEST_DIR=signal_analysis/manifests/phantoms_ogse \
  bash signal_analysis/run_dataset.sh phantom ogse grad_correction
```

## Outputs

Default output roots:

| Dataset | Sequence | Default `ANALYSIS_ROOT` |
|---|---|---|
| `brain` | `ogse` | `analysis/<SIGNAL_EXTRACTION_TAG>--brains/ogse_experiments` |
| `brain` | `nogse` | `analysis/<SIGNAL_EXTRACTION_TAG>--brains/nogse_experiments` |
| `phantom` | `ogse` | `analysis/<SIGNAL_EXTRACTION_TAG>--phantoms/ogse_experiments` |
| `phantom` | `nogse` | `analysis/<SIGNAL_EXTRACTION_TAG>--phantoms/nogse_experiments` |

Typical output layout:

```text
analysis/den_gr-topup--plain--brains/ogse_experiments/
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
│   ├── grad_correction/
│   └── signal_fit_monoexp_value_norm_vs_bvalue_g/
└── alpha_macro/
```

Important outputs:

| Output | Created by | Meaning |
|---|---|---|
| `master.long.parquet` | `ingest` | Canonical long-format master table |
| `master.last_points.long.parquet` | `filter_master_points` | Filtered master table |
| rotated signal rows | `rotate` | Tensor-rotated signal directions appended to the master table |
| `fits/signal_fit_*/*` | `fit_signal` | Monoexponential signal fit parameters, points used, and plots |
| `fits/grad_correction/*` | `grad_correction` | Correction tables, plots, and updated master factors |
| `plots-master/signal/*` | `plot_signal` | Signal-vs-gradient figures |
| `plots-master/monoexp_D_vs_time/*` | `plot_monoexp_d` | Monoexponential D-vs-time figures |
| `alpha_macro/*` | `alpha` | Macroscopic alpha tables and plots |
| `master.xlsx` | `export_master_xlsx` | Inspection-only Excel export |

## Runner Interface

All steps are run through:

```bash
bash signal_analysis/run_dataset.sh <type_subj> <type_seq> <step...>
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
bash signal_analysis/run_dataset.sh --type-subj brain --type-seq ogse ingest
```

Get help:

```bash
bash signal_analysis/run_dataset.sh --help
bash signal_analysis/run_dataset.sh brain ogse ingest --help
bash signal_analysis/run_dataset.sh brain ogse rotate --help
```

Step-specific settings are environment variables placed before the command:

```bash
VAR=value OTHER_VAR=value \
  bash signal_analysis/run_dataset.sh brain ogse <step>
```

Extra Python flags go through the corresponding `*_EXTRA_ARGS` variable:

```bash
GRAD_CORR_EXTRA_ARGS="--avg-N --no-fill-missing" \
  bash signal_analysis/run_dataset.sh brain ogse grad_correction
```

## Available Steps

| Step name | Purpose |
|---|---|
| `ingest` | Read `Results/*_results.xlsx` into `master.long.parquet` |
| `filter_master_points` | Create a filtered master table using last-points rules |
| `rotate` | Rotate signal tensor directions |
| `plot_signal` | Plot signal curves from the master table |
| `fit_signal` | Fit monoexponential signal curves, with `auto_fit_points` by default |
| `fit_signal_gradcorr` | Same as `fit_signal`, but applies embedded gradient-correction factors |
| `grad_correction` | Build and embed gradient-correction factors |
| `alpha` | Build alpha_macro summaries and D-vs-Delta plots |
| `plot_monoexp_d` | Plot monoexponential D vs `td_ms` or `Delta_app_ms` |
| `export_master_xlsx` | Export the selected master parquet to Excel |

The runner also exposes additional contrast and global-fit steps from the shared
code base. Use `--help` to see the full current list.

## How Variables Behave Across Steps

When multiple steps are passed in one command, for example:

```bash
bash signal_analysis/run_dataset.sh brain ogse ingest rotate grad_correction
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
nohup bash signal_analysis/run_dataset.sh brain ogse ingest rotate grad_correction \
  > logs/brain_core.log 2>&1 &
```

Brain OGSE, skeleton ROI:

```bash
DWI_LEVEL=den_gr-topup ROI_VARIANT=sket2 \
nohup bash signal_analysis/run_dataset.sh brain ogse ingest rotate grad_correction \
  > logs/brain_core_sket2.log 2>&1 &
```

Phantom OGSE:

```bash
PARAMS_XLSX=Data-BIDS/sequence_parameters_phantoms.xlsx \
nohup bash signal_analysis/run_dataset.sh phantom ogse ingest rotate grad_correction \
  > logs/phantom_core.log 2>&1 &
```

## Brain OGSE Commands

Unfiltered run:

```bash
mkdir -p logs
nohup bash signal_analysis/run_dataset.sh brain ogse ingest > logs/01_ingest.log 2>&1 &
ROTATE_EXTRA_ARGS="--solver solve" nohup bash signal_analysis/run_dataset.sh brain ogse rotate > logs/02_rotate.log 2>&1 &
PLOT_SIGNAL_YCOL=value nohup bash signal_analysis/run_dataset.sh brain ogse plot_signal > logs/03_plot_signal.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh brain ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
ALPHA_EXTRA_ARGS="--bvalmax 10 --roi-bvalmax Syringe=7 --roi-bvalmax Right-Lateral-Ventricle=5 --roi-bvalmax Left-Lateral-Ventricle=5" DPROJ_DIRS="long tra" nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/06_alpha.log 2>&1 &
```

All-in-one core chain:

```bash
nohup bash signal_analysis/run_dataset.sh brain ogse ingest rotate grad_correction > logs/brain_core.log 2>&1 &
```

`plot_signal`, `alpha`, and `plot_monoexp_d` are usually run separately because
each one often needs its own filters and extra arguments.

## Phantom OGSE Commands

Unfiltered run:

```bash
PARAMS_XLSX=Data-BIDS/sequence_parameters_phantoms.xlsx nohup bash signal_analysis/run_dataset.sh phantom ogse ingest > logs/01_ingest.log 2>&1 &
ROTATE_EXTRA_ARGS="--solver solve" nohup bash signal_analysis/run_dataset.sh phantom ogse rotate > logs/02_rotate.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh phantom ogse plot_signal > logs/03_plot_signal.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh phantom ogse grad_correction > logs/05_grad_correction.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh phantom ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh phantom ogse alpha > logs/06_alpha.log 2>&1 &
```

All-in-one core chain:

```bash
PARAMS_XLSX=Data-BIDS/sequence_parameters_phantoms.xlsx \
  nohup bash signal_analysis/run_dataset.sh phantom ogse ingest rotate grad_correction > logs/phantom_core.log 2>&1 &
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
nohup bash signal_analysis/run_dataset.sh brain ogse ingest > logs/01_ingest.log 2>&1 &
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" nohup bash signal_analysis/run_dataset.sh brain ogse filter_master_points > logs/00_filter.log 2>&1 &
MASTER_PARQUET=analysis/den_gr-topup--plain--brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh brain ogse rotate > logs/02_rotate.log 2>&1 &
MASTER_PARQUET=analysis/den_gr-topup--plain--brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh brain ogse plot_signal > logs/03_plot_signal.log 2>&1 &
MASTER_PARQUET=analysis/den_gr-topup--plain--brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh brain ogse grad_correction > logs/05_grad_correction.log 2>&1 &
MASTER_PARQUET=analysis/den_gr-topup--plain--brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh brain ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
MASTER_PARQUET=analysis/den_gr-topup--plain--brains/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/06_alpha.log 2>&1 &
```

Brain OGSE all-in-one filtered core chain:

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/run_dataset.sh brain ogse ingest filter_master_points rotate grad_correction > logs/brain_filtered.log 2>&1 &
```

Phantom OGSE:

```bash
nohup bash signal_analysis/run_dataset.sh phantom ogse ingest > logs/01_ingest.log 2>&1 &
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" nohup bash signal_analysis/run_dataset.sh phantom ogse filter_master_points > logs/00_filter.log 2>&1 &
MASTER_PARQUET=analysis/raw--manual--phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh phantom ogse rotate > logs/02_rotate.log 2>&1 &
MASTER_PARQUET=analysis/raw--manual--phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh phantom ogse plot_signal > logs/03_plot_signal.log 2>&1 &
MASTER_PARQUET=analysis/raw--manual--phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh phantom ogse grad_correction > logs/05_grad_correction.log 2>&1 &
MASTER_PARQUET=analysis/raw--manual--phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh phantom ogse plot_monoexp_d > logs/06b_plot_monoexp_d.log 2>&1 &
MASTER_PARQUET=analysis/raw--manual--phantoms/ogse_experiments/master.last_points.long.parquet nohup bash signal_analysis/run_dataset.sh phantom ogse alpha > logs/06_alpha.log 2>&1 &
```

Phantom OGSE all-in-one filtered core chain:

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/run_dataset.sh phantom ogse ingest filter_master_points rotate grad_correction > logs/phantom_filtered.log 2>&1 &
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
nohup bash signal_analysis/run_dataset.sh brain nogse ingest > logs/brain_nogse_01_ingest.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh brain nogse rotate > logs/brain_nogse_02_rotate.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh brain nogse plot_signal > logs/brain_nogse_03_plot_signal.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh brain nogse alpha > logs/brain_nogse_06_alpha.log 2>&1 &
```

Brain NOGSE all-in-one core chain:

```bash
nohup bash signal_analysis/run_dataset.sh brain nogse ingest rotate > logs/brain_nogse_core.log 2>&1 &
```

Brain NOGSE filtered core chain:

```bash
MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/run_dataset.sh brain nogse ingest filter_master_points rotate > logs/brain_nogse_filtered.log 2>&1 &
```

Phantom NOGSE:

```bash
NOGSE_ONEG=1 nohup bash signal_analysis/run_dataset.sh phantom nogse ingest > logs/phantom_nogse_01_ingest.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh phantom nogse rotate > logs/phantom_nogse_02_rotate.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh phantom nogse plot_signal > logs/phantom_nogse_03_plot_signal.log 2>&1 &
nohup bash signal_analysis/run_dataset.sh phantom nogse alpha > logs/phantom_nogse_06_alpha.log 2>&1 &
```

Use `NOGSE_ONEG=1` for direct-g NOGSE phantom results where each
`*_results.xlsx` file contributes one gradient point to a shared curve. Use
`NOGSE_ONEG=0` or omit it only when each results file already contains a full
multi-point curve. This flag is only needed during `ingest`; later steps read
the `one_g_per_sequence` marker from `master.long.parquet`.

Phantom NOGSE all-in-one core chain:

```bash
NOGSE_ONEG=1 nohup bash signal_analysis/run_dataset.sh phantom nogse ingest rotate > logs/phantom_nogse_core.log 2>&1 &
```

Phantom NOGSE filtered core chain:

```bash
NOGSE_ONEG=1 MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
  nohup bash signal_analysis/run_dataset.sh phantom nogse ingest filter_master_points rotate > logs/phantom_nogse_filtered.log 2>&1 &
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
  nohup bash signal_analysis/run_dataset.sh brain ogse plot_signal > logs/plot_signal.log 2>&1 &

PLOT_SUBJ=20220622_BRAIN PLOT_DIRECTION=long PLOT_SIGNAL_XCOL=g_thorsten \
  nohup bash signal_analysis/run_dataset.sh brain ogse plot_signal > logs/plot_signal.log 2>&1 &
```

### `fit_signal`

Fits signal curves from `master.long.parquet`. At the moment this step supports
only the monoexponential model:

```text
S(b) = M0 * exp(-b * D0)
```

The default input row kind is `signal_rotated`, so the usual order is:

```bash
bash signal_analysis/run_dataset.sh brain ogse ingest rotate fit_signal
```

By default, `fit_signal` uses `--auto_fit_points`. The automatic selection is
not an arbitrary subset search. For each curve it sorts points by `b_step`, fits
leading prefixes, and keeps the largest prefix whose `rmse_log` remains within
tolerance when the next point is added. The defaults are:

| Option | Default | Meaning |
|---|---:|---|
| `--auto_fit_min_points` | `3` | First prefix size tested |
| `--auto_fit_max_points` | `9` | Last prefix size tested |
| `--auto_fit_tol` | `0.05` | Relative tolerance when comparing consecutive prefixes |
| `--auto_fit_err_floor` | `0.005` | Minimum previous `rmse_log` used in the comparison |

To use a fixed number of leading points instead, clear the default auto-fit
args and pass `--fit_points`:

```bash
SIGNAL_FIT_AUTO_ARGS="" \
SIGNAL_FIT_EXTRA_ARGS="--fit_points 6" \
  bash signal_analysis/run_dataset.sh brain ogse fit_signal
```

Input selection can come from `signal_analysis/manifests/<dataset>/signal_fits.csv`.
If that manifest is missing, the step warns and fits the selected master rows
directly. A `signal_fits.csv` manifest may contain:

| Column | Meaning |
|---|---|
| `subj` | Subject/group selector, or `ALL` |
| `sheet` | Sequence-parameter sheet selector, or `ALL` |
| `roi` | ROI selector, or `ALL` |
| `direction` | Direction selector, or `ALL` |
| `td_ms` | Diffusion time selector |
| `N` | Oscillation count selector |
| `Hz` | Frequency selector |
| `model` | Must be `monoexp` for now |

Main variables:

| Variable | Default | Description |
|---|---|---|
| `FIT_SIGNAL_SCRIPT` | `$REPO_ROOT/scripts/fitting/fit_signal.py` | Python script |
| `SIGNAL_FIT_MANIFEST` | `$MANIFEST_DIR/signal_fits.csv` | Optional curve-selection manifest |
| `SIGNAL_FIT_MODEL` | `monoexp` | Only `monoexp` is supported for now |
| `SIGNAL_FIT_YCOL` | `value_norm` | Signal column to fit |
| `SIGNAL_FIT_B_AXIS` | OGSE: `bvalue_g`; NOGSE: `g` | Fit x-axis. `SIGNAL_FIT_G_TYPE` is accepted as an alias |
| `SIGNAL_FIT_AUTO_ARGS` | `--auto_fit_points` | Point-selection mode |
| `SIGNAL_FIT_EXTRA_ARGS` | none | Extra arguments for `fit_signal.py` |
| `SIGNAL_FIT_OUT_ROOT` | `$ANALYSIS_ROOT/fits/signal_fit_<model>_<ycol>_vs_<axis>` | Output directory |
| `SIGNAL_FIT_PLOT_DIR` | `$SIGNAL_FIT_OUT_ROOT/plots` | Fit-plot directory |

Outputs:

```text
$SIGNAL_FIT_OUT_ROOT/
├── signal_fit_params.monoexp.<axis>.<ycol>.parquet
├── signal_fit_params.monoexp.<axis>.<ycol>.xlsx
├── signal_fit_points.monoexp.<axis>.<ycol>.parquet
├── signal_fit_points.monoexp.<axis>.<ycol>.xlsx
└── plots/
```

`signal_fit_params` stores one row per fitted curve, including `D0_mm2_s`,
`M0`, `fit_points`, `n_fit`, `fit_strategy`, `auto_fit_score`, `rmse_log`, and
the fit status message. `signal_fit_points` stores one row per signal point and
marks `used_for_fit=True` for the points included in the selected prefix.

Examples:

```bash
nohup bash signal_analysis/run_dataset.sh brain ogse fit_signal \
  > logs/04_fit_signal.log 2>&1 &

SIGNAL_FIT_EXTRA_ARGS="--auto_fit_tol 0.10 --auto_fit_max_points 10" \
  nohup bash signal_analysis/run_dataset.sh brain ogse fit_signal \
  > logs/04_fit_signal.log 2>&1 &

SIGNAL_FIT_EXTRA_ARGS="--free_M0" \
  nohup bash signal_analysis/run_dataset.sh brain ogse fit_signal \
  > logs/04_fit_signal.log 2>&1 &
```

After `grad_correction` has embedded `grad_correction_factor` in the master
table, run the corrected preset:

```bash
nohup bash signal_analysis/run_dataset.sh brain ogse fit_signal_gradcorr \
  > logs/04_fit_signal_gradcorr.log 2>&1 &
```

That preset is equivalent to adding `--apply_grad_corr`; it scales the b-axis
by `correction_factor^2` before fitting. For OGSE it also defaults to
`--directions long tra` unless a direction filter is already present, because
the gradient-correction table is normally defined only for those rotated
directions.

### `plot_monoexp_d`

Plots monoexponential diffusion estimates from signal-fit parquets. The usual
input is the output of `fit_signal`, but externally generated monoexp fit
tables can also be used if they follow the expected schema and live under
`SIGNAL_FITS_ROOT`.

| Variable | Default | Description |
|---|---|---|
| `SIGNAL_FITS_ROOT` | `$ANALYSIS_ROOT/fits` | Root scanned for signal fits. Set this to a specific `fit_signal` output directory to narrow the scan |
| `MONOEXP_D_OUT_DIR` | `$ANALYSIS_ROOT/plots-master/monoexp_D_vs_time` | Output directory |
| `PLOT_MONOEXP_D_EXTRA_ARGS` | none | Extra arguments for `plot_monoexp_D_vs_time.py` |

Examples:

```bash
nohup bash signal_analysis/run_dataset.sh brain ogse plot_monoexp_d > logs/plot_monoexp_d.log 2>&1 &

SIGNAL_FITS_ROOT=analysis/den_gr-topup--plain--brains/ogse_experiments/fits/signal_fit_monoexp_value_norm_vs_bvalue_g \
  nohup bash signal_analysis/run_dataset.sh brain ogse plot_monoexp_d > logs/plot_monoexp_d.log 2>&1 &
```

### `grad_correction`

Fits `D0_nogse` and `D0_monoexp` on signal curves from the syringe or another
reference phantom ROI listed in the manifest, then computes:

```text
correction_factor = sqrt(D0_nogse / D0_monoexp_avg)
```

By default, the NOGSE fit uses `g` as the gradient axis and the monoexp fit
uses `bvalue_g` as the b-value axis.

How `D0_monoexp` is averaged:

- The monoexp fit uses `auto_fit_points` by default in this step. It follows the
  same leading-prefix rule described in `fit_signal`: sort by `b_step`, test
  increasing prefix sizes, and keep the largest prefix that remains within the
  `rmse_log` tolerance.
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
| `GRAD_CORR_AUTO_FIT_ARGS` | `--auto_fit_points` | Auto-fit flags for the monoexp `D0` used in the correction |
| `GRAD_CORR_EXTRA_ARGS` | none | Extra arguments for `make_grad_correction_table.py` |

Useful `GRAD_CORR_EXTRA_ARGS`: `--avg-N`, `--avg-N 4 8`,
`--no-fill-missing`, `--row-kind signal_rotated`, `--stat avg`, `--free-M0`,
`--gbase g`, `--bbase bvalue_g`, `--auto_fit_tol 0.05`,
`--auto_fit_min_points 3`, `--auto_fit_max_points 9`.

Examples:

```bash
nohup bash signal_analysis/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_EXTRA_ARGS="--avg-N" \
  nohup bash signal_analysis/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_EXTRA_ARGS="--avg-N 1 4 8" \
  nohup bash signal_analysis/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_EXTRA_ARGS="--no-fill-missing" \
  nohup bash signal_analysis/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_AUTO_FIT_ARGS="--no_auto_fit_points" \
  nohup bash signal_analysis/run_dataset.sh brain ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &

GRAD_CORR_ROI=Water2 \
  nohup bash signal_analysis/run_dataset.sh phantom ogse grad_correction \
  > logs/05_grad_correction.log 2>&1 &
```

### `alpha`

Computes macroscopic alpha from `D_proj` values in `signal_rotated` master rows.
Also generates D-vs-Delta_app plots.

| Variable | Default | Description |
|---|---|---|
| `ALPHA_N` | `1` | N selector for `make_alpha_macro_summary.py` |
| `ALPHA_OUT_DIR` | `$ANALYSIS_ROOT/alpha_macro/master` | Output directory |
| `ALPHA_PLOT_BSTEPS` | none | Candidate b-value positions used for alpha selection and D-vs-Delta plots |
| `ALPHA_PLOT_BVALUES` | none | Candidate rounded b-values used for alpha selection and D-vs-Delta plots |
| `ALPHA_EXTRA_ARGS` | none | Extra arguments for `make_alpha_macro_summary.py` |
| `DPROJ_N` | same as `ALPHA_N` | N selector for D-vs-Delta plots when it should differ from `ALPHA_N` |
| `DPROJ_DIRS` | none | Direction filter for D-vs-Delta plots |
| `DPROJ_ROIS` | none | ROI filter for D-vs-Delta plots |
| `PLOT_D0_EXTRA_ARGS` | none | Extra arguments only for `plot_D0_vs_Delta.py` |

Useful `ALPHA_EXTRA_ARGS`: `--bvalmax 5`, `--roi-bvalmax AntCC=7`,
`--dirs long tra`, `--subjs MBBL LUDG`, `--sheets 20230630_MBBL-3`.
Add `--no-annotate-master-alpha` when testing a subset and you do not want the
step to write `alpha_macro` back into `master.long.parquet`.

Subject filters:

- `--subjs` filters the master-table `subj` column. In the brain table these
  are family labels such as `ADBN`, `ARVE`, `BRAIN`, `LUDG`, `MBBL`, `SNVN`.
- `--sheets` filters the acquisition/session label, for example
  `20230630_MBBL-3` or `20230710_LUDG-3`.

Use `ALPHA_PLOT_BSTEPS` or `ALPHA_PLOT_BVALUES` when the b-values used for
alpha should be chosen from the same subset that appears in the D-vs-Delta
plots:

- `ALPHA_PLOT_BSTEPS="1 3 5"` keeps only those 1-based b-value positions after
  sorting b-values ascending within each group.
- `ALPHA_PLOT_BVALUES="500 980 1280 2000"` keeps only those rounded b-values.
  Values are compared after the `--bvalue-decimals` rounding used by the script.

After the candidate b-values are defined, `--bvalmax` and `--roi-bvalmax` choose
the b-value used to compute `alpha_macro` and the horizontal annotation. The
selector can be either a candidate bstep or a b-value:

- `--bvalmax 4` means “use the 4th candidate b-value”.
- `--bvalmax 2000` means “use b=2000”.
- `--roi-bvalmax Syringe=3` means “for Syringe, use the 3rd candidate b-value”.
- `--roi-bvalmax Syringe=1280` means “for Syringe, use b=1280”.

`PLOT_D0_EXTRA_ARGS="--plot-bsteps ..."` and
`PLOT_D0_EXTRA_ARGS="--plot-bvalues ..."` are still available for plotting-only
filters, but they do not affect the alpha calculation. Prefer
`ALPHA_PLOT_BSTEPS` / `ALPHA_PLOT_BVALUES` when alpha and plots should use the
same b-value subset.

Outputs under `$ALPHA_OUT_DIR/`:

- `summary_alpha_values.xlsx` - alpha_macro by subject, ROI, and direction
- `D_vs_delta_app.combined.xlsx` - aggregated D vs Delta_app
- `alpha_macro_vs_roi.png` - summary plot by ROI
- `<subj>/<roi>/<dir>_*.png` - individual curves

Examples:

```bash
ALPHA_N=1 ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Left-Lateral-Ventricle=5 --roi-bvalmax Right-Lateral-Ventricle=5 --roi-bvalmax Syringe=7 --dirs long tra" nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

ALPHA_N=1 ALPHA_EXTRA_ARGS="--subjs MBBL LUDG --bvalmax 4 --roi-bvalmax Left-Lateral-Ventricle=2 --roi-bvalmax Right-Lateral-Ventricle=2 --roi-bvalmax Syringe=3 --dirs long tra" nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

ALPHA_N=1 ALPHA_EXTRA_ARGS="--sheets 20230630_MBBL-3 20230710_LUDG-3 --bvalmax 4 --roi-bvalmax Left-Lateral-Ventricle=2 --roi-bvalmax Right-Lateral-Ventricle=2 --roi-bvalmax Syringe=3 --dirs long tra" nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

ALPHA_N=1 ALPHA_EXTRA_ARGS="--subjs MBBL --bvalmax 4 --dirs long tra" PLOT_D0_EXTRA_ARGS="--plot-bsteps 1 2 4" nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

ALPHA_N=1 ALPHA_EXTRA_ARGS="--subjs MBBL --bvalmax 4 --dirs long tra" PLOT_D0_EXTRA_ARGS="--plot-bvalues 250 500 750" nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

ALPHA_N=1 ALPHA_PLOT_BVALUES="500 980 1280 2000" ALPHA_EXTRA_ARGS="--subjs BRAIN MBBL LUDG --bvalmax 4 --roi-bvalmax Left-Lateral-Ventricle=1 --roi-bvalmax Right-Lateral-Ventricle=1 --roi-bvalmax Syringe=3 --dirs long tra --rois Left-Lateral-Ventricle Right-Lateral-Ventricle Syringe" nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

ALPHA_N=1 ALPHA_PLOT_BVALUES="500 980 1280 2000" ALPHA_EXTRA_ARGS="--subjs BRAIN MBBL LUDG --bvalmax 2000 --roi-bvalmax Left-Lateral-Ventricle=500 --roi-bvalmax Right-Lateral-Ventricle=500 --roi-bvalmax Syringe=1280 --dirs long tra --rois Left-Lateral-Ventricle Right-Lateral-Ventricle Syringe" nohup bash signal_analysis/run_dataset.sh brain ogse alpha > logs/alpha.log 2>&1 &

ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Syringe=7" nohup bash signal_analysis/run_dataset.sh phantom ogse alpha > logs/alpha.log 2>&1 &
```

`summary_alpha_values.xlsx` is a typical input for downstream tc fitting steps
outside this repository, for example `TC_METHOD=pseudohuber_fixed_macro`.

## Export Master Table to Excel

```bash
MASTER_PARQUET=analysis/den_gr-topup--plain--brains/ogse_experiments/master.long.parquet MASTER_XLSX=analysis/den_gr-topup--plain--brains/ogse_experiments/master.xlsx \
  nohup bash signal_analysis/run_dataset.sh brain ogse export_master_xlsx > logs/export_master_xlsx.log 2>&1 &

MASTER_PARQUET=analysis/raw--manual--phantoms/ogse_experiments/master.long.parquet MASTER_XLSX=analysis/raw--manual--phantoms/ogse_experiments/master.xlsx \
  nohup bash signal_analysis/run_dataset.sh phantom ogse export_master_xlsx > logs/export_master_xlsx.log 2>&1 &
```

## Global Dataset Variables

All variables are inherited from `master_table_common.sh` and can be set before
any command:

| Variable | Description |
|---|---|
| `PY` | Python interpreter |
| `SIGNALS_ROOT` | Root containing sequence-parameter Excel files; default `$PROJECT_ROOT/Data-BIDS` |
| `DWI_LEVEL` | `signal_extraction` DWI level used by `ingest`; defaults: brains `den_gr-topup`, phantoms `raw` |
| `ROI_VARIANT` | segmentation/ROI namespace for brains, manual-mask namespace for phantoms; defaults: brains `plain`, phantoms `manual` |
| `SIGNAL_EXTRACTION_TAG` | flat folder under `Data-BIDS/derivatives/signal_extraction`; default `$DWI_LEVEL--$ROI_VARIANT` |
| `SIGNAL_EXTRACTION_LAYOUT` | `flat` by default; set `nested` to read legacy `$DWI_LEVEL/$ROI_VARIANT` results |
| `RESULTS_ROOT` / `--results-root` repeatable | Results folder or folders to ingest |
| `PARAMS_XLSX` | Sequence-parameter Excel file; required for phantom `ingest` |
| `NOGSE_ONEG` | Set to `1` during NOGSE direct-g `ingest` when each results file is one gradient point |
| `ANALYSIS_TAG` | analysis namespace; default `$SIGNAL_EXTRACTION_TAG--<brains\|phantoms>` |
| `ANALYSIS_ROOT` | Analysis output root; default `$PROJECT_ROOT/analysis/$ANALYSIS_TAG/<experiment>` |
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
bash -n signal_analysis/run_dataset.sh
bash -n signal_analysis/steps/01_ingest_results.sh
bash -n signal_analysis/steps/02_rotate_signals.sh
bash -n signal_analysis/steps/05_make_grad_correction.sh
```
