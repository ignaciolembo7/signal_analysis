#!/usr/bin/env bash

pipeline_setup_common() {
    MASTER_HELPER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_ROOT="$(cd "$MASTER_HELPER_DIR/.." && pwd)"
    TEMPLATE_ROOT="$REPO_ROOT"
    PROJECT_ROOT="$(cd "$REPO_ROOT/.." && pwd)"

    export PYTHONPATH="$REPO_ROOT/src:${PYTHONPATH:-}"
    export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"

    local default_py="python"
    if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
        default_py="${CONDA_PREFIX}/bin/python"
    elif [[ -x "$HOME/.conda/envs/nogse_pipe_env/bin/python" ]]; then
        default_py="$HOME/.conda/envs/nogse_pipe_env/bin/python"
    elif command -v python3 >/dev/null 2>&1; then
        default_py="$(command -v python3)"
    fi
    PY="${PY:-$default_py}"
}

pipeline_normalize_type_subj() {
    case "$1" in
        brain|brains) echo "brains" ;;
        phantom|phantoms) echo "phantoms" ;;
        *)
            echo "ERROR: unknown type_subj '$1' (use brain or phantom)" >&2
            return 2
            ;;
    esac
}

pipeline_normalize_type_seq() {
    case "$1" in
        ogse|nogse) echo "$1" ;;
        *)
            echo "ERROR: unknown type_seq '$1' (use ogse or nogse)" >&2
            return 2
            ;;
    esac
}

pipeline_set_dataset_defaults() {
    local dataset
    local type_seq
    dataset="$(pipeline_normalize_type_subj "$1")" || exit 2
    type_seq="$(pipeline_normalize_type_seq "${2:-${TYPE_SEQ:-ogse}}")" || exit 2

    DATASET="$dataset"
    TYPE_SUBJ="${dataset%s}"
    TYPE_SEQ="$type_seq"
    EXPERIMENT_ROOT_NAME="${type_seq}_experiments"
    SIGNALS_ROOT="${SIGNALS_ROOT:-$PROJECT_ROOT/Data-BIDS}"
    DWI_LEVEL="${DWI_LEVEL:-den_gr-topup}"
    ROI_VARIANT="${ROI_VARIANT:-plain}"
    RESULTS_ROOT="${RESULTS_ROOT:-$PROJECT_ROOT/Data-BIDS/derivatives/signal_extraction/$DWI_LEVEL/$ROI_VARIANT}"

    case "$dataset" in
        brains)
            PARAMS_XLSX="${PARAMS_XLSX:-$SIGNALS_ROOT/sequence_parameters_brains.xlsx}"
            ANALYSIS_ROOT="${ANALYSIS_ROOT:-$PROJECT_ROOT/analysis/brains/$EXPERIMENT_ROOT_NAME}"
            ;;
        phantoms)
            PARAMS_XLSX="${PARAMS_XLSX:-$SIGNALS_ROOT/sequence_parameters_phantoms.xlsx}"
            ANALYSIS_ROOT="${ANALYSIS_ROOT:-$PROJECT_ROOT/analysis/phantoms/$EXPERIMENT_ROOT_NAME}"
            ;;
    esac

    MASTER_PARQUET="${MASTER_PARQUET:-$ANALYSIS_ROOT/master.long.parquet}"
    MANIFEST_DIR="${MANIFEST_DIR:-$TEMPLATE_ROOT/manifests/${dataset}_${type_seq}}"

    export DATASET TYPE_SUBJ TYPE_SEQ EXPERIMENT_ROOT_NAME DWI_LEVEL ROI_VARIANT
}

pipeline_require_file() {
    local path="$1"
    local label="$2"
    if [[ ! -f "$path" ]]; then
        echo "ERROR: $label not found: $path" >&2
        exit 1
    fi
}

pipeline_apply_master_last_points() {
    local spec="${MASTER_LAST_POINTS_BY_TD:-}"
    if [[ -z "${spec// }" || "${spec^^}" == "ALL" ]]; then
        return 0
    fi
    if [[ -n "${MASTER_LAST_POINTS_APPLIED:-}" ]]; then
        return 0
    fi

    local filter_script="$REPO_ROOT/scripts/data/filter_master_table.py"
    pipeline_require_file "$MASTER_PARQUET" "master table"
    pipeline_require_file "$filter_script" "master filter script"

    local out_dir
    out_dir="${MASTER_LAST_POINTS_DIR:-$ANALYSIS_ROOT}"
    mkdir -p "$out_dir"

    MASTER_PARQUET_ORIGINAL="${MASTER_PARQUET_ORIGINAL:-$MASTER_PARQUET}"
    MASTER_PARQUET="${MASTER_LAST_POINTS_PARQUET:-$out_dir/master.last_points.long.parquet}"

    "$PY" "$filter_script" "$MASTER_PARQUET_ORIGINAL" \
        --out-parquet "$MASTER_PARQUET" \
        --last-points-by-td-n "$spec"

    MASTER_LAST_POINTS_APPLIED=1
    export MASTER_PARQUET MASTER_PARQUET_ORIGINAL MASTER_LAST_POINTS_APPLIED
    echo "Using filtered master table: $MASTER_PARQUET"
    echo "Last points by td_ms/N: $spec"
}

pipeline_usage() {
    cat <<'EOF'
Usage:
  bash signal_analysis/run_dataset.sh <type_subj> <type_seq> <step...>

Run from:
  /mnt/storage/tier2/MUMI-EXT-001/mumi-data/Project-Balseiro-Microstructure

Subject and sequence:
  <type_subj>  brain | phantom
               Aliases accepted: brains, phantoms
  <type_seq>   ogse | nogse

Equivalent option form:
  bash signal_analysis/run_dataset.sh --type-subj brain --type-seq ogse <step...>

Runner arguments:
  --type-subj brain|phantom   Select subject type.
  --type_subj brain|phantom   Select subject type, underscore spelling.
  --dataset brain|phantom     Alias for --type-subj.
  --type-seq ogse|nogse       Select sequence type.
  --type_seq ogse|nogse       Select sequence type, underscore spelling.
  --results-root DIR          Add one Results input folder for ingest. Repeatable.
  -h, --help                  Show this help, or step help after a step name.

How step arguments work:
  The runner only parses <type_subj>, <type_seq>, step names, and --results-root.
  Step-specific settings are environment variables placed before the command:

    VAR=value OTHER_VAR=value \
      bash signal_analysis/run_dataset.sh brain ogse <step>

  Extra Python flags go through each step's *_EXTRA_ARGS variable:

    SIGNAL_FIT_EXTRA_ARGS="--fix_M0 1.0 --auto_fit_tol 0.05" \
      bash signal_analysis/run_dataset.sh brain ogse fit_signal

Available steps:
  Data import:
    ingest           Read Results/*_results.xlsx into master.long.parquet.

  Master-table construction:
    filter_master_points
                     Write an intermediate master table with first points by td_ms.
    rotate           Rotate signal tensor directions.
    export_master_xlsx
                     Export the selected master parquet to an Excel workbook.

  Plots:
    plot_signal      Plot signal curves from master.long.parquet.
    plot_monoexp_d   Plot monoexponential D vs td_ms/Delta_app_ms.

  Fits:
    fit_signal                 Fit signal curves using signal_fits.csv.
    fit_signal_gradcorr        fit_signal with gradient correction enabled.
    fit_global_signal          Fit mixed/global signal models directly.
    fit_contrast               Fit contrast rows.
    fit_contrast_free          fit_contrast with free-model defaults.
    fit_contrast_mixed_global  fit_contrast with ogse_mixed_global model (OGSE only).

  Contrast construction:
    contrast           Build direct contrast rows (point-by-point subtraction).
    contrast_resampled Build resampled contrasts from signal fit curves (requires fit_signal).

  Summaries and analysis:
    grad_correction  Build and embed gradient-correction factors.
    alpha            Build alpha_macro summaries and D vs Delta plots.
    extract_tc_peak  Consolidate fit_params, extract tc_peak table, and plot panels.
    tc               Fit tc-vs-td summaries.

Help for one step:
  bash signal_analysis/run_dataset.sh brain ogse ingest --help
  bash signal_analysis/run_dataset.sh brain ogse rotate --help
  bash signal_analysis/run_dataset.sh brain ogse fit_signal --help

Common environment variables:
  PY                 Python interpreter.
  SIGNALS_ROOT       Root containing sequence parameter workbooks. Default: Data-BIDS.
  DWI_LEVEL          signal_extraction DWI level. Default: den_gr-topup.
  ROI_VARIANT        signal_extraction ROI variant. Default: plain.
  PARAMS_XLSX        Sequence-parameter workbook.
  ANALYSIS_ROOT      Output analysis root.
  MASTER_PARQUET     Master table path.
  MANIFEST_DIR       Directory with contrasts.csv, signal_fits.csv, and grad_correction.csv.
  MASTER_LAST_POINTS_BY_TD
                    Optional TD=POINTS or TD:N=POINTS rules applied before post-ingest steps.
                    Keeps the last POINTS b_step values per (td_ms[, N]) group.
                    Example: MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8"
                    Unlisted td_ms/N combinations keep all points.
  MASTER_LAST_POINTS_PARQUET
                    Optional filtered master output path. Default:
                    $ANALYSIS_ROOT/master.last_points.long.parquet.
  TC_FIT_PARAMS     Fit-params parquet used by the tc step. Defaults to the
                    tc_peak_table.parquet produced by extract_tc_peak. Override
                    only if pointing to a non-default location.

Master table format:
  master.long.parquet is the canonical master table and the pipeline only
  appends to parquet. To inspect it in Excel, export a copy explicitly:

    python signal_analysis/scripts/data/export_master_table.py \
      analysis/brains/ogse_experiments/master.long.parquet \
      --out-xlsx analysis/brains/ogse_experiments/master.inspect.xlsx

Examples:
  # Ingest one signal_extraction Results folder.
  bash signal_analysis/run_dataset.sh brain ogse \
    --results-root Data-BIDS/derivatives/signal_extraction/den_gr-topup/sket1/sub-MBBL-3/ses-T0/Results \
    ingest

  # Or ingest all Results folders under one DWI/ROI namespace.
  DWI_LEVEL=den_gr-topup ROI_VARIANT=sket1 \
    bash signal_analysis/run_dataset.sh brain ogse ingest

  # Run the next core steps after ingest.
  bash signal_analysis/run_dataset.sh brain ogse rotate contrast

  # Filter a rotate step by subject/sheet.
  MASTER_SUBJ=BRAIN MASTER_SHEET=20220622_BRAIN \
    bash signal_analysis/run_dataset.sh brain ogse rotate

  # Plot one ROI/direction.
  PLOT_ROI=Left-Lateral-Ventricle PLOT_DIRECTION=long \
    bash signal_analysis/run_dataset.sh brain ogse plot_signal

  # Fit signals with extra Python options.
  SIGNAL_FIT_MODEL=monoexp \
  SIGNAL_FIT_EXTRA_ARGS="--fix_M0 1.0" \
    bash signal_analysis/run_dataset.sh brain ogse fit_signal

  # Run several steps from PIPELINE_STEPS instead of positional step names.
  PIPELINE_STEPS="rotate contrast fit_signal fit_contrast extract_tc_peak tc" \
    bash signal_analysis/run_dataset.sh brain ogse
EOF
}

pipeline_step_help() {
    local step="$1"
    case "$step" in
        ingest)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> ingest
  bash run_dataset.sh <type_subj> <type_seq> --results-root RESULTS_ROOT ingest

What it does:
  Reads *_results.xlsx files from BIDS-session Results folders, uses the
  sequence-parameter workbook, writes long signal rows, and appends
  row_kind='signal' rows to master.long.parquet.

Runner arguments:
  --results-root DIR  Add one input folder through run_dataset.sh. Repeat for multiple folders.

Variables for this step:
  RESULTS_ROOT        One Results folder/root if RESULTS_ROOTS is unset.
                      Default: Data-BIDS/derivatives/signal_extraction/$DWI_LEVEL/$ROI_VARIANT
  RESULTS_ROOTS       Space-separated Results roots.
  DWI_LEVEL           signal_extraction DWI level. Default: den_gr-topup.
  ROI_VARIANT         signal_extraction ROI variant. Default: plain.
  PARAMS_XLSX         Sequence-parameter workbook.
  RESULTS_GLOB        Results filename pattern. Default: *_results.xlsx
  MASTER_PARQUET      Master table output.
  PROCESS_SCRIPT      Python script override.
  PROCESS_OUT_ROOT    Per-file table output root. Default: $ANALYSIS_ROOT/data/tables

Examples:
  bash signal_analysis/run_dataset.sh brain ogse \
    --results-root Data-BIDS/derivatives/signal_extraction/den_gr-topup/sket1/sub-MBBL-3/ses-T0/Results \
    ingest

  RESULTS_ROOTS="Data-BIDS/derivatives/signal_extraction/den_gr-topup/sket1/sub-MBBL-3/ses-T0/Results Data-BIDS/derivatives/signal_extraction/den_gr-topup/sket1/sub-LUDG-3/ses-T0/Results" \
    bash signal_analysis/run_dataset.sh brain ogse ingest

  DWI_LEVEL=den_gr-topup ROI_VARIANT=sket1 \
    bash signal_analysis/run_dataset.sh brain ogse ingest

  PARAMS_XLSX=Data-BIDS/sequence_parameters_brains.xlsx \
  RESULTS_GLOB="*_results.xlsx" \
    bash signal_analysis/run_dataset.sh brain ogse ingest
EOF
            ;;
        filter_master_points)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> filter_master_points

What it does:
  Writes an intermediate master table that keeps only the last N b_step values
  per (td_ms, N) group. The original master table is not modified.

  Rules have the form TD=POINTS (all N values at that td_ms) or TD:N=POINTS
  (specific td_ms AND N). More-specific TD:N rules take precedence over TD-only
  rules when both match the same row.

Variables for this step:
  MASTER_PARQUET              Input master table.
  MASTER_LAST_POINTS_BY_TD    TD=POINTS or TD:N=POINTS rules.
                              Example: "120:8=6,120:4=4,210=8,90=ALL".
                              Unlisted td_ms/N combinations keep all points.
  FILTERED_MASTER_PARQUET     Output filtered master table path.
                              Also accepted as MASTER_LAST_POINTS_PARQUET.
                              Default: $ANALYSIS_ROOT/master.last_points.long.parquet
  MASTER_LAST_POINTS_DIR      Output directory override. Default: $ANALYSIS_ROOT
  FILTER_MASTER_SCRIPT        Python script override.

Note:
  After running filter_master_points, pass the output path as MASTER_PARQUET
  to subsequent steps:

Examples:
  # Keep last 6 points at td=120 N=8, last 4 at td=120 N=4, last 8 at td=210 (all N)
  MASTER_LAST_POINTS_BY_TD="120:8=6,120:4=4,210=8" \
    bash signal_analysis/run_dataset.sh brain ogse filter_master_points

  # Keep all points at td=90, last 6 at td=120 (all N)
  MASTER_LAST_POINTS_BY_TD="120=6,90=ALL" \
    bash signal_analysis/run_dataset.sh brain ogse filter_master_points

  MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet \
    bash signal_analysis/run_dataset.sh brain ogse rotate contrast fit_signal
EOF
            ;;
        export_master_xlsx)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> export_master_xlsx

What it does:
  Exports MASTER_PARQUET to an inspection .xlsx workbook. It does not modify the
  parquet master table.

Variables for this step:
  MASTER_PARQUET           Input master table.
  MASTER_XLSX              Output Excel path. Default: MASTER_PARQUET with .xlsx suffix.
  EXPORT_MASTER_ROW_KIND   Optional row_kind filter. Space- or comma-separated.
  EXPORT_MASTER_HEAD       Optional number of first rows to export.
  EXPORT_MASTER_SCRIPT     Python script override.

Examples:
  bash signal_analysis/run_dataset.sh brain ogse export_master_xlsx

  MASTER_PARQUET=analysis/brains/ogse_experiments/master.last_points.long.parquet \
  MASTER_XLSX=analysis/brains/ogse_experiments/master.last_points.xlsx \
    bash signal_analysis/run_dataset.sh brain ogse export_master_xlsx
EOF
            ;;
        rotate)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> rotate

What it does:
  Selects row_kind='signal' rows from master.long.parquet, rotates tensor directions,
  and appends row_kind='signal_rotated' rows. The rotated rows include D_proj.

Variables for this step:
  MASTER_PARQUET      Input/output master table.
  MASTER_SUBJ         Optional selector for the master column 'subj'.
                      Example values in brain data may be BRAIN, LUDG, MBBL.
  MASTER_SHEET        Optional selector for the master column 'sheet'.
                      Example values may be 20220622_BRAIN or 20230619_BRAIN-3.
  DIRS_TXT            Direction table. Default: assets/dirs/dirs_6.txt
  TRA_AXES            Comma-separated axes averaged into direction='tra'.
                      Default: y,z for brain; x,y for phantom.
  LONG_AXES           Comma-separated axes averaged into direction='long'.
                      Default: x for brain; z for phantom.
  ROTATE_SCRIPT       Python script override.
  ROTATE_EXTRA_ARGS   Extra rotate_ogse_tensor.py options.

Useful ROTATE_EXTRA_ARGS:
  --solver lstsq|solve
  --s0_mode dir1|mean
  --b_col bvalue

Examples:
  bash signal_analysis/run_dataset.sh brain ogse rotate

  MASTER_SUBJ=BRAIN \
    bash signal_analysis/run_dataset.sh brain ogse rotate

  MASTER_SHEET=20220622_BRAIN \
    bash signal_analysis/run_dataset.sh brain ogse rotate

  MASTER_SUBJ=BRAIN \
  MASTER_SHEET=20220622_BRAIN \
    bash signal_analysis/run_dataset.sh brain ogse rotate

  ROTATE_EXTRA_ARGS="--s0_mode mean" \
    bash signal_analysis/run_dataset.sh brain ogse rotate
EOF
            ;;
        contrast)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> contrast

What it does:
  Reads declarative contrast selectors from manifests/<type_subj>_<type_seq>/contrasts.csv,
  selects two signal_rotated groups from master, subtracts them point-by-point, and appends
  row_kind='contrast' rows to MASTER_PARQUET. Also produces figures of the two signals and
  the resulting contrast.

  To build resampled contrasts from pre-fitted signal curves instead, use contrast_resampled.

Variables for this step:
  CONTRAST_MANIFEST         CSV contrast manifest.
  MASTER_PARQUET            Input/output master table.
  MAKE_CONTRAST_SCRIPT      Python script override.
  CONTRAST_OUT_ROOT         Output root. Default: $ANALYSIS_ROOT/contrast-data-master
  MAKE_CONTRAST_EXTRA_ARGS  Extra make_contrast.py options (see below).

Manifest columns:
  subj, sheet, roi, direction, td_ms, N_1, N_2, Hz_1, Hz_2

Useful MAKE_CONTRAST_EXTRA_ARGS:
  --apply_grad_corr
      Correct gradient axes by the per-direction grad_correction_factor stored in
      master rows before computing the contrast. Requires step grad_correction first.
      Mutually exclusive with --no_grad_corr.
  --no_grad_corr
      Explicitly skip gradient correction (default behaviour).
  --master-rotated / --no-master-rotated
      Select rotated or raw signal rows from master. Default: --master-rotated.

Examples:
  bash signal_analysis/run_dataset.sh brain ogse contrast

  CONTRAST_MANIFEST=signal_analysis/manifests/brains_ogse/contrasts.csv \
    bash signal_analysis/run_dataset.sh brain ogse contrast

  # Direct contrast with gradient correction
  MAKE_CONTRAST_EXTRA_ARGS="--apply_grad_corr" \
    bash signal_analysis/run_dataset.sh brain ogse contrast
EOF
            ;;
        contrast_resampled)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> contrast_resampled

What it does:
  Builds resampled contrast tables and figures from pre-fitted signal curves. Uses the
  individual signal fit parameters produced by fit_signal or fit_global_signal (step 04)
  to reconstruct each signal curve, then subtracts the two fitted curves to obtain a
  smooth, resampled contrast on a common gradient grid.

  Requires step 04 (fit_signal or fit_global_signal) to have been run first.
  The direct contrast step (contrast) must also have been run so that CONTRAST_ROOT exists.

Variables for this step:
  SIGNAL_FIT_OUT_ROOT         Signal fits directory from step 04. Required — no default.
                              Example: $ANALYSIS_ROOT/fits/ogse_signal_ogse_mixed_offset
  CONTRAST_ROOT               Direct contrast tables root from step 07 (contrast).
                              Default: $ANALYSIS_ROOT/contrast-data-master
  CONTRAST_RESAMPLED_OUT_DIR  Output directory for resampled contrast tables and figures.
                              Default: $ANALYSIS_ROOT/contrast-data-resampled
  EXPORT_RESAMPLED_SCRIPT     Python script override.
  EXPORT_GRID_MIN_MTM         Minimum gradient in mT/m for the common grid. Default: 0
  EXPORT_GRID_MAX_MODE        How to set the upper grid limit: observed_pair_max|fixed.
                              Default: observed_pair_max
  EXPORT_GRID_N               Number of points in the common gradient grid. Default: 1000
  EXPORT_CONTRAST_MODE        signal_fit|rest_only. Default: signal_fit
  EXPORT_MODELS               Optional model filter (space-separated).
  EXPORT_RESAMPLED_EXTRA_ARGS Extra export_ogse_resampled_contrasts_from_fits.py options.

Examples:
  SIGNAL_FIT_OUT_ROOT=analysis/brains/ogse_experiments/fits/ogse_signal_ogse_mixed_offset \
    bash signal_analysis/run_dataset.sh brain ogse contrast_resampled
EOF
            ;;
        plot_signal)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> plot_signal

What it does:
  Plots signal curves selected from master.long.parquet using the script for TYPE_SEQ.

Variables for this step:
  MASTER_PARQUET          Input master table.
  PLOT_SIGNAL_SCRIPT      Python script override.
  PLOT_OUT_ROOT           Output root. Default: $ANALYSIS_ROOT/plots-master/signal
  PLOT_ROW_KIND           signal_rotated|signal. Default: signal_rotated
  PLOT_SUBJ               Optional subj selector.
  PLOT_SHEET              Optional sheet selector.
  PLOT_ROI                Optional ROI selector.
  PLOT_DIRECTION          Optional direction selector, e.g. long|tra|x|y|z.
  PLOT_TD_MS              Optional td_ms selector.
  PLOT_N                  Optional N selector.
  PLOT_SIGNAL_YCOL        value|value_norm. Default: value_norm
  PLOT_SIGNAL_XCOL        Gradient column (x axis).
                          OGSE: g|g_max|g_lin_max|g_thorsten|bvalue|bvalue_g|bvalue_thorsten. Default: g_thorsten
                          NOGSE: g|g_max|g_lin_max|g_thorsten|bvalue|bvalue_g|bvalue_thorsten. Default: g
  PLOT_SIGNAL_G_TYPE      Backward-compatible alias for PLOT_SIGNAL_XCOL.
  PLOT_STAT               avg|std. Default: avg
  PLOT_SIGNAL_EXTRA_ARGS  Extra plot_<type_seq>_signal_vs_g.py options.

Examples:
  bash signal_analysis/run_dataset.sh brain ogse plot_signal

  PLOT_ROI=Left-Lateral-Ventricle PLOT_DIRECTION=long \
    bash signal_analysis/run_dataset.sh brain ogse plot_signal

  PLOT_SUBJ=20220622_BRAIN PLOT_DIRECTION="long" PLOT_SIGNAL_XCOL=g_thorsten \
    bash signal_analysis/run_dataset.sh brain ogse plot_signal
EOF
            ;;
        fit_signal|fit_signal_gradcorr)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> fit_signal
  bash run_dataset.sh <type_subj> <type_seq> fit_signal_gradcorr

What it does:
  Fits signal curves selected from master according to manifests/<type_subj>_<type_seq>/signal_fits.csv.
  Currently this step supports only model=monoexp. Manifest rows may include a
  "model" column, but non-monoexp values are rejected.
  fit_signal_gradcorr is a preset that adds --apply_grad_corr to every fit automatically.

Variables for this step:
  MASTER_PARQUET          Input master table.
  FIT_SIGNAL_SCRIPT       Python script override.
  SIGNAL_FIT_MANIFEST     CSV signal-fit manifest.
  SIGNAL_FIT_OUT_ROOT     Output root. Default: $ANALYSIS_ROOT/fits/signal_fit_<model>_<ycol>_vs_<axis>
  SIGNAL_FIT_MODEL        monoexp only for now. Default: monoexp
  SIGNAL_FIT_B_AXIS       B-value/gradient axis for monoexp fitting.
                          OGSE: g|g_max|g_lin_max|g_thorsten|bvalue|bvalue_g|bvalue_thorsten. Default: bvalue_g
                          NOGSE: g|g_max|g_lin_max|g_thorsten|bvalue|bvalue_g|bvalue_thorsten. Default: g
                          SIGNAL_FIT_G_TYPE is accepted as a backwards-compatible alias.
  SIGNAL_FIT_YCOL         value|value_norm. Default: value_norm
  SIGNAL_FIT_AUTO_ARGS    Point-selection flags. Default: --auto_fit_points
                          Set to "" and use SIGNAL_FIT_EXTRA_ARGS="--fit_points N" for fixed-prefix fits.
  SIGNAL_FIT_EXTRA_ARGS   Extra fit_signal.py options (see below).

Manifest columns:
  subj, sheet, roi, direction, td_ms, N, Hz, model

Useful SIGNAL_FIT_EXTRA_ARGS:

  Gradient correction
  -------------------
  --apply_grad_corr
      Scale the gradient axis by the per-direction grad_correction_factor in master rows
      before fitting. Requires grad_correction to have been run.
      Equivalent to running fit_signal_gradcorr instead of fit_signal.
      Mutually exclusive with --no_grad_corr.
  --no_grad_corr
      Explicitly skip gradient correction (default behaviour).

  Fit point selection  [mutually exclusive]
  -----------------------------------------
  --fit_points N
      Fixed number of leading b_step points to include in the monoexp fit. Default: 6
      if neither --fit_points nor --auto_fit_points is passed.
  --auto_fit_points
      Automatically keep the largest leading b_step prefix whose rmse_log remains
      within tolerance when adding points.
  --auto_fit_tol F        Relative tolerance for --auto_fit_points. Default: 0.05.
  --auto_fit_err_floor F  Absolute RMSE floor before comparing k values. Default: 0.005.
  --auto_fit_min_points N First k tested by --auto_fit_points. Default: 3.
  --auto_fit_max_points N Last k tested by --auto_fit_points. Default: 9.

  Parameter fixing
  ----------------
  --fix_M0 F
      Fix M0 to a constant. Default: 1.0.
  --free_M0
      Fit M0 freely (releases the --fix_M0 default).
  --D0_init F
      Initial D0 seed (mm2/s). Default: 0.0023.

Examples:
  bash signal_analysis/run_dataset.sh brain ogse fit_signal

  SIGNAL_FIT_MANIFEST=signal_analysis/manifests/brains_ogse/signal_fits.csv \
    bash signal_analysis/run_dataset.sh brain ogse fit_signal

  # Monoexp with automatic point selection
  SIGNAL_FIT_EXTRA_ARGS="--auto_fit_points --auto_fit_tol 0.05" \
    bash signal_analysis/run_dataset.sh brain ogse fit_signal

  # Monoexp with gradient correction
  SIGNAL_FIT_EXTRA_ARGS="--apply_grad_corr --D0_init 0.0023" \
    bash signal_analysis/run_dataset.sh brain ogse fit_signal

  bash signal_analysis/run_dataset.sh brain ogse fit_signal_gradcorr
EOF
            ;;
        grad_correction)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> grad_correction

What it does:
  Builds gradient-correction factors from manifests/<type_subj>_<type_seq>/grad_correction.csv.
  For each manifest curve it fits NOGSE-free D0 and monoexp D0, then writes:

    correction_factor = sqrt(D0_nogse / D0_monoexp_avg)

  The monoexp D0 uses auto_fit_points by default in this step.

Variables for this step:
  MASTER_PARQUET             Input/output master table.
  GRAD_CORR_MANIFEST         CSV grad-correction manifest.
  GRAD_CORR_ROI              Reference ROI override. Default: Syringe for brains, Water1 for phantoms.
  GRAD_CORR_OUT_DIR          Output directory. Default: $ANALYSIS_ROOT/fits/grad_correction
  GRAD_CORR_PLOT_DIR         Plot directory. Default: $GRAD_CORR_OUT_DIR/plots
  GRAD_CORR_AUTO_FIT_ARGS    Auto-fit flags for monoexp D0. Default: --auto_fit_points
                             Set to "--no_auto_fit_points" to use all valid monoexp points.
  GRAD_CORR_EXTRA_ARGS       Extra make_grad_correction_table.py options.

Useful GRAD_CORR_EXTRA_ARGS:
  --bbase bvalue_g|bvalue_thorsten|bvalue   Default: bvalue_g
  --gbase g|g_lin_max|g_thorsten            Default: g
  --avg-N 4 8
  --auto_fit_tol 0.05
  --auto_fit_min_points 3
  --auto_fit_max_points 9

Examples:
  bash signal_analysis/run_dataset.sh brain ogse grad_correction

  GRAD_CORR_EXTRA_ARGS="--bbase bvalue_g --avg-N 4 8" \
    bash signal_analysis/run_dataset.sh phantom ogse grad_correction

  GRAD_CORR_AUTO_FIT_ARGS="--no_auto_fit_points" \
    bash signal_analysis/run_dataset.sh brain ogse grad_correction
EOF
            ;;
        fit_contrast|fit_contrast_free|fit_contrast_mixed_global)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> fit_contrast
  bash run_dataset.sh <type_subj> <type_seq> fit_contrast_free
  bash run_dataset.sh <type_subj> <type_seq> fit_contrast_mixed_global

What it does:
  Fits contrast rows from master.long.parquet. Saves per-experiment tables to subdirs and
  accumulates all results in a descriptively named fit_params parquet in FIT_OUT_ROOT.

  fit_contrast_free and fit_contrast_mixed_global are convenience presets. They call
  the same step script as fit_contrast and only set FIT_MODEL when FIT_MODEL is not
  already provided:
    fit_contrast_free         -> ogse_free (OGSE) / nogse_free (NOGSE)
    fit_contrast_mixed_global -> ogse_mixed_global

  Therefore, these are equivalent:
    bash run_dataset.sh brain ogse fit_contrast_free
    FIT_MODEL=ogse_free bash run_dataset.sh brain ogse fit_contrast

  After fitting, run extract_tc_peak (step 09) to consolidate the tc_peak table and
  plot panels, then tc (step 10) for the tc-vs-td analysis.

Variables for this step:
  MASTER_PARQUET        Input master table.
  CONTRAST_PARQUET      Contrast table to fit. Default: $ANALYSIS_ROOT/contrast-data-master/master_contrast.parquet
                        For resampled contrasts (from contrast_resampled step) set this to
                        $ANALYSIS_ROOT/contrast-data-resampled/master_contrast_resampled.parquet
  FIT_CONTRAST_SCRIPT   Python script override.
  FIT_OUT_ROOT          Output root. Default: $ANALYSIS_ROOT/fits/<master_name>/<type_seq>_<ycol>_vs_<gtype>_<model>
  FIT_MODEL             OGSE: ogse_free|ogse_tort|ogse_rest|ogse_rest_offset|ogse_mixed|rest_offset_globC. Default: ogse_free
                        NOGSE: nogse_free|nogse_free_grad_offset|nogse_tort|nogse_rest. Default: nogse_free
                        fit_contrast_free         → presets ogse_free (OGSE) / nogse_free (NOGSE).
                        fit_contrast_mixed_global → presets ogse_mixed_global.
  FIT_GBASE             Gradient column for the x-axis.
                        g|g_lin_max|g_max|g_thorsten|bvalue|bvalue_g|bvalue_thorsten. Default: g_lin_max
  FIT_YCOL              value|value_norm. Default: value_norm
  FIT_STAT              avg|std|ALL. Default: avg
  FIT_EXTRA_ARGS        Extra fit_<type_seq>_contrast_vs_g.py options (see below).

Useful FIT_EXTRA_ARGS (OGSE and NOGSE):

  Gradient correction
  -------------------
  --apply_grad_corr
      Apply the per-direction grad_correction_factor embedded in master contrast rows.
      Only meaningful if the contrast step was run with --apply_grad_corr (direct mode)
      or if the factor columns grad_correction_factor_1/2 are present in the rows.
      Requires grad_correction to have been run first.
      Mutually exclusive with --no_grad_corr.
  --no_grad_corr
      Explicitly skip gradient correction (default behaviour).
  --corr_td_ms F
      Override the td_ms used when applying the gradient correction factor.

  Parameter fixing  [each pair is mutually exclusive]
  ----------------------------------------------------
  --fix_M0 F / --free_M0 [SEED]
      Fix or free M0. Default: free, initial seed 1.0.
  --fix_D0 F / --free_D0 [SEED]
      Fix or free D0 (m²/ms; e.g. 3.2e-12 = 0.0032 mm²/s). Default: free, seed 2.3e-12.
  --fix_tc F / --free_tc [SEED]
      Fix or free tc (ms). Default: free, seed from --tc_init (5.0).
      NOGSE: only applies to model=nogse_rest.
  --tc_init F
      Initial tc seed (ms) when --fix_tc/--free_tc are not given. Default: 5.0.
  --fix_C F / --free_C [SEED]
      Fix or free rest-offset offset parameter C. Default: free, seed 0.0.
      OGSE: applies to model=ogse_rest_offset only.

  NOGSE-specific parameter fixing
  --------------------------------
  --fix_g0 F / --free_g0 [SEED]
      Fix or free gradient offset g0 (mT/m). Default: free, seed 0.0.
      Applies to model=nogse_free_grad_offset only.

  Parameter bounds
  ----------------
  --tc_bounds MIN MAX     Search bounds for tc (ms). Default: 0.1 1000.0.
  --M0_bounds MIN MAX     Search bounds for M0. Default: 0.0 5.0.
  --D0_bounds MIN MAX     Search bounds for D0 (m²/ms). Default: 1e-16 1e-10.
  --C_bounds  MIN MAX     Search bounds for C. Default: 0.0 1.0.
  --g0_bounds MIN MAX     Search bounds for g0 (mT/m). NOGSE. Default: -20.0 20.0.

  Fitted peak
  -----------
  --peak_grid_n N
      Grid resolution for the fitted-peak search. Default: 1000.
  --peak_D0_fix F
      Fixed D0 (m²/ms) used to convert the gradient peak into tc_peak_ms. Default: 3.2e-12.
  --peak_g_max_mTm F
      Upper gradient limit (mT/m) for the peak search.
  --peak_resample_gradient
      Also compute tc_peak_resampled_ms by rebuilding the fitted contrast on a common
      gradient grid. OGSE only.
  --peak_resample_g_max_corr_mTm F
      Upper gradient (mT/m) for the resampled-peak grid (after correction). OGSE only.

  Global / shared-parameter fits  [OGSE only]
  --------------------------------------------
  --global_params PARAM [PARAM...]
      Fit all td curves jointly per subj/ROI/direction and share the listed parameters.
      Unlisted free parameters remain local to each td curve.
      Examples: --global_params tc_ms   or   --global_params tc_ms D0_m2_ms
  --alpha_table PATH
      Fixed-alpha table for model=ogse_mixed_global.
  --alpha_col COL
      Alpha column in --alpha_table. Default: alpha, then alpha_macro.
  --alpha_td_tol_ms F
      Tolerance (ms) for matching alpha rows by td_ms. Default: 0.001.

  Filtering and output
  --------------------
  --directions DIR [DIR...]   Filter by direction (e.g. long tra). Default: all.
  --rois ROI [ROI...]         Filter by ROI. Default: all.
  --n_fit N                   Use only the first N points sorted by x.
  --no_plots                  Skip generating fit plots.
  --oneg                      Allow one-g-per-sequence contrast tables.

Examples:
  bash signal_analysis/run_dataset.sh brain ogse fit_contrast_free

  FIT_MODEL=ogse_free \
    bash signal_analysis/run_dataset.sh brain ogse fit_contrast

  FIT_MODEL=ogse_mixed FIT_GBASE=g_lin_max FIT_YCOL=value_norm \
    bash signal_analysis/run_dataset.sh brain ogse fit_contrast

  # Fit with gradient correction applied during fitting
  FIT_EXTRA_ARGS="--apply_grad_corr" \
    bash signal_analysis/run_dataset.sh brain ogse fit_contrast

  # Fix D0 and restrict tc range
  FIT_EXTRA_ARGS="--fix_D0 3.2e-12 --tc_bounds 0.5 200.0" \
    bash signal_analysis/run_dataset.sh brain ogse fit_contrast_free

  # Jointly fit tc across all td curves
  FIT_EXTRA_ARGS="--global_params tc_ms" \
    bash signal_analysis/run_dataset.sh brain ogse fit_contrast_free
EOF
            ;;
        fit_global_signal)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> fit_global_signal

What it does:
  Fits a global/mixed signal model pooling all subjects from master.long.parquet.
  Parameters that are "global_td" are shared across subjects at the same td_ms.
  Parameters that are "global_contrast" are shared across the N pair defined by the manifest
  (e.g. N=8 and N=4 at the same td_ms/roi/direction), NOT across all N values.
  A manifest CSV controls which curves enter the fit:
    - When global_contrast is active, GLOBAL_SIGNAL_MANIFEST defaults to contrasts.csv
      (columns N_1, N_2) which selects exactly those two N values per row.
    - Otherwise it defaults to signal_fits.csv (column N) for individual curve selection.
  Missing gradient-correction factors default to the mean of other subjects (avg_others).

Available models (set via GLOBAL_SIGNAL_MODEL):
  OGSE: ogse_free | ogse_rest | ogse_rest_offset | ogse_mixed | ogse_mixed_offset
  NOGSE: nogse_free | nogse_free_offset | nogse_rest | nogse_rest_offset | nogse_mixed | nogse_mixed_offset
  Parameter coverage per model:
    ogse_free / nogse_free*          →  M0, D0
    ogse_rest / nogse_rest           →  tc, M0, D0
    ogse_rest_offset                 →  tc, M0, D0, C
    ogse_mixed / nogse_mixed         →  tc, alpha, M0, D0
    ogse_mixed_offset / nogse_mixed_offset → tc, alpha, RN, M0, D0, C
  Mode variables for parameters not in the chosen model are silently ignored.

Variables for this step:
  MASTER_PARQUET                 Input master table.
  FIT_GLOBAL_SIGNAL_SCRIPT       Python script override.
  GLOBAL_SIGNAL_OUT_ROOT         Output root.
                                 Default: $ANALYSIS_ROOT/fits/<type_seq>_signal_<model>
                                 Override explicitly when running multiple models or correction modes
                                 so outputs don't overwrite each other.
  GLOBAL_SIGNAL_ROW_KIND         signal_rotated|signal. Default: signal_rotated
  GLOBAL_SIGNAL_MODEL            See model list above.
                                 Default: ogse_mixed_offset (OGSE) | nogse_mixed_offset (NOGSE)
  GLOBAL_SIGNAL_YCOL             value|value_norm. Default: value
  GLOBAL_SIGNAL_G_TYPE           Gradient column.
                                 OGSE brain: g_thorsten (default) | g_lin_max | bvalue_thorsten | ...
                                 OGSE phantom / NOGSE: g (default) | g_lin_max | g_thorsten | ...
  GLOBAL_SIGNAL_STAT             avg|std. Default: avg
  GLOBAL_SIGNAL_MIN_POINTS       Minimum points per group to attempt a fit. Default: 4
  Parameter modes and fixed values
  ---------------------------------
  Each parameter has a _MODE variable and a _FIXED variable. They work as a pair:
    _MODE controls how the parameter is treated across subjects and td_ms values.
    _FIXED sets the numeric value when _MODE=fixed.
  Mode and fixed vars for parameters not present in the chosen model are silently ignored.

  GLOBAL_SIGNAL_TC_MODE          fixed|free|global_td|global_contrast. Default: global_td
  GLOBAL_SIGNAL_TC_FIXED         Fixed tc value [ms]. Required when TC_MODE=fixed.
                                 (tc: rest, rest_offset, mixed, mixed_offset only)
  GLOBAL_SIGNAL_ALPHA_MODE       fixed|free|global_td|global_contrast. Default: global_td
  GLOBAL_SIGNAL_ALPHA_FIXED      Fixed alpha value [0–1]. Required when ALPHA_MODE=fixed.
                                 (alpha: mixed, mixed_offset only)
  GLOBAL_SIGNAL_RN_MODE          fixed|free|global_td|global_contrast. Default: global_td
  GLOBAL_SIGNAL_RN_FIXED         Fixed RN value. Required when RN_MODE=fixed.
                                 (RN: mixed_offset, nogse_mixed_offset only)
  GLOBAL_SIGNAL_M0_MODE          fixed|free|global_td|global_contrast. Default: global_contrast
  GLOBAL_SIGNAL_M0_FIXED         Fixed M0 value. Required when M0_MODE=fixed.
  GLOBAL_SIGNAL_C_MODE           fixed|free|global_td|global_contrast. Default: global_contrast
  GLOBAL_SIGNAL_C_FIXED          Fixed C value. Required when C_MODE=fixed.
                                 (C: rest_offset, mixed_offset only)
  GLOBAL_SIGNAL_D0_MODE          fixed|free|global_td|global_contrast. Default: fixed
  GLOBAL_SIGNAL_D0_FIXED         Fixed D0 in m²/s. Default: 3.2e-12 (brain) | 2.3e-12 (phantom)
  GLOBAL_SIGNAL_DIRECTIONS       Space- or comma-separated directions, or ALL.
                                 Default: long tra (OGSE brain) | ALL (others)
  GLOBAL_SIGNAL_ROIS             Space- or comma-separated ROIs, or ALL. Default: ALL
  GLOBAL_SIGNAL_SUBJS            Space- or comma-separated subjects, or ALL. Default: ALL
  GLOBAL_SIGNAL_APPLY_GRAD_CORR  true|false. Default: true
                                 Missing factors are handled at the grad_correction step
                                 (cross-subject fill). By the time fit_global_signal runs,
                                 all signal rows should have a factor or grad_correction
                                 was not run.
  Manifest (curve selection)
  --------------------------
  GLOBAL_SIGNAL_MANIFEST         Path to a manifest CSV that selects which curves to fit.
                                 Auto-detected when not set:
                                   · global_contrast mode → $MANIFEST_DIR/contrasts.csv
                                   · otherwise            → $MANIFEST_DIR/signal_fits.csv
                                 If the file does not exist the step warns and fits all rows.
                                 Set to "none" to always fit all rows regardless.
                                 Contrast manifest (columns N_1, N_2): selects exactly those two
                                   N values per row, linking them as a pair for global_contrast.
                                 Signal manifest (column N): selects individual N values.
                                 Values of ALL in any column are wildcards.
  GLOBAL_SIGNAL_EXTRA_ARGS       Extra fit_global_signal.py options (e.g. --tc_init 5.0).

Examples:
  bash signal_analysis/run_dataset.sh brain ogse fit_global_signal

  # ogse_mixed_offset with RN fixed at 10 (grad correction on by default):
  GLOBAL_SIGNAL_RN_MODE=fixed \
  GLOBAL_SIGNAL_RN_FIXED=10 \
    bash signal_analysis/run_dataset.sh brain ogse fit_global_signal

  # ogse_rest (no alpha, no RN, no C): only TC and M0 modes apply
  GLOBAL_SIGNAL_MODEL=ogse_rest \
  GLOBAL_SIGNAL_TC_MODE=global_td \
    bash signal_analysis/run_dataset.sh brain ogse fit_global_signal

  # run without grad correction
  GLOBAL_SIGNAL_APPLY_GRAD_CORR=false \
  GLOBAL_SIGNAL_OUT_ROOT="analysis/brains/ogse_experiments/fits/ogse_signal_ogse_mixed_offset_raw" \
    bash signal_analysis/run_dataset.sh brain ogse fit_global_signal

  # use a custom manifest (e.g. only N=8/N=4 pairs at td=90 for specific ROIs)
  GLOBAL_SIGNAL_MANIFEST="my_manifests/subset_contrasts.csv" \
    bash signal_analysis/run_dataset.sh brain ogse fit_global_signal

  # fit all signal rows without manifest filtering
  GLOBAL_SIGNAL_MANIFEST=none \
    bash signal_analysis/run_dataset.sh brain ogse fit_global_signal

  # restrict to specific ROIs/directions
  GLOBAL_SIGNAL_DIRECTIONS="long tra" GLOBAL_SIGNAL_ROIS="Left-Lateral-Ventricle AntCC" \
    bash signal_analysis/run_dataset.sh brain ogse fit_global_signal
EOF
            ;;
        plot_monoexp_d)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> plot_monoexp_d

What it does:
  Builds monoexp D vs td_ms and Delta_app_ms plots from signal fit outputs.

Variables for this step:
  PLOT_MONOEXP_D_SCRIPT     Python script override.
  SIGNAL_FITS_ROOT          Root scanned for signal fit parquets.
                            Default: $ANALYSIS_ROOT/fits (scans all subdirs).
                            Set explicitly to the specific fit folder, e.g.:
                            $ANALYSIS_ROOT/fits/master/ogse_value_norm_vs_bvaluethorsten_monoexp
  MONOEXP_D_OUT_DIR         Output plot/table directory. Default: $ANALYSIS_ROOT/plots-master/monoexp_D_vs_time
  PLOT_MONOEXP_D_EXTRA_ARGS Extra plot_monoexp_D_vs_time.py options.

Examples:
  SIGNAL_FITS_ROOT=analysis/brains/ogse_experiments/fits/master/ogse_value_norm_vs_bvaluethorsten_monoexp \
    bash signal_analysis/run_dataset.sh brain ogse plot_monoexp_d
EOF
            ;;
        alpha)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> alpha

What it does:
  1. Computes alpha_macro from D_proj values in signal_rotated rows and writes
     summary_alpha_values.xlsx (make_alpha_macro_summary.py).
  2. Adds alpha_macro back into matching MASTER_PARQUET rows by subj/roi/direction.
  3. Generates per-group D vs Delta_app_ms plots with a dotted alpha annotation
     at the selected bvalue (plot_D0_vs_Delta.py). The bvalue per ROI is read
     automatically from the summary written in step 1, so --roi-bvalmax is
     respected in the plots without extra configuration.

Variables for this step:
  MASTER_PARQUET       Input master table. Must contain signal_rotated rows with D_proj.
  ALPHA_MACRO_SCRIPT   Python script override for make_alpha_macro_summary.py.
  PLOT_D0_SCRIPT       Python script override for plot_D0_vs_Delta.py.
  ALPHA_N              N selector passed to make_alpha_macro_summary.py. Default: 1
  ALPHA_OUT_DIR        Output directory for all outputs. Default: $ANALYSIS_ROOT/alpha_macro/master
  ALPHA_PLOT_BSTEPS    Space-separated bvalue positions to use as candidates for
                       alpha selection and to draw in D-vs-Delta plots.
  ALPHA_PLOT_BVALUES   Space-separated rounded bvalues to use as candidates for
                       alpha selection and to draw in D-vs-Delta plots.
  ALPHA_EXTRA_ARGS     Extra options forwarded to make_alpha_macro_summary.py (see below).
  DPROJ_N              N selector for the D-vs-Delta plots (if different from ALPHA_N).
  DPROJ_HZ             Hz selector for the D-vs-Delta plots.
  DPROJ_ROIS           Space-separated ROI list for the D-vs-Delta plots.
  DPROJ_DIRS           Space-separated direction list for the D-vs-Delta plots.
                       Example: DPROJ_DIRS="long tra"
  PLOT_D0_EXTRA_ARGS   Extra options forwarded only to plot_D0_vs_Delta.py.

Note:
  This step always passes --no-master-fit-params, so alpha_macro results are NOT
  appended to any cumulative fit-params table. Instead, the scalar alpha_macro column
  is written directly onto matching master rows; matching rows with the same
  subj/roi/direction receive repeated alpha_macro values as expected.

Useful ALPHA_EXTRA_ARGS (passed to make_alpha_macro_summary.py):
  --bvalmax X              Use bstep X among the candidate bvalues, or use bvalue X
                           directly when X is outside the candidate bstep range.
                           Default: highest candidate bvalue.
  --roi-bvalmax ROI=X      Per-ROI override. X can be a candidate bstep or a bvalue.
                           Example: --roi-bvalmax CSF=3 or --roi-bvalmax CSF=1280
  --dirs DIR1 DIR2         Restrict summary to specific directions.
  --rois ROI1 ROI2         Restrict summary to specific ROIs.
  --subjs S1 S2            Restrict summary to specific subj values, e.g. MBBL LUDG.
  --sheets SHEET1 SHEET2   Restrict summary to specific sheet/session names,
                           e.g. 20230630_MBBL-3 20230710_LUDG-3.
  --reference-D0 F         Reference D0 used to compute alpha_macro. Default: 0.0032
  --no-annotate-master-alpha
                           Do not write alpha_macro back to MASTER_PARQUET.
                           Useful for testing subset commands.
  --out-plot PATH          Override output path for alpha_macro_vs_roi.png.

Useful PLOT_D0_EXTRA_ARGS (passed only to plot_D0_vs_Delta.py):
  --bvalmax N              Fallback bstep for groups absent in summary_alpha_values.xlsx.
  --plot-bsteps N1 N2      Plot only these 1-based bvalue positions per group.
                           Example: --plot-bsteps 1 3 5
  --plot-bvalues B1 B2     Plot only these rounded bvalues per group.
                           Example: --plot-bvalues 300 600
  --reference-D0 F         Reference D0 used for the horizontal annotation. Default: 0.0032

Outputs:
  $MASTER_PARQUET                                updated with alpha_macro column
  $ALPHA_OUT_DIR/summary_alpha_values.xlsx      alpha_macro per subj/roi/direction
  $ALPHA_OUT_DIR/D_vs_delta_app.combined.xlsx   aggregated D vs Delta_app table
  $ALPHA_OUT_DIR/alpha_macro_vs_roi.png         bar summary plot
  $ALPHA_OUT_DIR/<subj>/<roi>/<dir>_*.png       per-group D vs Delta_app curves

Examples:
  bash signal_analysis/run_dataset.sh brain ogse alpha

  ALPHA_N=1 \
  ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax AntCC=7 --roi-bvalmax CSF=3 --dirs long tra" \
  DPROJ_DIRS="long tra" \
    bash signal_analysis/run_dataset.sh brain ogse alpha

  # Restrict plots to a subset of ROIs while keeping the full summary
  ALPHA_EXTRA_ARGS="--bvalmax 5 --roi-bvalmax Syringe=7" \
  DPROJ_DIRS="long tra" \
  DPROJ_ROIS="AntCC MidCC" \
    bash signal_analysis/run_dataset.sh brain ogse alpha
EOF
            ;;
        extract_tc_peak)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> extract_tc_peak

What it does:
  1. Consolidates fit_params from step 08 (fit_contrast) into a canonical tc_peak table
     using run_tc_pipeline.py, writing tc_peak_table.parquet and tc_peak_table.xlsx.
  2. (OGSE only) Plots tc_peak panels using plot_ogse-contrast_tc_peak_panels.py.

  The tc_peak_table.parquet produced here is the default input for step tc (analyze_tc_vs_td).

Variables for this step:
  FIT_OUT_ROOT              Contrast fits directory from step 08. Derived automatically
                            using the same FIT_MODEL/FIT_GBASE/FIT_YCOL defaults.
  TC_PIPELINE_SCRIPT        Python script override for run_tc_pipeline.py.
  TC_PEAK_DIR               Output directory. Default: $FIT_OUT_ROOT/tc_peak
  TC_PEAK_MODELS            Optional model filter (space-separated).
  TC_PEAK_SUBJS             Optional subject filter (space-separated).
  TC_PEAK_ROIS              Optional ROI filter (space-separated).
  TC_PEAK_DIRECTIONS        Optional direction filter (space-separated).
  TC_PIPELINE_EXTRA_ARGS    Extra run_tc_pipeline.py options.

  OGSE panel plot variables:
  PLOT_TC_PEAKS_SCRIPT      Python script override for plot_ogse-contrast_tc_peak_panels.py.
  TC_PEAKS_OUT_DIR          Output directory for panels. Default: $TC_PEAK_DIR/panels
  TC_PEAKS_CONTRAST_ROOT    Root with resampled contrast tables (from contrast_resampled).
                            Default: $ANALYSIS_ROOT/contrast-data-resampled
  TC_PEAKS_CONTRAST_SOURCE  direct|fitted_resampled|auto. Default: fitted_resampled
  TC_PEAKS_PEAK_SOURCE      standard|resampled|both. Default: resampled
  TC_PEAKS_MODELS           Optional model filter for panels.
  TC_PEAKS_SUBJS            Optional subject filter for panels.
  TC_PEAKS_ROIS             Optional ROI filter for panels.
  TC_PEAKS_DIRECTIONS       Optional direction filter for panels.
  TC_PEAKS_N1               Keep only fits with this first OGSE N value.
  TC_PEAKS_N2               Keep only fits with this second OGSE N value.
  TC_PEAKS_X_VARS           X-axis variables for panels. Default: g Ld lcf lcf_a tc
  TC_PEAKS_SHOW_RESAMPLED_FIT Set to 1 to overlay the resampled fit curve.
  TC_PEAKS_HIDE_DATA_POINTS   Set to 1 to hide experimental contrast points.
  TC_PEAKS_EXCLUDE_FITRESAMP  Set to 1 to skip fit_params from fitted/resampled contrasts.
  TC_PEAKS_EXTRA_ARGS       Extra plot_ogse-contrast_tc_peak_panels.py options.

Examples:
  bash signal_analysis/run_dataset.sh brain ogse extract_tc_peak

  # Filter to specific models/ROIs
  TC_PEAK_MODELS="ogse_free ogse_tort" \
  TC_PEAKS_ROIS="AntCC MidCC" \
    bash signal_analysis/run_dataset.sh brain ogse extract_tc_peak
EOF
            ;;
        tc)
            cat <<'EOF'
Usage:
  bash run_dataset.sh <type_subj> <type_seq> tc

What it does:
  Fits tc-vs-td summaries from the tc_peak table produced by extract_tc_peak.

Variables for this step:
  TC_FIT_PARAMS      Tc_peak parquet from step 09 (extract_tc_peak). Derived automatically
                     from TC_PEAK_DIR/tc_peak_table.parquet using the same defaults.
                     Override only if pointing to a non-default location.
  TC_VS_TD_SCRIPT    Python script override.
  TC_OUT_DIR         Output root. Default: $ANALYSIS_ROOT/fits/tc_vs_td_master
                     Each method writes to its own subdirectory: $TC_OUT_DIR/$TC_METHOD
  TC_METHOD          pseudohuber_free|pseudohuber_fixed_macro|linear. Default: pseudohuber_fixed_macro
                     pseudohuber_free: fits c, delta, alpha_macro freely.
                     pseudohuber_fixed_macro: fits c and delta; alpha_macro fixed from summary.
                       Requires: TC_EXTRA_ARGS="--summary-alpha path/to/summary_alpha_values.xlsx"
                     linear: fits tc(Td) = c + slope * Td. Useful for the large-Td regime.
  TC_Y_COL           tc_peak_ms|tc_peak_resampled_ms. Default: tc_peak_ms
                     Override TC_OUT_DIR when changing TC_Y_COL to avoid output conflicts.
  TC_EXTRA_ARGS      Extra run_tc_vs_td.py options.

Examples:
  # pseudohuber with fixed alpha_macro from summary
  TC_FIT_PARAMS="fits/master/ogse_.../fit_params.*.parquet" \
  TC_METHOD=pseudohuber_fixed_macro \
  TC_EXTRA_ARGS="--summary-alpha analysis/brains/ogse_experiments/alpha_macro/master/summary_alpha_values.xlsx" \
    bash signal_analysis/run_dataset.sh brain ogse tc

  # linear fit on resampled tc_peak
  TC_FIT_PARAMS="fits/master/ogse_.../fit_params.*.parquet" \
  TC_METHOD=linear \
  TC_Y_COL=tc_peak_resampled_ms \
  TC_OUT_DIR="analysis/brains/ogse_experiments/fits/tc_vs_td_resampled_master" \
    bash signal_analysis/run_dataset.sh brain ogse tc
EOF
            ;;
        *)
            echo "ERROR: no help available for unknown step '$step'" >&2
            return 2
            ;;
    esac
}

pipeline_maybe_step_help() {
    local step="$1"
    shift || true
    if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        pipeline_step_help "$step"
        exit 0
    fi
}

pipeline_step_script() {
    case "$1" in
        ingest) echo "$TEMPLATE_ROOT/steps/01_ingest_results.sh" ;;
        filter_master_points) echo "$TEMPLATE_ROOT/steps/00_filter_master_points.sh" ;;
        export_master_xlsx) echo "$TEMPLATE_ROOT/steps/99_export_master_xlsx.sh" ;;
        rotate) echo "$TEMPLATE_ROOT/steps/02_rotate_signals.sh" ;;
        plot_signal) echo "$TEMPLATE_ROOT/steps/03_plot_signals.sh" ;;
        fit_signal|fit_signal_gradcorr) echo "$TEMPLATE_ROOT/steps/04_fit_signals.sh" ;;
        fit_global_signal) echo "$TEMPLATE_ROOT/steps/04_fit_signals_global.sh" ;;
        grad_correction) echo "$TEMPLATE_ROOT/steps/05_make_grad_correction.sh" ;;
        alpha) echo "$TEMPLATE_ROOT/steps/06_alpha_macro.sh" ;;
        plot_monoexp_d) echo "$TEMPLATE_ROOT/steps/06b_plot_monoexp_D_vs_time.sh" ;;
        contrast) echo "$TEMPLATE_ROOT/steps/07_make_contrast.sh" ;;
        contrast_resampled) echo "$TEMPLATE_ROOT/steps/07b_make_contrast_resampled.sh" ;;
        fit_contrast|fit_contrast_free|fit_contrast_mixed_global) echo "$TEMPLATE_ROOT/steps/08_fit_contrast.sh" ;;
        extract_tc_peak) echo "$TEMPLATE_ROOT/steps/09_extract_tc_peak.sh" ;;
        tc) echo "$TEMPLATE_ROOT/steps/10_analyze_tc_vs_td.sh" ;;
        *) return 1 ;;
    esac
}

pipeline_prepare_step_env() {
    case "$1" in
        fit_signal_gradcorr)
            SIGNAL_FIT_MODEL="${SIGNAL_FIT_MODEL:-monoexp}"
            if [[ "$TYPE_SEQ" == "nogse" ]]; then
                SIGNAL_FIT_G_TYPE="${SIGNAL_FIT_G_TYPE:-g}"
            else
                SIGNAL_FIT_G_TYPE="${SIGNAL_FIT_G_TYPE:-bvalue_g}"
                if [[ " ${SIGNAL_FIT_EXTRA_ARGS:-} " != *" --directions "* && " ${SIGNAL_FIT_EXTRA_ARGS:-} " != *" --direction "* ]]; then
                    SIGNAL_FIT_EXTRA_ARGS="${SIGNAL_FIT_EXTRA_ARGS:-} --directions long tra"
                fi
            fi
            SIGNAL_FIT_EXTRA_ARGS="${SIGNAL_FIT_EXTRA_ARGS:-} --apply_grad_corr"
            export SIGNAL_FIT_MODEL SIGNAL_FIT_G_TYPE SIGNAL_FIT_EXTRA_ARGS
            ;;
        fit_contrast_free)
            if [[ "$TYPE_SEQ" == "nogse" ]]; then
                FIT_MODEL="${FIT_MODEL:-nogse_free}"
            else
                FIT_MODEL="${FIT_MODEL:-ogse_free}"
            fi
            export FIT_MODEL
            ;;
        fit_contrast_mixed_global)
            FIT_MODEL="${FIT_MODEL:-ogse_mixed_global}"
            export FIT_MODEL
            ;;
    esac
}

pipeline_run_steps() {
    if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        pipeline_usage
        return 0
    fi
    if [[ $# -eq 2 && ( "${2:-}" == "-h" || "${2:-}" == "--help" ) ]]; then
        pipeline_step_help "$1"
        return 0
    fi
    local raw_steps="${PIPELINE_STEPS:-$*}"
    raw_steps="${raw_steps//,/ }"
    if [[ -z "${raw_steps// }" ]]; then
        pipeline_usage
        return 0
    fi

    local step script
    for step in $raw_steps; do
        script="$(pipeline_step_script "$step")" || {
            echo "ERROR: unknown step '$step'" >&2
            pipeline_usage >&2
            exit 2
        }
        if [[ "$step" != "ingest" && "$step" != "filter_master_points" && "$step" != "export_master_xlsx" ]]; then
            pipeline_apply_master_last_points
        fi
        pipeline_prepare_step_env "$step"
        echo
        echo "==> [$TYPE_SUBJ/$TYPE_SEQ] $step"
        DATASET="$DATASET" TYPE_SUBJ="$TYPE_SUBJ" TYPE_SEQ="$TYPE_SEQ" bash "$script"
    done
}
