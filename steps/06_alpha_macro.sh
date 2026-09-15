#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help alpha "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}"

ALPHA_MACRO_SCRIPT="${ALPHA_MACRO_SCRIPT:-$REPO_ROOT/scripts/summary/make_alpha_macro_summary.py}"
PLOT_D0_SCRIPT="${PLOT_D0_SCRIPT:-$REPO_ROOT/scripts/plotting/plot_D0_vs_Delta.py}"
ALPHA_OUT_DIR="${ALPHA_OUT_DIR:-$ANALYSIS_ROOT/alpha_macro/master}"
SUMMARY_ALPHA="$ALPHA_OUT_DIR/summary_alpha_values.xlsx"

pipeline_require_file "$ALPHA_MACRO_SCRIPT" "alpha macro script"
pipeline_require_file "$MASTER_PARQUET" "master table"
mkdir -p "$ALPHA_OUT_DIR"

alpha_plot_args=()
[[ -n "${ALPHA_PLOT_BSTEPS:-}" ]]  && alpha_plot_args+=(--plot-bsteps ${ALPHA_PLOT_BSTEPS})
[[ -n "${ALPHA_PLOT_BVALUES:-}" ]] && alpha_plot_args+=(--plot-bvalues ${ALPHA_PLOT_BVALUES})

"$PY" "$ALPHA_MACRO_SCRIPT" \
    --master-parquet "$MASTER_PARQUET" \
    --no-master-fit-params \
    --N "${ALPHA_N:-1}" \
    --reference-D0 "${ALPHA_REFERENCE_D0_MM2_S:-0.0032}" \
    --reference-D0-error "${ALPHA_REFERENCE_D0_ERROR_MM2_S:-0.0000283512}" \
    --out-summary "$SUMMARY_ALPHA" \
    --out-avg "$ALPHA_OUT_DIR/D_vs_delta_app.combined.xlsx" \
    "${alpha_plot_args[@]}" \
    ${ALPHA_EXTRA_ARGS:-}

# Generate D vs Delta_app plots with alpha annotation per ROI.
# Reads selected_bstep_map from the summary just written so each ROI gets
# the correct bvalue annotated even when --roi-bvalmax was used above.
if [[ -f "$PLOT_D0_SCRIPT" ]]; then
    plot_args=(
        --master-parquet "$MASTER_PARQUET"
        --out-dir "$ALPHA_OUT_DIR"
        --summary-alpha "$SUMMARY_ALPHA"
        --reference-D0 "${ALPHA_REFERENCE_D0_MM2_S:-0.0032}"
        --reference-D0-error "${ALPHA_REFERENCE_D0_ERROR_MM2_S:-0.0000283512}"
    )
    [[ -n "${DPROJ_N:-}" ]]    && plot_args+=(--N "$DPROJ_N")
    [[ -n "${DPROJ_HZ:-}" ]]   && plot_args+=(--Hz "$DPROJ_HZ")
    [[ -n "${DPROJ_ROIS:-}" ]] && plot_args+=(--rois ${DPROJ_ROIS})
    [[ -n "${DPROJ_DIRS:-}" ]] && plot_args+=(--dirs ${DPROJ_DIRS})
    [[ -n "${ALPHA_PLOT_BSTEPS:-}" ]]  && plot_args+=(--plot-bsteps ${ALPHA_PLOT_BSTEPS})
    [[ -n "${ALPHA_PLOT_BVALUES:-}" ]] && plot_args+=(--plot-bvalues ${ALPHA_PLOT_BVALUES})
    "$PY" "$PLOT_D0_SCRIPT" "${plot_args[@]}" ${PLOT_D0_EXTRA_ARGS:-}
fi
