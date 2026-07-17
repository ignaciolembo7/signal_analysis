#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help plot_signal "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}"

if [[ "$TYPE_SEQ" == "nogse" ]]; then
    PLOT_SIGNAL_SCRIPT="${PLOT_SIGNAL_SCRIPT:-$REPO_ROOT/scripts/plotting/plot_nogse_signal_vs_g.py}"
    DEFAULT_PLOT_SIGNAL_XCOL="${DEFAULT_PLOT_SIGNAL_XCOL:-g}"
else
    PLOT_SIGNAL_SCRIPT="${PLOT_SIGNAL_SCRIPT:-$REPO_ROOT/scripts/plotting/plot_ogse_signal_vs_g.py}"
    DEFAULT_PLOT_SIGNAL_XCOL="${DEFAULT_PLOT_SIGNAL_XCOL:-g_thorsten}"
fi
PLOT_OUT_ROOT="${PLOT_OUT_ROOT:-$ANALYSIS_ROOT/plots-master/signal}"

pipeline_require_file "$PLOT_SIGNAL_SCRIPT" "plot signal script"
pipeline_require_file "$MASTER_PARQUET" "master table"
mkdir -p "$PLOT_OUT_ROOT"

args=(--master-parquet "$MASTER_PARQUET" --out_root "$PLOT_OUT_ROOT" --row-kind "${PLOT_ROW_KIND:-signal_rotated}")
[[ "${PLOT_SUBJ:-ALL}" != "ALL" ]] && args+=(--subj "$PLOT_SUBJ")
[[ "${PLOT_SHEET:-ALL}" != "ALL" ]] && args+=(--sheet "$PLOT_SHEET")
[[ "${PLOT_ROI:-ALL}" != "ALL" ]] && args+=(--roi "$PLOT_ROI")
[[ "${PLOT_DIRECTION:-ALL}" != "ALL" ]] && args+=(--direction "$PLOT_DIRECTION")
[[ -n "${PLOT_TD_MS:-}" ]] && args+=(--td_ms "$PLOT_TD_MS")
[[ -n "${PLOT_N:-}" ]] && args+=(--N "$PLOT_N")

"$PY" "$PLOT_SIGNAL_SCRIPT" \
    "${args[@]}" \
    --ycol "${PLOT_SIGNAL_YCOL:-value_norm}" \
    --xcol "${PLOT_SIGNAL_XCOL:-${PLOT_SIGNAL_G_TYPE:-$DEFAULT_PLOT_SIGNAL_XCOL}}" \
    --stat "${PLOT_STAT:-avg}" \
    ${PLOT_SIGNAL_EXTRA_ARGS:-}
