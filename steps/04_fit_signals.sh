#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help fit_signal "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}" "${TYPE_SEQ:-ogse}"

FIT_SIGNAL_SCRIPT="${FIT_SIGNAL_SCRIPT:-$REPO_ROOT/scripts/fitting/fit_signal.py}"
SIGNAL_FIT_MANIFEST="${SIGNAL_FIT_MANIFEST:-$MANIFEST_DIR/signal_fits.csv}"
SIGNAL_FIT_MODEL="${SIGNAL_FIT_MODEL:-monoexp}"
SIGNAL_FIT_YCOL="${SIGNAL_FIT_YCOL:-value_norm}"

if [[ "$SIGNAL_FIT_MODEL" == "ogse_free" ]]; then
    SIGNAL_FIT_B_AXIS="${SIGNAL_FIT_B_AXIS:-${SIGNAL_FIT_G_TYPE:-g}}"
else
    case "$TYPE_SEQ" in
        ogse)  SIGNAL_FIT_B_AXIS="${SIGNAL_FIT_B_AXIS:-${SIGNAL_FIT_G_TYPE:-bvalue_g}}" ;;
        nogse) SIGNAL_FIT_B_AXIS="${SIGNAL_FIT_B_AXIS:-${SIGNAL_FIT_G_TYPE:-g}}" ;;
    esac
fi

SIGNAL_FIT_OUT_ROOT="${SIGNAL_FIT_OUT_ROOT:-$ANALYSIS_ROOT/fits/signal_fit_${SIGNAL_FIT_MODEL}_${SIGNAL_FIT_YCOL}_vs_${SIGNAL_FIT_B_AXIS}}"
SIGNAL_FIT_PLOT_DIR="${SIGNAL_FIT_PLOT_DIR:-$SIGNAL_FIT_OUT_ROOT/plots}"
SIGNAL_FIT_AUTO_ARGS="${SIGNAL_FIT_AUTO_ARGS:---auto_fit_points}"

pipeline_require_file "$FIT_SIGNAL_SCRIPT" "signal fit script"
pipeline_require_file "$MASTER_PARQUET" "master table"
mkdir -p "$SIGNAL_FIT_OUT_ROOT"

manifest_args=()
if [[ -f "$SIGNAL_FIT_MANIFEST" ]]; then
    manifest_args=(--manifest "$SIGNAL_FIT_MANIFEST")
else
    echo "WARNING: signal-fit manifest not found: $SIGNAL_FIT_MANIFEST"
    echo "         Fitting selected master rows directly. Use MASTER_* selectors or SIGNAL_FIT_EXTRA_ARGS to narrow the run."
fi

correction_args=()
case "${SIGNAL_FIT_CORRECTION_MODE:-}" in
    raw)       correction_args=(--no-grad-corr) ;;
    corrected) correction_args=(--apply-grad-corr) ;;
    "") ;;
    *)
        echo "ERROR: SIGNAL_FIT_CORRECTION_MODE must be raw or corrected." >&2
        exit 2
        ;;
esac

"$PY" "$FIT_SIGNAL_SCRIPT" \
    --master-parquet "$MASTER_PARQUET" \
    --row-kind "${SIGNAL_FIT_ROW_KIND:-signal_rotated}" \
    "${manifest_args[@]}" \
    --model "$SIGNAL_FIT_MODEL" \
    --out-root "$SIGNAL_FIT_OUT_ROOT" \
    --plot-dir "$SIGNAL_FIT_PLOT_DIR" \
    --ycol "$SIGNAL_FIT_YCOL" \
    --b-axis "$SIGNAL_FIT_B_AXIS" \
    "${correction_args[@]}" \
    ${SIGNAL_FIT_AUTO_ARGS:-} \
    ${SIGNAL_FIT_EXTRA_ARGS:-}
