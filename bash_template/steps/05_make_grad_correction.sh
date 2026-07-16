#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help grad_correction "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}"

GRAD_CORR_SCRIPT="${GRAD_CORR_SCRIPT:-$REPO_ROOT/scripts/data/make_grad_correction_table.py}"
GRAD_CORR_MANIFEST="${GRAD_CORR_MANIFEST:-$MANIFEST_DIR/grad_correction.csv}"
GRAD_CORR_OUT_DIR="${GRAD_CORR_OUT_DIR:-$ANALYSIS_ROOT/fits/grad_correction}"
GRAD_CORR_PLOT_DIR="${GRAD_CORR_PLOT_DIR:-$GRAD_CORR_OUT_DIR/plots}"

case "$DATASET" in
    brains)   GRAD_CORR_ROI="${GRAD_CORR_ROI:-Syringe}" ;;
    phantoms) GRAD_CORR_ROI="${GRAD_CORR_ROI:-Water1}" ;;
esac

pipeline_require_file "$GRAD_CORR_SCRIPT" "gradient correction script"
pipeline_require_file "$GRAD_CORR_MANIFEST" "gradient correction manifest"
pipeline_require_file "$MASTER_PARQUET" "master table"
mkdir -p "$GRAD_CORR_OUT_DIR"

"$PY" "$GRAD_CORR_SCRIPT" \
    --manifest "$GRAD_CORR_MANIFEST" \
    --master-parquet "$MASTER_PARQUET" \
    --out-xlsx "$GRAD_CORR_OUT_DIR/grad_correction.xlsx" \
    --out-csv "$GRAD_CORR_OUT_DIR/grad_correction.csv" \
    --plot-dir "$GRAD_CORR_PLOT_DIR" \
    --roi "$GRAD_CORR_ROI" \
    ${GRAD_CORR_EXTRA_ARGS:-}
