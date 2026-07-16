#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help plot_monoexp_d "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}"

PLOT_MONOEXP_D_SCRIPT="${PLOT_MONOEXP_D_SCRIPT:-$REPO_ROOT/scripts/plotting/plot_monoexp_D_vs_time.py}"
SIGNAL_FITS_ROOT="${SIGNAL_FITS_ROOT:-$ANALYSIS_ROOT/fits}"
MONOEXP_D_OUT_DIR="${MONOEXP_D_OUT_DIR:-$ANALYSIS_ROOT/plots-master/monoexp_D_vs_time}"

pipeline_require_file "$PLOT_MONOEXP_D_SCRIPT" "monoexp D plot script"
pipeline_require_file "$SIGNAL_FITS_ROOT" "signal fits root"
mkdir -p "$MONOEXP_D_OUT_DIR"

"$PY" "$PLOT_MONOEXP_D_SCRIPT" \
    --fits-root "$SIGNAL_FITS_ROOT" \
    --out-dir "$MONOEXP_D_OUT_DIR" \
    ${PLOT_MONOEXP_D_EXTRA_ARGS:-}
