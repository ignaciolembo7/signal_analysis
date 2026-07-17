#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help export_master_xlsx "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}"

EXPORT_MASTER_SCRIPT="${EXPORT_MASTER_SCRIPT:-$REPO_ROOT/scripts/data/export_master_table.py}"
MASTER_XLSX="${MASTER_XLSX:-${MASTER_PARQUET%.parquet}.xlsx}"

pipeline_require_file "$EXPORT_MASTER_SCRIPT" "master export script"
pipeline_require_file "$MASTER_PARQUET" "master table"
mkdir -p "$(dirname "$MASTER_XLSX")"

args=("$MASTER_PARQUET" --out-xlsx "$MASTER_XLSX")
if [[ -n "${EXPORT_MASTER_ROW_KIND:-}" ]]; then
    for row_kind in ${EXPORT_MASTER_ROW_KIND//,/ }; do
        args+=(--row-kind "$row_kind")
    done
fi
if [[ -n "${EXPORT_MASTER_HEAD:-}" ]]; then
    args+=(--head "$EXPORT_MASTER_HEAD")
fi

"$PY" "$EXPORT_MASTER_SCRIPT" "${args[@]}"
