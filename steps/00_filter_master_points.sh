#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help filter_master_points "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}"

FILTER_MASTER_SCRIPT="${FILTER_MASTER_SCRIPT:-$REPO_ROOT/scripts/data/filter_master_table.py}"
FILTERED_MASTER_PARQUET="${FILTERED_MASTER_PARQUET:-${MASTER_LAST_POINTS_PARQUET:-$ANALYSIS_ROOT/master.last_points.long.parquet}}"

pipeline_require_file "$FILTER_MASTER_SCRIPT" "master filter script"
pipeline_require_file "$MASTER_PARQUET" "master table"
if [[ -z "${MASTER_LAST_POINTS_BY_TD:-}" || "${MASTER_LAST_POINTS_BY_TD^^}" == "ALL" ]]; then
    echo "ERROR: set MASTER_LAST_POINTS_BY_TD, for example MASTER_LAST_POINTS_BY_TD=\"120:8=6,120:4=4,210=8\"." >&2
    exit 2
fi

mkdir -p "$(dirname "$FILTERED_MASTER_PARQUET")"

"$PY" "$FILTER_MASTER_SCRIPT" "$MASTER_PARQUET" \
    --out-parquet "$FILTERED_MASTER_PARQUET" \
    --last-points-by-td-n "$MASTER_LAST_POINTS_BY_TD"

echo "Filtered master table: $FILTERED_MASTER_PARQUET"
