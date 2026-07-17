#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help ingest "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}"

PROCESS_SCRIPT="${PROCESS_SCRIPT:-$REPO_ROOT/scripts/data/process_one_results.py}"
PROCESS_OUT_ROOT="${PROCESS_OUT_ROOT:-$ANALYSIS_ROOT/data/tables}"
RESULTS_GLOB="${RESULTS_GLOB:-*_results.xlsx}"

pipeline_require_file "$PROCESS_SCRIPT" "process script"
pipeline_require_file "$PARAMS_XLSX" "sequence params"
mkdir -p "$PROCESS_OUT_ROOT" "$(dirname "$MASTER_PARQUET")"

echo "Results root : $RESULTS_ROOT"
if [[ -n "${RESULTS_ROOTS:-}" ]]; then
    echo "Results roots: $RESULTS_ROOTS"
fi
echo "Master       : $MASTER_PARQUET"

mapfile -t sheet_names < <("$PY" -c '
import sys
import pandas as pd
for name in pd.ExcelFile(sys.argv[1]).sheet_names:
    print(name)
' "$PARAMS_XLSX")

if [[ ${#sheet_names[@]} -eq 0 ]]; then
    echo "ERROR: No sheets found in PARAMS_XLSX: $PARAMS_XLSX" >&2
    exit 1
fi
echo "Sheets       : ${sheet_names[*]}"

is_known_sheet() {
    local candidate="$1" s
    for s in "${sheet_names[@]}"; do
        [[ "$candidate" == "$s" ]] && return 0
    done
    return 1
}

count=0
roots="${RESULTS_ROOTS:-$RESULTS_ROOT}"
for root in $roots; do
    if [[ ! -d "$root" ]]; then
        echo "ERROR: Results root not found: $root" >&2
        exit 1
    fi

    # The pipeline contract is a BIDS-session Results folder containing
    # *_results.xlsx files directly.
    target_dirs=()
    if find "$root" -maxdepth 1 -type f -name "$RESULTS_GLOB" | read -r _; then
        target_dirs+=("$root")
    elif [[ -d "$root/Results" ]] && find "$root/Results" -maxdepth 1 -type f -name "$RESULTS_GLOB" | read -r _; then
        target_dirs+=("$root/Results")
    else
        while IFS= read -r d; do
            [[ -z "$d" ]] && continue
            target_dirs+=("$d")
        done < <(find "$root" -path "*/Results" -type d | while read -r results_dir; do
            find "$results_dir" -maxdepth 1 -type f -name "$RESULTS_GLOB" | read -r _ && printf '%s\n' "$results_dir"
        done | sort)
    fi

    if [[ ${#target_dirs[@]} -eq 0 ]]; then
        echo "WARNING: no pipeline Results folders with $RESULTS_GLOB found under $root" >&2
        continue
    fi

    for d in "${target_dirs[@]}"; do
        echo "Dataset dir  : $d"
        while read -r file; do
            [[ -z "$file" ]] && continue
            count=$((count + 1))
            "$PY" "$PROCESS_SCRIPT" "$file" "$PARAMS_XLSX" \
                --out_dir "$PROCESS_OUT_ROOT" \
                --master-parquet "$MASTER_PARQUET"
        done < <(find "$d" -maxdepth 1 -type f -name "$RESULTS_GLOB" | sort)
    done
done

echo "Done ingest. Files processed: $count"
