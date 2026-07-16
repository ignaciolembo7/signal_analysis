#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/helpers/master_table_common.sh"
pipeline_maybe_step_help rotate "$@"
pipeline_setup_common
pipeline_set_dataset_defaults "${TYPE_SUBJ:-${DATASET:?TYPE_SUBJ or DATASET is required}}"

ROTATE_SCRIPT="${ROTATE_SCRIPT:-$REPO_ROOT/scripts/data/rotate_ogse_tensor.py}"
DIRS_TXT="${DIRS_TXT:-$REPO_ROOT/assets/dirs/dirs_6.txt}"

if [[ "$TYPE_SUBJ" == "phantom" ]]; then
    TRA_AXES="${TRA_AXES:-x,y}"
    LONG_AXES="${LONG_AXES:-z}"
else
    TRA_AXES="${TRA_AXES:-y,z}"
    LONG_AXES="${LONG_AXES:-x}"
fi

pipeline_require_file "$ROTATE_SCRIPT" "rotate script"
pipeline_require_file "$MASTER_PARQUET" "master table"
pipeline_require_file "$DIRS_TXT" "dirs txt"

count=0
while IFS=$'\t' read -r subj sheet td_ms n hz; do
    count=$((count + 1))
    echo "Rotating: subj=$subj sheet=$sheet td_ms=$td_ms N=$n Hz=$hz"
    "$PY" "$ROTATE_SCRIPT" \
        --master-parquet "$MASTER_PARQUET" \
        --subj "$subj" \
        --sheet "$sheet" \
        --td_ms "$td_ms" \
        --N "$n" \
        --Hz "$hz" \
        --dirs_txt "$DIRS_TXT" \
        --tra_axes "$TRA_AXES" \
        --long_axes "$LONG_AXES" \
        --no-legacy-output \
        ${ROTATE_EXTRA_ARGS:-}
done < <("$PY" - "$MASTER_PARQUET" "${MASTER_SUBJ:-ALL}" "${MASTER_SHEET:-ALL}" <<'PY'
import pandas as pd
import sys

path, subj_sel, sheet_sel = sys.argv[1:4]
df = pd.read_parquet(path)
sig = df[df["row_kind"].astype(str).eq("signal")].copy()
if subj_sel != "ALL":
    sig = sig[sig["subj"].astype(str).eq(subj_sel)]
if sheet_sel != "ALL":
    sig = sig[sig["sheet"].astype(str).eq(sheet_sel)]
cols = [c for c in ["subj", "sheet", "td_ms", "N", "Hz"] if c in sig.columns]
for row in sig[cols].drop_duplicates().sort_values(cols).itertuples(index=False):
    print("\t".join(str(x) for x in row))
PY
)

if [[ "$count" -eq 0 ]]; then
    echo "ERROR: rotate selected no signal rows." >&2
    echo "       MASTER_SUBJ=${MASTER_SUBJ:-ALL}" >&2
    echo "       MASTER_SHEET=${MASTER_SHEET:-ALL}" >&2
    echo "       Check the exact values in master.long.parquet columns 'subj' and 'sheet'." >&2
    exit 1
fi

echo "Done rotate. Groups processed: $count"
