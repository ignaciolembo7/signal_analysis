#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/helpers/master_table_common.sh"

pipeline_setup_common

type_subj="${TYPE_SUBJ:-${DATASET:-}}"
type_seq="${TYPE_SEQ:-}"
steps=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            if [[ ${#steps[@]} -gt 0 ]]; then
                steps+=("$1")
                shift
            else
                pipeline_usage
                exit 0
            fi
            ;;
        --type-subj|--type_subj|--dataset)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires brain or phantom" >&2; exit 2; }
            type_subj="$2"
            shift 2
            ;;
        --type-seq|--type_seq)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires ogse or nogse" >&2; exit 2; }
            type_seq="$2"
            shift 2
            ;;
        --results-root)
            [[ $# -ge 2 ]] || { echo "ERROR: --results-root requires a path" >&2; exit 2; }
            RESULTS_ROOTS="${RESULTS_ROOTS:+$RESULTS_ROOTS }$2"
            export RESULTS_ROOTS
            shift 2
            ;;
        brain|brains|phantom|phantoms)
            type_subj="$1"
            shift
            ;;
        ogse|nogse)
            type_seq="$1"
            shift
            ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do
                steps+=("$1")
                shift
            done
            ;;
        *)
            steps+=("$1")
            shift
            ;;
    esac
done

if [[ -z "${type_subj// }" ]]; then
    echo "ERROR: missing type_subj (use brain or phantom)" >&2
    pipeline_usage >&2
    exit 2
fi

pipeline_set_dataset_defaults "$type_subj" "${type_seq:-ogse}"
pipeline_run_steps "${steps[@]}"
