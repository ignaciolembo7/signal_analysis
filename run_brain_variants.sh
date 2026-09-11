#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
LOG_DIR="${LOG_DIR:-$PROJECT_ROOT/logs/signal_analysis}"

DWI_LEVEL="${DWI_LEVEL:-den_gr-topup}"
ROI_VARIANTS="${ROI_VARIANTS:-plain erode1 sket1}"

mkdir -p "$LOG_DIR"

for roi_variant in $ROI_VARIANTS; do
    log_file="$LOG_DIR/${DWI_LEVEL}--${roi_variant}--brains.log"
    echo "Running brain analysis for ${DWI_LEVEL}--${roi_variant}"
    echo "Log: $log_file"

    (
        export DWI_LEVEL
        export ROI_VARIANT="$roi_variant"

        # Prevent unrelated overrides inherited from the calling shell from
        # routing different variants into the same input or output folder.
        unset SIGNAL_EXTRACTION_TAG ANALYSIS_TAG ANALYSIS_ROOT
        unset RESULTS_ROOT RESULTS_ROOTS MASTER_PARQUET MASTER_XLSX PLOT_OUT_ROOT

        bash "$SCRIPT_DIR/run_dataset.sh" brain ogse \
            ingest rotate grad_correction

        # Use physical b-values instead of b-step indices so this selection
        # works for both the legacy 10-point and newer 4-point acquisitions.
        ALPHA_N=1 \
        DPROJ_DIRS="long tra" \
        ALPHA_EXTRA_ARGS="--bvalmax 2000 --roi-bvalmax Syringe=980 --roi-bvalmax Right-Lateral-Ventricle=500 --roi-bvalmax Left-Lateral-Ventricle=500" \
            bash "$SCRIPT_DIR/run_dataset.sh" brain ogse alpha

        PLOT_ROW_KIND=signal \
            bash "$SCRIPT_DIR/run_dataset.sh" brain ogse plot_signal

        PLOT_ROW_KIND=signal_rotated \
            bash "$SCRIPT_DIR/run_dataset.sh" brain ogse plot_signal

        bash "$SCRIPT_DIR/run_dataset.sh" brain ogse export_master_xlsx
    ) >"$log_file" 2>&1

    echo "Completed ${DWI_LEVEL}--${roi_variant}"
done

echo "All requested brain variants completed."


# CLI command

# nohup bash repos/signal_analysis/run_brain_variants.sh >/tmp/run_brain_variants.nohup.log 2>&1 &
