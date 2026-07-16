from __future__ import annotations

import repo_bootstrap  # noqa: F401

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from data_processing.io import infer_layout_from_filename, read_result_xls, write_table_outputs
from data_processing.match_params import parse_results_filename, select_params_row
from data_processing.master_table import append_master_rows
from data_processing.params import read_sequence_params_xlsx
from data_processing.reshape import to_long
from data_processing.result_signals import (
    add_clean_sequence_params,
    add_direct_g_derivatives,
    add_g_and_b_derivatives,
    add_s0_and_value_norm,
    aggregate_g_results,
    detect_gradient_input_kind,
    direction_from_vector_path,
    extract_clean_sequence_params,
    find_matching_gradient_vector,
    first_present,
    is_single_point_results,
    master_analysis_id,
    merge_into_group_curve,
    output_stem_from_results_file,
    output_strip_tokens,
    sha256_file,
    truthy_env,
)
from data_processing.schema import finalize_clean_signal_long


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    ap.add_argument("results_file", type=Path)
    ap.add_argument("params_xlsx", type=Path)
    ap.add_argument("--out_dir", type=Path, default=Path("analysis/ogse_experiments/data/tables"))
    ap.add_argument("--gamma", type=float, default=267.5221900, help="1/(ms*mT)")
    ap.add_argument(
        "--oneg",
        action="store_true",
        default=truthy_env("NOGSE_ONEG"),
        help="Treat one-g-per-sequence result files as points of the same direct-g curve.",
    )
    ap.add_argument(
        "--strip-output-token",
        "--output-stem-strip",
        dest="strip_output_tokens",
        action="append",
        default=[],
        help="Token to remove from generated output stems. Can be repeated or comma/space separated.",
    )
    ap.add_argument(
        "--master-parquet",
        type=Path,
        default=None,
        help="Optional master table parquet to append the processed signal rows to.",
    )
    return ap


def result_tables_to_long(stats: dict[str, pd.DataFrame], results_file: Path) -> tuple[pd.DataFrame, str]:
    gradient_input_kind = detect_gradient_input_kind(stats)
    meta = parse_results_filename(results_file)

    if gradient_input_kind == "b":
        layout = infer_layout_from_filename(results_file)
        ndirs = meta.ndirs or layout.ndirs
        nbvals = meta.nbvals or layout.nbvals
        if ndirs is None or nbvals is None:
            raise SystemExit(f"Could not infer ndirs/nbvals from filename: {results_file.name}")
        return to_long(stats, ndirs=ndirs, nbvals=nbvals, source_file=results_file.name), gradient_input_kind

    if is_single_point_results(stats):
        vector_path = find_matching_gradient_vector(results_file)
        direction_label = direction_from_vector_path(vector_path)
        return (
            aggregate_g_results(stats, direction_label=direction_label, source_file=results_file.name),
            gradient_input_kind,
        )

    layout = infer_layout_from_filename(results_file)
    ndirs = meta.ndirs or layout.ndirs
    nbvals = meta.nbvals or layout.nbvals
    if ndirs is None or nbvals is None:
        raise SystemExit(f"Could not infer ndirs/nbvals from filename: {results_file.name}")

    any_df = next(iter(stats.values()))
    gcol = first_present([str(c) for c in any_df.columns], ["gvalues", "gval", "g"])
    if gcol is None:
        raise SystemExit(f"Could not find a g column in: {list(any_df.columns)}")

    df_long = to_long(
        stats,
        ndirs=ndirs,
        nbvals=nbvals,
        bcol=gcol,
        out_col="g",
        source_file=results_file.name,
    )
    df_long["bvalue"] = np.nan
    return df_long, gradient_input_kind


def attach_features(
    df_long: pd.DataFrame,
    *,
    gradient_input_kind: str,
    clean_params,
    gamma: float,
) -> pd.DataFrame:
    if gradient_input_kind == "b":
        out = add_g_and_b_derivatives(
            df_long,
            gamma=gamma,
            N=int(clean_params.N),
            delta_ms=float(clean_params.delta_ms),
            Delta_app_ms=float(clean_params.Delta_app_ms),
            g_thorsten=clean_params.g_thorsten,
        )
        return add_s0_and_value_norm(out)

    out = df_long.copy()
    if "g_thorsten" not in out.columns:
        out["g_thorsten"] = pd.NA
    out["bvalue_g"] = pd.NA
    out["bvalue_g_lin_max"] = pd.NA
    out["bvalue_thorsten"] = pd.NA
    out = add_direct_g_derivatives(out, gamma=gamma)
    out["S0"] = pd.NA
    out["value_norm"] = pd.NA
    zero_mask = np.isclose(pd.to_numeric(out["g"], errors="coerce"), 0.0, rtol=0.0, atol=1e-9)
    out.loc[zero_mask, "S0"] = pd.to_numeric(out.loc[zero_mask, "value"], errors="coerce")
    out.loc[zero_mask, "value_norm"] = 1.0
    return out


def save_processed_rows(
    df_long: pd.DataFrame,
    *,
    args: argparse.Namespace,
    clean_params,
    stats: dict[str, pd.DataFrame],
    gradient_input_kind: str,
) -> tuple[pd.DataFrame, Path]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    exp_dir = args.out_dir / clean_params.sheet
    exp_dir.mkdir(parents=True, exist_ok=True)

    grouped_g_curve = bool(args.oneg) and gradient_input_kind == "g" and is_single_point_results(stats)
    if args.oneg and not grouped_g_curve:
        raise ValueError("--oneg expects a direct-g results file with one row per stat sheet.")

    out_stem = output_stem_from_results_file(
        args.results_file,
        grouped_g_curve=grouped_g_curve,
        strip_tokens=output_strip_tokens(args.strip_output_tokens),
    )
    out_path = exp_dir / f"{out_stem}.long.parquet"

    if grouped_g_curve:
        existing = pd.read_parquet(out_path) if out_path.exists() else None
        df_to_save = merge_into_group_curve(existing, df_long, gamma=float(args.gamma))
    else:
        df_to_save = df_long

    write_table_outputs(df_to_save, out_path, xlsx_path=out_path.with_suffix(".xlsx"))
    return df_to_save, out_path


def append_to_master(df_to_save: pd.DataFrame, *, results_file: Path, master_path: Path) -> None:
    master_rows = df_to_save.copy()
    master_rows["source_path"] = str(results_file.resolve())
    master_rows["source_hash"] = sha256_file(results_file)
    analysis_id = master_analysis_id(master_rows, row_kind="signal")
    append_master_rows(
        master_path if master_path.exists() else None,
        master_rows,
        row_kind="signal",
        analysis_id=analysis_id,
        out_path=master_path,
    )


def print_selected_params(clean_params) -> None:
    print("Selected params (clean):")
    print(
        pd.Series(
            {
                "sheet": clean_params.sheet,
                "subj": clean_params.subj,
                "protocol": clean_params.protocol,
                "sequence": clean_params.sequence,
                "group": clean_params.group,
                "G": clean_params.G,
                "TN": clean_params.TN,
                "x": clean_params.x,
                "y": clean_params.y,
                "type": clean_params.type,
                "Hz": clean_params.Hz,
                "bmax": clean_params.bmax,
                "N": clean_params.N,
                "delta_ms": clean_params.delta_ms,
                "Delta_app_ms": clean_params.Delta_app_ms,
                "max_dur_ms": clean_params.max_dur_ms,
                "tm_ms": clean_params.tm_ms,
                "td_ms": clean_params.td_ms,
                "TE": clean_params.TE,
                "TR": clean_params.TR,
            }
        ).to_string()
    )


def main() -> None:
    args = build_parser().parse_args()
    stats = read_result_xls(args.results_file)
    meta = parse_results_filename(args.results_file)

    df_long, gradient_input_kind = result_tables_to_long(stats, args.results_file)

    params = read_sequence_params_xlsx(args.params_xlsx)
    row = select_params_row(params, meta)
    if row is None:
        logger.warning(
            f"Skipping {args.results_file.name}: No matching parameters found in Excel. "
            f"sheet={meta.sheet!r}, seq={meta.seq}, Hz={meta.Hz}, d_ms={meta.d_ms}, "
            f"delta_ms={meta.delta_ms}, Delta_ms={meta.Delta_ms}"
        )
        return

    clean_params = extract_clean_sequence_params(row, meta)
    df_long = add_clean_sequence_params(df_long, clean_params)
    df_long["gradient_axis_kind"] = gradient_input_kind
    df_long = attach_features(
        df_long,
        gradient_input_kind=gradient_input_kind,
        clean_params=clean_params,
        gamma=float(args.gamma),
    )
    df_long["one_g_per_sequence"] = bool(args.oneg)
    df_long = finalize_clean_signal_long(df_long)
    df_long = df_long.sort_values(["stat", "roi", "direction", "b_step"], kind="stable").reset_index(drop=True)

    df_to_save, out_path = save_processed_rows(
        df_long,
        args=args,
        clean_params=clean_params,
        stats=stats,
        gradient_input_kind=gradient_input_kind,
    )

    if args.master_parquet is not None:
        append_to_master(df_to_save, results_file=args.results_file, master_path=args.master_parquet)

    print_selected_params(clean_params)
    print("Saved:", out_path)
    if args.master_parquet is not None:
        print("Appended master:", args.master_parquet)


if __name__ == "__main__":
    main()
