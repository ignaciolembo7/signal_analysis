from __future__ import annotations

import repo_bootstrap  # noqa: F401

import argparse
from pathlib import Path

from data_processing.master_table import build_analysis_id_from_columns, load_master_table, select_plot_signal, split_selector_values
from plotting.ogse.signal_vs_g import plot_ogse_signal_summary


def _master_selectors(args: argparse.Namespace) -> dict[str, object]:
    selectors: dict[str, object] = {}
    for arg_name, col_name in [
        ("analysis_id", "analysis_id"),
        ("subj", "subj"),
        ("sheet", "sheet"),
        ("roi", "roi"),
        ("direction", "direction"),
    ]:
        values = split_selector_values(getattr(args, arg_name, None))
        if values is not None:
            selectors[col_name] = values
    for col in ("td_ms", "N", "Hz"):
        value = getattr(args, col, None)
        if value is not None:
            selectors[col] = float(value)
    return selectors


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot OGSE signal curves from master-table signal rows.")
    ap.add_argument("--master-parquet", type=Path, required=True, help="Read signal rows from the master table.")
    ap.add_argument("--row-kind", choices=["signal", "signal_rotated"], default="signal_rotated")
    ap.add_argument("--analysis-id", action="append", default=None)
    ap.add_argument("--subj", action="append", default=None)
    ap.add_argument("--sheet", action="append", default=None)
    ap.add_argument("--roi", action="append", default=None)
    ap.add_argument("--direction", action="append", default=None)
    ap.add_argument("--td_ms", type=float, default=None)
    ap.add_argument("--N", type=float, default=None)
    ap.add_argument("--Hz", type=float, default=None)
    ap.add_argument("--out_root", "--out_dir", dest="out_root", type=Path, default=Path("plots/ogse_vs_g"))
    ap.add_argument("--ycol", "--y_col", dest="ycol", default="value_norm")
    ap.add_argument("--xcol", default="g_thorsten")
    ap.add_argument("--stat", default="avg")
    ap.add_argument("--no_ylim", action="store_true")
    args = ap.parse_args()

    master = load_master_table(args.master_parquet)
    df = select_plot_signal(
        master,
        rotated=args.row_kind == "signal_rotated",
        **_master_selectors(args),
    )
    if df.empty:
        raise ValueError("No master signal rows matched the requested selectors.")
    try:
        exp_id = build_analysis_id_from_columns(
            df,
            columns=("subj", "sheet", "td_ms", "N", "Hz"),
            prefix=args.row_kind,
        )
    except ValueError:
        exp_id = f"{args.row_kind}_master_selection"
    out_dir = args.out_root / exp_id
    ylim = None if args.no_ylim else (0.0, 1.0)

    outputs = plot_ogse_signal_summary(
        df,
        out_dir,
        xcol=str(args.xcol),
        ycol=str(args.ycol),
        stat=str(args.stat),
        ylim=ylim,
    )

    print("Saved plots to:", out_dir)
    print("Generated:", len(outputs))


if __name__ == "__main__":
    main()
