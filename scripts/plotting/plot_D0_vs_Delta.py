from __future__ import annotations

import repo_bootstrap  # noqa: F401

import argparse
from pathlib import Path

import pandas as pd

from data_processing.master_table import filter_master_rows, load_master_table, split_selector_values
from data_processing.io import write_xlsx_csv_outputs
from monoexp_fitting.plot_D0_vs_Delta import (
    load_all_measurements,
    load_selected_bstep_map,
    plot_all_groups,
)
from tc_fittings.alpha_macro_summary import load_dproj_measurements_from_table


def _pick_column(df: pd.DataFrame, candidates: list[str]) -> str:
    lower_to_name = {str(col).strip().lower(): str(col) for col in df.columns}
    for candidate in candidates:
        col = lower_to_name.get(candidate.lower())
        if col is not None:
            return col
    raise KeyError(f"Missing one of columns {candidates}. Available columns: {list(df.columns)}")


def _read_summary_alpha_groups(path: Path) -> set[tuple[str, str, str]]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        summary = pd.read_csv(path)
    elif suffix in {".xlsx", ".xls"}:
        summary = pd.read_excel(path)
    elif suffix == ".parquet":
        summary = pd.read_parquet(path)
    else:
        raise ValueError(f"Unsupported summary format: {path}")

    subj_col = _pick_column(summary, ["subj", "brain"])
    roi_col = _pick_column(summary, ["roi", "region"])
    direction_col = _pick_column(summary, ["direction", "direccion"])
    groups: set[tuple[str, str, str]] = set()
    for _, row in summary[[subj_col, roi_col, direction_col]].dropna().iterrows():
        groups.add(
            (
                str(row[subj_col]).strip(),
                str(row[roi_col]).strip().replace("_norm", "").lower(),
                str(row[direction_col]).strip(),
            )
        )
    return groups


def _filter_to_summary_groups(df: pd.DataFrame, groups: set[tuple[str, str, str]]) -> pd.DataFrame:
    keys = pd.DataFrame(
        {
            "subj": df["subj"].astype(str).str.strip(),
            "roi": df["roi"].astype(str).str.strip().str.replace("_norm", "", regex=False).str.lower(),
            "direction": df["direction"].astype(str).str.strip(),
        },
        index=df.index,
    )
    row_keys = list(zip(keys["subj"], keys["roi"], keys["direction"]))
    return df.loc[[key in groups for key in row_keys]].copy()


def _clean_old_d_vs_delta_pngs(out_dir: Path) -> None:
    for path in out_dir.glob("D_vs_delta_app_*.png"):
        path.unlink()


def _master_selectors(args: argparse.Namespace, *, N: float | None) -> dict[str, object]:
    selectors: dict[str, object] = {"row_kind": "signal_rotated"}
    for arg_name, col_name in [
        ("analysis_id", "analysis_id"),
        ("subjs", "subj"),
        ("sheets", "sheet"),
        ("rois", "roi"),
        ("dirs", "direction"),
    ]:
        values = split_selector_values(getattr(args, arg_name, None))
        if values is not None:
            selectors[col_name] = values
    if N is not None:
        selectors["N"] = float(N)
    if args.Hz is not None:
        selectors["Hz"] = float(args.Hz)
    return selectors


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build D vs Delta_app_ms curves from *.Dproj.long.parquet tables."
    )
    ap.add_argument("--dproj-root", default=None, help="Root folder with *.Dproj.long.parquet tables.")
    ap.add_argument("--master-parquet", type=Path, default=None, help="Read D_proj values from signal_rotated rows in the master table.")
    ap.add_argument("--analysis-id", action="append", default=None, help="Master analysis_id selector.")
    ap.add_argument("--sheets", nargs="+", default=None, help="Master sheet selector.")
    ap.add_argument("--pattern", default="**/*.Dproj.long.parquet", help="Relative glob inside dproj-root.")
    ap.add_argument("--out-dir", required=True, help="Output folder for plots and the combined table.")
    ap.add_argument("--subjs", nargs="+", default=None, help="Subjects/phantoms to include, for example: BRAIN-3 LUDG-2 PHANTOM3.")
    ap.add_argument("--rois", nargs="+", default=None, help="ROIs to include.")
    ap.add_argument("--dirs", nargs="+", default=None, help="Directions to include. Default: all directions, or summary_alpha groups when --summary-alpha is available.")

    selector = ap.add_mutually_exclusive_group()
    selector.add_argument("--N", type=float, default=None, help="Filter by N.")
    selector.add_argument("--Hz", type=float, default=None, help="Filter by Hz.")

    ap.add_argument("--bvalue-decimals", type=int, default=1, help="Decimals used to round bvalue before grouping.")
    ap.add_argument(
        "--bvalmax",
        type=int,
        default=None,
        help=(
            "Bstep (1-based) to use for the horizontal line and plot alpha. "
            "If omitted, the highest bvalue is used."
        ),
    )
    ap.add_argument("--reference-D0", type=float, default=0.0032, help="Reference value used to annotate alpha in the plot.")
    ap.add_argument("--reference-D0-error", type=float, default=0.0000283512, help="Reference value error.")
    ap.add_argument(
        "--summary-alpha",
        default=None,
        help=(
            "Optional path to summary_alpha_values (xlsx/csv/parquet). "
            "If omitted and <out-dir>/summary_alpha_values.xlsx exists, it is used automatically."
        ),
    )
    args = ap.parse_args()

    N = 1.0 if args.N is None and args.Hz is None else args.N

    if args.master_parquet is not None:
        master = load_master_table(args.master_parquet)
        dproj = filter_master_rows(master, **_master_selectors(args, N=N))
        if "D_proj" in dproj.columns:
            dproj = dproj[dproj["D_proj"].notna()].copy()
        df = load_dproj_measurements_from_table(
            dproj,
            directions=args.dirs,
            rois=args.rois,
            subjs=args.subjs,
            N=N,
            Hz=args.Hz,
            bvalue_decimals=int(args.bvalue_decimals),
        )
    else:
        if not args.dproj_root:
            raise ValueError("Pass --master-parquet or --dproj-root.")
        df = load_all_measurements(
            args.dproj_root,
            pattern=args.pattern,
            dirs=args.dirs,
            rois=args.rois,
            subjs=args.subjs,
            N=N,
            Hz=args.Hz,
            bvalue_decimals=int(args.bvalue_decimals),
        )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    default_summary = out_dir / "summary_alpha_values.xlsx"
    summary_alpha_path = Path(args.summary_alpha) if args.summary_alpha is not None else (default_summary if default_summary.exists() else None)

    selected_bstep_by_group = None
    if summary_alpha_path is not None:
        selected_bstep_by_group = load_selected_bstep_map(summary_alpha_path)
        if args.dirs is None:
            summary_groups = _read_summary_alpha_groups(summary_alpha_path)
            before = len(df)
            df = _filter_to_summary_groups(df, summary_groups)
            if df.empty:
                raise ValueError(
                    f"No D-vs-Delta groups matched summary_alpha={summary_alpha_path}. "
                    "Check that subj/roi/direction labels agree."
                )
            print(f"[INFO] Filtered plot data to summary_alpha groups: {len(df)}/{before} rows.")
        print(
            f"[INFO] Loaded selected_bstep map from {summary_alpha_path} "
            f"({len(selected_bstep_by_group)} group entries)."
        )
        if args.bvalmax is not None:
            print(
                "[INFO] --bvalmax is used only as fallback for groups missing in summary_alpha."
            )

    _clean_old_d_vs_delta_pngs(out_dir)
    plot_all_groups(
        df,
        out_dir=out_dir,
        selected_bstep=args.bvalmax,
        selected_bstep_by_group=selected_bstep_by_group,
        reference_D0=float(args.reference_D0),
        reference_D0_error=float(args.reference_D0_error),
    )

    write_xlsx_csv_outputs(
        df,
        out_dir / "D_vs_delta_app.combined.xlsx",
        csv_path=out_dir / "D_vs_delta_app.combined.csv",
    )
    print(f"[OK] Plots + tables in: {out_dir}")


if __name__ == "__main__":
    main()
