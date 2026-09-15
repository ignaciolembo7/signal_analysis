from __future__ import annotations

import repo_bootstrap  # noqa: F401

import argparse
from pathlib import Path

import pandas as pd

from data_processing.io import write_table_outputs
from data_processing.master_table import (
    filter_master_rows,
    load_fit_params_table,
    load_master_table,
    split_selector_values,
    write_master_table,
)
from tc_fittings.alpha_macro_summary import (
    compute_alpha_macro_summary,
    load_dproj_measurements,
    load_dproj_measurements_from_table,
    parse_direction_aliases,
    plot_alpha_macro_vs_roi,
    write_alpha_macro_outputs,
)


def _safe_plot_alpha_macro_vs_roi(*args, **kwargs) -> None:
    try:
        plot_alpha_macro_vs_roi(*args, **kwargs)
    except ValueError as exc:
        msg = str(exc)
        if "No data available for plotting" in msg:
            print(f"[INFO] Skipping plot: {msg}")
            return
        raise


def _ordered_unique(series: pd.Series) -> list[str]:
    return list(dict.fromkeys(series.dropna().astype(str).tolist()))


def _canonical_roi_name(value: str) -> str:
    return str(value).strip().replace("_norm", "").lower()


def _parse_roi_bvalmax(
    items: list[str] | None,
) -> tuple[dict[str, float], dict[tuple[str, str], float]]:
    if not items:
        return {}, {}
    roi_values: dict[str, float] = {}
    roi_direction_values: dict[tuple[str, str], float] = {}
    for raw in items:
        token = str(raw).strip()
        if "=" not in token:
            raise ValueError(
                f"Invalid --roi-bvalmax value {raw!r}. Use ROI=X or ROI:DIRECTION=X, "
                "for example AntCC=2000 or fiber1:long=500."
            )
        selector, bstep = token.split("=", 1)
        selector = selector.strip()
        bstep = bstep.strip()
        if ":" in selector:
            roi, direction = (part.strip() for part in selector.rsplit(":", 1))
        else:
            roi, direction = selector, None
        if not roi:
            raise ValueError(f"Invalid ROI in --roi-bvalmax {raw!r}.")
        if direction is not None and not direction:
            raise ValueError(f"Invalid direction in --roi-bvalmax {raw!r}.")
        try:
            value = float(bstep)
        except ValueError as exc:
            raise ValueError(
                f"Invalid BSTEP/BVALUE in --roi-bvalmax {raw!r}. It must be a number >= 1."
            ) from exc
        if value < 1:
            raise ValueError(
                f"Invalid BSTEP/BVALUE in --roi-bvalmax {raw!r}. It must be a number >= 1."
            )
        if direction is None:
            roi_values[roi] = value
        else:
            roi_direction_values[(roi, direction)] = value
    return roi_values, roi_direction_values


def _load_measurements_from_args(args: argparse.Namespace) -> pd.DataFrame:
    if args.master_parquet is not None:
        master = load_master_table(args.master_parquet)
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
        if args.N is not None:
            selectors["N"] = float(args.N)
        if args.Hz is not None:
            selectors["Hz"] = float(args.Hz)
        df = filter_master_rows(master, **selectors)
        if "D_proj" in df.columns:
            df = df[pd.to_numeric(df["D_proj"], errors="coerce").notna()].copy()
        return load_dproj_measurements_from_table(
            df,
            subjs=args.subjs,
            rois=args.rois,
            directions=args.dirs,
            N=args.N,
            Hz=args.Hz,
            bvalue_decimals=int(args.bvalue_decimals),
        )

    if args.combined_table is not None:
        path = Path(args.combined_table)
        if not path.exists():
            raise FileNotFoundError(path)
        if path.suffix.lower() == ".csv":
            df = pd.read_csv(path)
        elif path.suffix.lower() in {".xlsx", ".xls"}:
            df = pd.read_excel(path)
        elif path.suffix.lower() == ".parquet":
            df = pd.read_parquet(path)
        else:
            raise ValueError(f"Unsupported format for combined_table={path}")

        if args.subjs is not None and "subj" in df.columns:
            df = df[df["subj"].astype(str).isin([str(x) for x in args.subjs])]
        if args.rois is not None and "roi" in df.columns:
            df = df[df["roi"].astype(str).isin([str(x) for x in args.rois])]
        if args.dirs is not None and "direction" in df.columns:
            df = df[df["direction"].astype(str).isin([str(x) for x in args.dirs])]
        if args.N is not None and "N" in df.columns:
            df = df[df["N"].astype(float).sub(float(args.N)).abs() <= 1e-6]
        if args.Hz is not None and "Hz" in df.columns:
            df = df[df["Hz"].astype(float).sub(float(args.Hz)).abs() <= 1e-6]
        return df

    if not args.dproj_root:
        raise ValueError("Pass --combined-table or --dproj-root.")

    return load_dproj_measurements(
        args.dproj_root,
        pattern=args.pattern,
        subjs=args.subjs,
        rois=args.rois,
        directions=args.dirs,
        N=args.N,
        Hz=args.Hz,
        bvalue_decimals=int(args.bvalue_decimals),
    )


def _default_master_fit_params(args: argparse.Namespace) -> Path | None:
    if args.no_master_fit_params:
        return None
    if args.master_fit_params is not None:
        return Path(args.master_fit_params)
    if args.master_parquet is not None:
        return Path(args.master_parquet).with_name("master_fit_params.parquet")
    return None


def _append_alpha_macro_fit_params(df_summary: pd.DataFrame, args: argparse.Namespace) -> Path | None:
    out_path = _default_master_fit_params(args)
    if out_path is None:
        return None

    rows = df_summary.copy()
    rows["fit_kind"] = "alpha_macro_summary"
    rows["fit_model"] = "Dproj_reference_D0"
    rows["fit_source"] = str(args.master_parquet or args.combined_table or args.dproj_root or "")
    rows["reference_D0"] = float(args.reference_D0)
    rows["reference_D0_error"] = float(args.reference_D0_error)
    if args.N is not None and "N" not in rows.columns:
        rows["N"] = float(args.N)
    if args.Hz is not None and "Hz" not in rows.columns:
        rows["Hz"] = float(args.Hz)

    if out_path.exists():
        existing = load_fit_params_table(out_path)
        rows = pd.concat([existing, rows], ignore_index=True, sort=False).drop_duplicates().reset_index(drop=True)
    else:
        rows = rows.drop_duplicates().reset_index(drop=True)
    write_table_outputs(rows, out_path, xlsx_path=out_path.with_suffix(".xlsx"))
    return out_path


def _annotate_master_alpha_macro(df_summary: pd.DataFrame, args: argparse.Namespace) -> tuple[Path, int] | None:
    if args.master_parquet is None or args.no_annotate_master_alpha:
        return None

    master = load_master_table(args.master_parquet)
    if "alpha_macro" not in master.columns:
        master["alpha_macro"] = pd.NA

    alpha = df_summary[["subj", "roi", "direction", "alpha_macro"]].copy()
    alpha["alpha_macro"] = pd.to_numeric(alpha["alpha_macro"], errors="coerce")
    alpha = alpha.dropna(subset=["alpha_macro"]).copy()
    alpha["_subj_key"] = alpha["subj"].astype(str).str.strip()
    alpha["_roi_key"] = alpha["roi"].astype(str).str.strip().str.replace("_norm", "", regex=False).str.lower()
    alpha["_direction_key"] = alpha["direction"].astype(str).str.strip()
    alpha = alpha.groupby(["_subj_key", "_roi_key", "_direction_key"], as_index=False).agg(alpha_macro=("alpha_macro", "mean"))

    eligible = pd.Series(True, index=master.index)
    if "value" in master.columns:
        eligible &= pd.to_numeric(master["value"], errors="coerce").notna()

    keys = pd.DataFrame(
        {
            "_row_index": master.index,
            "_subj_key": master["subj"].astype(str).str.strip(),
            "_roi_key": master["roi"].astype(str).str.strip().str.replace("_norm", "", regex=False).str.lower(),
            "_direction_key": master["direction"].astype(str).str.strip(),
        }
    ).loc[eligible]
    matched = keys.merge(alpha, on=["_subj_key", "_roi_key", "_direction_key"], how="left").dropna(subset=["alpha_macro"])
    if not matched.empty:
        values = matched.set_index("_row_index")["alpha_macro"]
        master.loc[values.index, "alpha_macro"] = values.astype(float)

    write_master_table(master, args.master_parquet)
    return args.master_parquet, int(len(matched))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Compute alpha_macro = <D0>/0.0032 for a given N/Hz and plot alpha_macro vs ROI."
    )
    ap.add_argument("--combined-table", type=Path, default=None, help="Combined table generated by plot_D0_vs_Delta.py.")
    ap.add_argument("--master-parquet", type=Path, default=None, help="Read D_proj values from signal_rotated rows in the master table.")
    ap.add_argument("--master-fit-params", type=Path, default=None, help="Append alpha_macro rows to this cumulative fit params table.")
    ap.add_argument("--no-master-fit-params", action="store_true", help="Do not append alpha_macro rows to master_fit_params.parquet.")
    ap.add_argument("--no-annotate-master-alpha", action="store_true", help="Do not write alpha_macro back to --master-parquet.")
    ap.add_argument("--dproj-root", default=None, help="Root with *.Dproj.long.parquet files when --combined-table is not passed.")
    ap.add_argument("--pattern", default="**/*.Dproj.long.parquet", help="Relative glob inside dproj-root.")
    ap.add_argument("--analysis-id", action="append", default=None, help="Master analysis_id selector. Can be repeated or comma-separated.")
    ap.add_argument("--subjs", nargs="+", default=None, help="Subjects/phantoms to include (for example BRAIN-3 LUDG-2 PHANTOM3).")
    ap.add_argument("--sheets", nargs="+", default=None, help="Master sheet selector.")
    ap.add_argument("--rois", nargs="+", default=None, help="ROIs to include.")
    ap.add_argument("--dirs", nargs="+", default=None, help="Raw directions to include. When omitted, no direction filter is applied.")

    selector = ap.add_mutually_exclusive_group()
    selector.add_argument("--N", type=float, default=None, help="Filter by N.")
    selector.add_argument("--Hz", type=float, default=None, help="Filter by Hz.")

    ap.add_argument("--bvalue-decimals", type=int, default=1, help="Decimals used to round bvalue before grouping.")
    ap.add_argument(
        "--bvalmax",
        type=float,
        default=None,
        help=(
            "Bstep or bvalue to use for alpha_macro. "
            "If the value is within the available candidate bstep range it is treated as a 1-based bstep; "
            "otherwise it is matched as a rounded bvalue. Example: --bvalmax 7 or --bvalmax 2000."
        ),
    )
    ap.add_argument(
        "--plot-bsteps",
        "--plot_bsteps",
        nargs="+",
        type=int,
        default=None,
        help=(
            "Candidate 1-based bvalue positions used both for alpha selection and, "
            "when passed to plot_D0_vs_Delta.py, for choosing curves to draw."
        ),
    )
    ap.add_argument(
        "--plot-bvalues",
        "--plot_bvalues",
        nargs="+",
        type=float,
        default=None,
        help=(
            "Candidate rounded bvalues used both for alpha selection and, "
            "when passed to plot_D0_vs_Delta.py, for choosing curves to draw."
        ),
    )
    ap.add_argument(
        "--roi-bvalmax",
        action="append",
        default=None,
        help=(
            "Per-ROI or per-ROI/direction bstep/bvalue override, repeatable. Use ROI=X or "
            "ROI:DIRECTION=X. Example: --roi-bvalmax Syringe=980 "
            "--roi-bvalmax fiber1:long=500. A direction-specific value takes precedence over "
            "the ROI value, then --bvalmax, then the highest bvalue."
        ),
    )
    ap.add_argument("--reference-D0", type=float, default=0.0032, help="Reference value used for alpha_macro.")
    ap.add_argument("--reference-D0-error", type=float, default=0.0000283512, help="Reference value error.")
    ap.add_argument(
        "--direction-alias",
        action="append",
        default=None,
        help=(
            "Optional legacy alias raw=grouped. Repeatable. No aliases are applied by default; "
            "modern masters already contain direct long/tra rows."
        ),
    )
    ap.add_argument("--out-summary", type=Path, default=Path("plots/summary_alpha_values.xlsx"), help="Output summary_alpha_values.xlsx")
    ap.add_argument(
        "--out-avg",
        type=Path,
        default=None,
        help="Optional output for the aggregated D vs Delta_app table. If omitted, it is not rewritten because it matches D_vs_delta_app.combined.",
    )
    ap.add_argument("--plot-rois", nargs="+", default=None, help="ROIs to show in the alpha_macro vs ROI plot. Defaults to --rois or all.")
    ap.add_argument("--plot-directions", nargs="+", default=None, help="Directions to show in the plot. Defaults to --dirs or all.")
    ap.add_argument("--out-plot", type=Path, default=None, help="Output PNG for alpha_macro vs ROI. Default: <out-summary-dir>/alpha_macro_vs_roi.png")
    args = ap.parse_args()

    if args.N is None and args.Hz is None and args.combined_table is None:
        args.N = 1.0

    df_avg = _load_measurements_from_args(args)
    if "roi" in df_avg.columns:
        df_avg["roi"] = df_avg["roi"].astype(str).str.strip().str.replace("_norm", "", regex=False)
    roi_bvalmax, roi_direction_bvalmax = _parse_roi_bvalmax(args.roi_bvalmax)
    if roi_bvalmax or roi_direction_bvalmax:
        present_rois = [str(x) for x in df_avg["roi"].dropna().astype(str).tolist()]
        roi_canon_to_actual: dict[str, str] = {}
        for roi_name in present_rois:
            key = _canonical_roi_name(roi_name)
            if key and key not in roi_canon_to_actual:
                roi_canon_to_actual[key] = roi_name

        present_directions = [str(x) for x in df_avg["direction"].dropna().astype(str).tolist()]
        direction_canon_to_actual: dict[str, str] = {}
        for direction_name in present_directions:
            key = direction_name.strip().lower()
            if key and key not in direction_canon_to_actual:
                direction_canon_to_actual[key] = direction_name

        roi_bvalmax_resolved: dict[str, float] = {}
        roi_direction_bvalmax_resolved: dict[tuple[str, str], float] = {}
        unknown_rois: list[str] = []
        for roi_name, bstep in roi_bvalmax.items():
            key = _canonical_roi_name(roi_name)
            actual = roi_canon_to_actual.get(key)
            if actual is None:
                unknown_rois.append(roi_name)
                continue
            roi_bvalmax_resolved[actual] = bstep

        unknown_selectors: list[str] = []
        for (roi_name, direction_name), bstep in roi_direction_bvalmax.items():
            actual_roi = roi_canon_to_actual.get(_canonical_roi_name(roi_name))
            actual_direction = direction_canon_to_actual.get(direction_name.strip().lower())
            if actual_roi is None or actual_direction is None:
                unknown_selectors.append(f"{roi_name}:{direction_name}")
                continue
            roi_direction_bvalmax_resolved[(actual_roi, actual_direction)] = bstep

        if unknown_rois:
            print(
                "[WARN] ROIs in --roi-bvalmax are not present in the filtered data; ignoring: "
                + ", ".join(sorted(unknown_rois))
            )
        if unknown_selectors:
            print(
                "[WARN] ROI/direction selectors in --roi-bvalmax are not present in the filtered data; ignoring: "
                + ", ".join(sorted(unknown_selectors))
            )
        roi_bvalmax = roi_bvalmax_resolved
        roi_direction_bvalmax = roi_direction_bvalmax_resolved
    aliases = parse_direction_aliases(args.direction_alias)
    df_avg, df_summary = compute_alpha_macro_summary(
        df_avg,
        reference_D0=float(args.reference_D0),
        reference_D0_error=float(args.reference_D0_error),
        selected_bstep=args.bvalmax,
        roi_selected_bsteps=roi_bvalmax or None,
        roi_direction_selected_bsteps=roi_direction_bvalmax or None,
        candidate_bsteps=args.plot_bsteps,
        candidate_bvalues=args.plot_bvalues,
        direction_aliases=aliases,
    )

    write_alpha_macro_outputs(
        df_avg,
        df_summary,
        out_summary_xlsx=args.out_summary,
        out_avg_xlsx=args.out_avg,
    )
    master_fit_params = _append_alpha_macro_fit_params(df_summary, args)
    master_alpha = _annotate_master_alpha_macro(df_summary, args)

    out_dir = args.out_summary.parent
    roi_order = args.plot_rois or args.rois or _ordered_unique(df_summary["roi"])
    directions = args.plot_directions or args.dirs or _ordered_unique(df_summary["direction"])
    out_plot = args.out_plot or (out_dir / "alpha_macro_vs_roi.png")
    _safe_plot_alpha_macro_vs_roi(
        df_summary,
        out_png=out_plot,
        roi_order=roi_order,
        directions=directions,
        subjs=args.subjs,
        title_prefix=rf"$\alpha_{{macro}}$ | N={args.N:g}" if args.N is not None else r"$\alpha_{macro}$",
    )

    print(f"[OK] summary: {args.out_summary}")
    if args.out_avg is not None:
        print(f"[OK] avg:     {args.out_avg}")
    if master_fit_params is not None:
        print(f"[OK] master fit params: {master_fit_params}")
    if master_alpha is not None:
        path, n_rows = master_alpha
        print(f"[OK] master alpha_macro: {path} ({n_rows} rows)")
    print(f"[OK] plot:    {out_plot}")


if __name__ == "__main__":
    main()
