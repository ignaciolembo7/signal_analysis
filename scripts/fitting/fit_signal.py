from __future__ import annotations

import argparse
from pathlib import Path

import repo_bootstrap  # noqa: F401

import numpy as np
import pandas as pd

from data_processing.master_table import append_master_rows, filter_master_rows, load_master_table
from fitting.cli_common import add_master_source_args, load_master_input, split_cli_values
from monoexp_fitting.fit_signal import (
    AUTO_FIT_ERR_FLOOR,
    AUTO_FIT_MAX_POINTS,
    AUTO_FIT_MIN_POINTS,
    AUTO_FIT_REL_TOL,
    fit_signal_monoexp,
    write_fit_outputs,
)
from ogse_fitting.fit_ogse_free_signal import fit_ogse_free_signal


def _split_optional(values: list[str] | None) -> list[str] | None:
    values = split_cli_values(values)
    if values is None or (len(values) == 1 and values[0].upper() == "ALL"):
        return None
    return [str(v) for v in values]


def _filter_close(df: pd.DataFrame, col: str, value: object, *, atol: float = 1e-6) -> pd.DataFrame:
    if col not in df.columns or pd.isna(value):
        return df
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return df[df[col].astype(str).str.strip() == str(value).strip()].copy()
    vals = pd.to_numeric(df[col], errors="coerce")
    return df[np.isclose(vals.to_numpy(dtype=float), numeric, atol=atol, equal_nan=False)].copy()


def _filter_text(df: pd.DataFrame, col: str, value: object) -> pd.DataFrame:
    if col not in df.columns or pd.isna(value):
        return df
    text = str(value).strip()
    if not text or text.upper() == "ALL":
        return df
    return df[df[col].astype(str).str.strip() == text].copy()


def _selected_from_manifest(
    master: pd.DataFrame,
    manifest_row: pd.Series,
    *,
    default_row_kind: str,
    requested_model: str,
) -> pd.DataFrame:
    model = str(manifest_row.get("model", "monoexp") or "monoexp").strip()
    if model and model.lower() not in {requested_model.lower(), "all"}:
        raise ValueError(
            f"Signal-fit manifest requests model={model!r}, but the command requests {requested_model!r}."
        )

    selectors: dict[str, object] = {"row_kind": default_row_kind}
    for col in ["subj", "sheet", "roi", "direction"]:
        value = manifest_row.get(col)
        if value is not None and not pd.isna(value) and str(value).strip().upper() != "ALL":
            selectors[col] = str(value).strip()
    selected = filter_master_rows(master, **selectors)
    for col in ["td_ms", "N", "Hz"]:
        selected = _filter_close(selected, col, manifest_row.get(col))
    return selected.reset_index(drop=True)


def _run_one_selection(
    df: pd.DataFrame,
    *,
    args: argparse.Namespace,
    directions: list[str] | None = None,
    rois: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    common = dict(
        df=df,
        directions=directions,
        rois=rois,
        ycol=args.ycol,
        free_M0=bool(args.free_M0),
        fix_M0=float(args.fix_M0),
        stat_keep=str(args.stat or "avg"),
        correction_factor_col="grad_correction_factor" if bool(args.apply_grad_corr) else None,
        outdir_plots=Path(args.plot_dir) if args.plot_dir is not None else None,
    )
    if args.model == "ogse_free":
        outputs = fit_ogse_free_signal(
            **common,
            g_axis=args.b_axis,
            D0_init_mm2_s=float(args.D0_init),
        )
    else:
        outputs = fit_signal_monoexp(
            **common,
            b_axis=args.b_axis,
            fit_points=args.fit_points,
            auto_fit_points=bool(args.auto_fit_points),
            auto_fit_min_points=int(args.auto_fit_min_points),
            auto_fit_max_points=args.auto_fit_max_points,
            auto_fit_rel_tol=float(args.auto_fit_tol),
            auto_fit_err_floor=float(args.auto_fit_err_floor),
            D0_init=float(args.D0_init),
            gamma=float(args.gamma),
        )
    return outputs.fit_params, outputs.fit_points


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Fit signal curves with a monoexponential or OGSE free-diffusion model. "
            "Monoexponential automatic mode tests leading b_step prefixes."
        )
    )
    ap.add_argument("table", type=Path, nargs="?", help="Optional legacy signal table parquet/xlsx/csv.")
    ap.add_argument("--manifest", type=Path, default=None, help="Optional signal_fits.csv manifest.")
    ap.add_argument("--model", default="monoexp", choices=["monoexp", "ogse_free"])
    ap.add_argument("--out-root", "--out_root", dest="out_root", type=Path, required=True)
    ap.add_argument("--plot-dir", "--plot_dir", dest="plot_dir", type=Path, default=None)
    ap.add_argument("--ycol", default="value_norm", choices=["value", "value_norm"])
    ap.add_argument("--b-axis", "--b_axis", "--g_type", dest="b_axis", default="bvalue_g")
    ap.add_argument("--gamma", type=float, default=267.5221900)
    ap.add_argument("--stat", default="avg")

    fit_group = ap.add_mutually_exclusive_group()
    fit_group.add_argument("--fit-points", "--fit_points", dest="fit_points", type=int, default=None)
    fit_group.add_argument("--auto-fit-points", "--auto_fit_points", dest="auto_fit_points", action="store_true")
    ap.add_argument(
        "--auto-fit-tol", "--auto_fit_tol", dest="auto_fit_tol", type=float,
        default=AUTO_FIT_REL_TOL,
        help="Maximum accepted rmse_log for automatic prefix selection (default: 0.05).",
    )
    ap.add_argument(
        "--auto-fit-err-floor",
        "--auto_fit_err_floor",
        dest="auto_fit_err_floor",
        type=float,
        default=AUTO_FIT_ERR_FLOOR,
    )
    ap.add_argument(
        "--auto-fit-min-points",
        "--auto_fit_min_points",
        dest="auto_fit_min_points",
        type=int,
        default=AUTO_FIT_MIN_POINTS,
    )
    ap.add_argument(
        "--auto-fit-max-points",
        "--auto_fit_max_points",
        dest="auto_fit_max_points",
        type=int,
        default=AUTO_FIT_MAX_POINTS,
    )
    ap.add_argument("--D0-init", "--D0_init", dest="D0_init", type=float, default=0.0023)

    corr_group = ap.add_mutually_exclusive_group()
    corr_group.add_argument("--apply-grad-corr", "--apply_grad_corr", dest="apply_grad_corr", action="store_true")
    corr_group.add_argument("--no-grad-corr", "--no_grad_corr", dest="no_grad_corr", action="store_true")

    m0_group = ap.add_mutually_exclusive_group()
    m0_group.add_argument("--fix-M0", "--fix_M0", dest="fix_M0", type=float, default=1.0)
    m0_group.add_argument("--free-M0", "--free_M0", dest="free_M0", action="store_true")

    ap.add_argument("--directions", nargs="+", default=None)
    ap.add_argument("--rois", nargs="+", default=None)
    ap.add_argument("--append-fit-params-to-master", action="store_true")
    ap.add_argument("--append-fit-points-to-master", action="store_true")
    add_master_source_args(ap, default_row_kind="signal_rotated", include_stat=False, include_td_ms=False, include_N=False)
    args = ap.parse_args()

    if args.fit_points is not None and args.fit_points <= 0:
        raise ValueError("--fit-points must be > 0.")
    if args.auto_fit_tol < 0:
        raise ValueError("--auto-fit-tol must be >= 0.")
    if args.auto_fit_err_floor < 0:
        raise ValueError("--auto-fit-err-floor must be >= 0.")
    if args.auto_fit_min_points < 1:
        raise ValueError("--auto-fit-min-points must be >= 1.")
    if args.auto_fit_max_points is not None and args.auto_fit_max_points < args.auto_fit_min_points:
        raise ValueError("--auto-fit-max-points must be >= --auto-fit-min-points.")

    params_tables: list[pd.DataFrame] = []
    points_tables: list[pd.DataFrame] = []

    if args.master_parquet is not None:
        master = load_master_table(args.master_parquet)
        default_row_kind = str(args.row_kind or "signal_rotated")
        if args.manifest is not None and args.manifest.exists():
            manifest = pd.read_csv(args.manifest, comment="#")
            if manifest.empty:
                raise ValueError(f"Manifest is empty: {args.manifest}")
            for _, row in manifest.iterrows():
                selected = _selected_from_manifest(
                    master,
                    row,
                    default_row_kind=default_row_kind,
                    requested_model=str(args.model),
                )
                if selected.empty:
                    print(f"WARNING: no rows matched signal-fit manifest row:\n{row.to_string()}")
                    continue
                dirs = [str(row["direction"]).strip()] if "direction" in row and str(row.get("direction", "")).strip().upper() != "ALL" else None
                rois = [str(row["roi"]).strip()] if "roi" in row and str(row.get("roi", "")).strip().upper() != "ALL" else None
                fp, ft = _run_one_selection(selected, args=args, directions=dirs, rois=rois)
                params_tables.append(fp)
                points_tables.append(ft)
        else:
            selected = load_master_input(args, default_row_kind=default_row_kind, signal_rotated=default_row_kind == "signal_rotated")
            fp, ft = _run_one_selection(
                selected,
                args=args,
                directions=_split_optional(args.directions),
                rois=_split_optional(args.rois),
            )
            params_tables.append(fp)
            points_tables.append(ft)
    else:
        if args.table is None:
            raise ValueError("Pass --master-parquet or a legacy table path.")
        suffix = args.table.suffix.lower()
        if suffix == ".parquet":
            df = pd.read_parquet(args.table)
        elif suffix in {".xlsx", ".xls"}:
            df = pd.read_excel(args.table, sheet_name=0)
        elif suffix == ".csv":
            df = pd.read_csv(args.table)
        else:
            raise ValueError(f"Unsupported input table: {args.table}")
        fp, ft = _run_one_selection(df, args=args, directions=_split_optional(args.directions), rois=_split_optional(args.rois))
        params_tables.append(fp)
        points_tables.append(ft)

    fit_params = pd.concat([t for t in params_tables if not t.empty], ignore_index=True, sort=False) if params_tables else pd.DataFrame()
    fit_points = pd.concat([t for t in points_tables if not t.empty], ignore_index=True, sort=False) if points_tables else pd.DataFrame()

    from monoexp_fitting.fit_signal import FitOutputs

    params_path, points_path = write_fit_outputs(
        FitOutputs(fit_params=fit_params, fit_points=fit_points),
        Path(args.out_root),
        b_axis=str(args.b_axis),
        ycol=str(args.ycol),
        model=str(args.model),
    )
    if params_path is not None:
        print("OK fit_params:", params_path)
    else:
        print("WARNING: no fit_params rows were produced.")
    if points_path is not None:
        print("OK fit_points:", points_path)

    if args.master_parquet is not None and args.append_fit_params_to_master and not fit_params.empty:
        append_master_rows(args.master_parquet, fit_params, row_kind="fit_params", out_path=args.master_parquet)
        print("Appended fit_params to master:", args.master_parquet)
    if args.master_parquet is not None and args.append_fit_points_to_master and not fit_points.empty:
        append_master_rows(args.master_parquet, fit_points, row_kind="fit_points", out_path=args.master_parquet)
        print("Appended fit_points to master:", args.master_parquet)


if __name__ == "__main__":
    main()
