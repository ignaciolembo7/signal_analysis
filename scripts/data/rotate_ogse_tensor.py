from __future__ import annotations

import repo_bootstrap  # noqa: F401

from pathlib import Path
import argparse
import pandas as pd

from data_processing.io import write_table_outputs
from data_processing.master_table import (
    append_master_rows,
    build_analysis_id_from_columns,
    load_master_table,
)
from fitting.cli_common import add_master_source_args
from pipeline.recipe import selected_rows_or_legacy_dataframe
from signal_rotation.rotation_tensor import rotate_signals_tensor


def _infer_exp_dir(df: pd.DataFrame, long_parquet: Path) -> str:
    if "sheet" in df.columns:
        vals = pd.Series(df["sheet"]).dropna().astype(str).str.strip().unique().tolist()
        if len(vals) == 1 and vals[0]:
            return vals[0]

    parent = long_parquet.parent.name
    if parent and parent != ".":
        return parent

    stem = long_parquet.stem.replace(".long", "")
    if "_ep2d" in stem:
        return stem.split("_ep2d")[0]
    return stem


def _analysis_id_from_columns(df: pd.DataFrame, *, row_kind: str) -> str:
    preferred = (
        "subj",
        "sheet",
        "type",
        "protocol",
        "group",
        "TN",
        "N",
        "td_ms",
        "Hz",
        "G",
        "sequence",
    )
    unique_cols = []
    for col in preferred:
        if col not in df.columns:
            continue
        values = pd.Series(df[col]).dropna().unique()
        if len(values) == 1:
            unique_cols.append(col)
    return build_analysis_id_from_columns(df, columns=unique_cols, prefix=row_kind)


def _drop_existing_rotation_rows(master: pd.DataFrame, rotated: pd.DataFrame) -> pd.DataFrame:
    """Remove previous rotated/dproj rows for the same acquisition group."""
    out = master.copy()
    group_cols = [c for c in ["subj", "sheet", "td_ms", "N", "Hz"] if c in rotated.columns and c in out.columns]
    if not group_cols:
        return out

    groups = rotated[group_cols].drop_duplicates()
    remove = out["row_kind"].astype(str).isin({"signal_rotated", "dproj"})
    for col in group_cols:
        values = groups[col].dropna().unique()
        if len(values) == 0:
            continue
        if pd.api.types.is_numeric_dtype(rotated[col]) or pd.api.types.is_numeric_dtype(out[col]):
            col_values = pd.to_numeric(out[col], errors="coerce")
            col_match = False
            for value in values:
                col_match = col_match | col_values.sub(float(value)).abs().le(1e-6)
            remove = remove & col_match
        else:
            remove = remove & out[col].astype(str).isin({str(v) for v in values})
    return out.loc[~remove].reset_index(drop=True)


def _normalize_key_value(value: object) -> str:
    if pd.isna(value):
        return ""
    try:
        num = float(value)
    except (TypeError, ValueError):
        return str(value)
    if abs(num - round(num)) < 1e-9:
        return str(int(round(num)))
    return f"{num:.12g}"


def _row_keys(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    normalized = df[cols].apply(lambda col: col.map(_normalize_key_value))
    return normalized.agg("\x1f".join, axis=1)


def _replace_signal_rows(master: pd.DataFrame | None, signal_rows: pd.DataFrame) -> pd.DataFrame | None:
    if master is None:
        return None
    key_cols = [
        c
        for c in ["subj", "sheet", "td_ms", "N", "Hz", "stat", "roi", "direction", "b_step"]
        if c in master.columns and c in signal_rows.columns
    ]
    if not key_cols:
        return master
    signal_keys = set(_row_keys(signal_rows, key_cols))
    master_keys = _row_keys(master, key_cols)
    remove = master["row_kind"].astype(str).eq("signal") & master_keys.isin(signal_keys)
    return master.loc[~remove].reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("long_parquet", type=Path, nargs="?", help="Legacy input clean signal .long.parquet.")
    ap.add_argument("--out_dir", type=Path, default=Path("analysis/ogse_experiments/data-rotated/tables"))
    add_master_source_args(ap, default_row_kind="signal")
    ap.add_argument(
        "--no-legacy-output",
        action="store_true",
        help="Do not write legacy *.rot_tensor.long.parquet and *.Dproj.long.parquet outputs.",
    )
    ap.add_argument("--solver", type=str, default="lstsq", choices=["lstsq", "solve"])
    ap.add_argument("--s0_mode", type=str, default="dir1", choices=["dir1", "mean"])
    ap.add_argument("--b_col", type=str, default="bvalue")
    ap.add_argument(
        "--g_type",
        type=str,
        default="g",
        choices=["g", "g_lin_max", "g_thorsten"],
        help=(
            "Gradient basis used to compute the b-value that drives the tensor fit "
            "(and that gets stamped on the rotated/derived directions' 'bvalue' column). "
            "'g' uses the actual measured gradient per b_step (matches the original signal "
            "rows' bvalue exactly); 'g_lin_max'/'g_thorsten' substitute an idealized gradient "
            "ramp, which only coincides with the real bvalue if the acquisition itself used "
            "a linear ramp."
        ),
    )
    ap.add_argument(
        "--tra_axes",
        type=str,
        default="y,z",
        help="Comma-separated axes averaged to build the 'tra' direction, e.g. 'x,y'.",
    )
    ap.add_argument(
        "--long_axes",
        type=str,
        default="x",
        help="Comma-separated axes averaged to build the 'long' direction, e.g. 'z'.",
    )
    ap.add_argument(
        "--dirs_txt",
        type=Path,
        default=None,
        help="Nx3 TXT with directions and no header. Defaults to assets/dirs/dirs_{ndirs}.txt.",
    )
    ap.add_argument("--dirs_csv", type=Path, default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.dirs_txt is not None and args.dirs_csv is not None:
        raise SystemExit("Use only one of --dirs_txt or --dirs_csv.")

    selected = selected_rows_or_legacy_dataframe(
        args,
        legacy_path=args.long_parquet,
        default_row_kind=str(args.row_kind or "signal"),
        signal_rotated=False,
    )
    df = selected.df
    if df is None:
        raise SystemExit("No input rows were selected.")
    input_kind = selected.source
    dirs_file = args.dirs_txt if args.dirs_txt is not None else args.dirs_csv
    tra_axes = tuple(a.strip().lower() for a in args.tra_axes.split(",") if a.strip())
    long_axes = tuple(a.strip().lower() for a in args.long_axes.split(",") if a.strip())

    res = rotate_signals_tensor(
        df,
        solver=args.solver,
        s0_mode=args.s0_mode,
        b_col=args.b_col,
        g_type=args.g_type,
        dirs_file=dirs_file,
        tra_axes=tra_axes,
        long_axes=long_axes,
    )

    if args.master_parquet is not None:
        master = load_master_table(args.master_parquet) if args.master_parquet.exists() else None
        if master is not None:
            master = _replace_signal_rows(master, res.signal_long_with_dproj)
            master = _drop_existing_rotation_rows(master, res.rotated_signal_long)
        signal_analysis_id = _analysis_id_from_columns(res.signal_long_with_dproj, row_kind="signal")
        master = append_master_rows(
            master,
            res.signal_long_with_dproj,
            row_kind="signal",
            analysis_id=signal_analysis_id,
        )
        print("Updated original signal D_proj in master:", args.master_parquet)
        rotated_analysis_id = _analysis_id_from_columns(res.rotated_signal_long, row_kind="signal_rotated")
        append_master_rows(
            master,
            res.rotated_signal_long,
            row_kind="signal_rotated",
            analysis_id=rotated_analysis_id,
            out_path=args.master_parquet,
        )
        print("Appended rotated signals to master:", args.master_parquet)

    if not args.no_legacy_output:
        if args.long_parquet is not None:
            exp_dir = args.out_dir / _infer_exp_dir(df, args.long_parquet)
            stem = args.long_parquet.stem.replace(".long", "")
        else:
            exp_label = _analysis_id_from_columns(df, row_kind=input_kind)
            exp_dir = args.out_dir / "master"
            stem = exp_label
        exp_dir.mkdir(parents=True, exist_ok=True)

        out_rot = exp_dir / f"{stem}.rot_tensor.long.parquet"
        out_dpr = exp_dir / f"{stem}.rot_tensor.Dproj.long.parquet"

        write_table_outputs(res.rotated_signal_long, out_rot, xlsx_path=out_rot.with_suffix(".xlsx"))
        write_table_outputs(res.dproj_long, out_dpr, xlsx_path=out_dpr.with_suffix(".xlsx"))

        print("Saved rotated signals:", out_rot)
        print("Saved Dproj:", out_dpr)


if __name__ == "__main__":
    main()
