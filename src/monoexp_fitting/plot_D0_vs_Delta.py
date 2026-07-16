from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd

from tc_fittings.alpha_macro_summary import (
    load_dproj_measurements,
    plot_d_vs_delta_curves,
)


def _pick_column(df: pd.DataFrame, candidates: Sequence[str]) -> str | None:
    lower_to_name = {str(col).strip().lower(): str(col) for col in df.columns}
    for candidate in candidates:
        col = lower_to_name.get(candidate.lower())
        if col is not None:
            return col
    return None


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported summary format: {path}")


def _norm_roi(value: object) -> str:
    return str(value).strip().replace("_norm", "").lower()


def load_selected_bstep_map(summary_alpha: str | Path) -> dict[tuple[str, str, str], int]:
    path = Path(summary_alpha)
    if not path.exists():
        raise FileNotFoundError(path)

    df = _read_table(path)
    subj_col = _pick_column(df, ["subj", "brain"])
    roi_col = _pick_column(df, ["roi", "region"])
    direction_col = _pick_column(df, ["direction", "direccion"])
    bstep_col = _pick_column(df, ["selected_bstep"])
    if subj_col is None or roi_col is None or direction_col is None:
        raise KeyError(
            f"summary_alpha missing key columns in {path}. Expected subj/roi/direction aliases, got {list(df.columns)}"
        )
    if bstep_col is None:
        return {}

    extra_cols = [c for c in ["direction_kind", "direction_components"] if c in df.columns]
    out = df[[subj_col, roi_col, direction_col, bstep_col] + extra_cols].copy()
    out[subj_col] = out[subj_col].astype(str).str.strip()
    out[roi_col] = out[roi_col].astype(str).str.strip()
    out[direction_col] = out[direction_col].astype(str).str.strip()
    out[bstep_col] = pd.to_numeric(out[bstep_col], errors="coerce")
    out = out.dropna(subset=[bstep_col])
    out = out[out[bstep_col] >= 1]
    if out.empty:
        return {}

    if "direction_kind" in out.columns:
        out["_pref"] = out["direction_kind"].astype(str).str.lower().map(lambda v: 0 if v == "derived" else 1)
    else:
        out["_pref"] = 0

    out = out.sort_values([subj_col, roi_col, direction_col, "_pref"], kind="stable")
    out = out.drop_duplicates(subset=[subj_col, roi_col, direction_col], keep="first")

    mapping: dict[tuple[str, str, str], int] = {}
    for _, row in out.iterrows():
        subj = str(row[subj_col])
        roi = str(row[roi_col])
        direction = str(row[direction_col])
        bstep = int(round(float(row[bstep_col])))
        mapping[(subj, roi, direction)] = bstep
        mapping[(subj, _norm_roi(roi), direction)] = bstep
        if "direction_components" in out.columns:
            for component in str(row.get("direction_components", "")).split("|"):
                component = component.strip()
                if component:
                    mapping[(subj, roi, component)] = bstep
                    mapping[(subj, _norm_roi(roi), component)] = bstep
    return mapping


def load_all_measurements(
    dproj_root: str | Path,
    *,
    pattern: str = "**/*.Dproj.long.parquet",
    dirs: Sequence[str] | None = None,
    rois: Sequence[str] | None = None,
    subjs: Sequence[str] | None = None,
    N: float | None = None,
    Hz: float | None = None,
    bvalue_decimals: int = 1,
) -> pd.DataFrame:
    return load_dproj_measurements(
        dproj_root,
        pattern=pattern,
        subjs=subjs,
        rois=rois,
        directions=dirs,
        N=N,
        Hz=Hz,
        bvalue_decimals=bvalue_decimals,
    )


def plot_all_groups(
    df_avg: pd.DataFrame,
    *,
    out_dir: str | Path,
    selected_bstep: int | None = None,
    selected_bstep_by_group: Mapping[tuple[str, str, str], int] | None = None,
    reference_D0: float | None = None,
    reference_D0_error: float | None = None,
) -> list[Path]:
    return plot_d_vs_delta_curves(
        df_avg,
        out_dir=out_dir,
        selected_bstep=selected_bstep,
        selected_bstep_by_group=selected_bstep_by_group,
        reference_D0=reference_D0,
        reference_D0_error=reference_D0_error,
    )
