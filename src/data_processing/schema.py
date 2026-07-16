from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

from fitting.b_from_g import b_from_g
from tools.strict_columns import raise_on_unrecognized_column_names

CLEAN_SIGNAL_LONG_COLUMNS = [
    "stat","roi","direction","b_step",
    "bvalue","bvalue_g","bvalue_g_lin_max","bvalue_thorsten","bmax",
    "g","g_max","g_lin_max","g_thorsten",
    "gradient_axis_kind",
    "one_g_per_sequence",
    "grad_correction_factor",
    "alpha_macro",
    "value","value_norm","S0","D_proj",
    "TE","TR","td_ms","Hz",
    "max_dur_ms","Delta_app_ms","delta_ms","tm_ms", 
    "TN","N","G","x", "y",
    "subj","type","protocol","sequence","group","sheet","source_file",
]

CLEAN_DPROJ_LONG_COLUMNS = [
    "roi", "direction", "b_step",
    "bvalue", "bvalue_g", "bvalue_g_lin_max", "bvalue_thorsten", "bmax",
    "g", "g_max", "g_lin_max", "g_thorsten",
    "gradient_axis_kind",
    "one_g_per_sequence",
    "D_proj",
    "TE","TR","td_ms","Hz", 
    "max_dur_ms","Delta_app_ms","delta_ms","tm_ms",
    "TN","N","G","x","y",
    "subj","type","protocol","sequence","group","sheet","source_file",
]

SIGNAL_LONG_PREFIX = [
    "stat","roi","direction","b_step","bvalue",
    "bvalue_orig","bvalue_g","bvalue_g_lin_max","bvalue_thorsten","bmax",
    "g","g_max","g_lin_max","g_thorsten",
    "gradient_axis_kind",
    "one_g_per_sequence",
    "grad_correction_factor",
    "alpha_macro",
    "value","value_norm","S0","D_proj",
    "TE","TR","td_ms","Hz",
    "max_dur_ms","Delta_app_ms","delta_ms","tm_ms", 
    "TN","N","G","x","y",
    "subj","type","sequence","group","protocol","sheet","source_file",
]

DPROJ_LONG_PREFIX = [
    "roi","direction","b_step","bvalue",
    "bvalue_orig","bvalue_g","bvalue_g_lin_max","bvalue_thorsten","bmax",
    "g","g_max","g_lin_max","g_thorsten",
    "gradient_axis_kind",
    "one_g_per_sequence",
    "D_proj",
    "TE","TR","td_ms","Hz",
    "max_dur_ms","Delta_app_ms","delta_ms","tm_ms",
    "TN","N","G","x","y",
    "subj","type","sequence","group","protocol","sheet","source_file",
]

RENAME_MAP = {
    "bvalues": "bvalue",
    "bval": "bvalue",
    "b": "bvalue",
    "Signal": "value",
    "signal": "value",
}

def _unique_scalar(df: pd.DataFrame, col: str) -> float | None:
    if col not in df.columns:
        return None
    u = pd.Series(df[col]).dropna().unique()
    if len(u) != 1:
        return None
    val = pd.to_numeric(u[0], errors="coerce")
    return float(val) if pd.notna(val) else None

def _apply_renames(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename = {k: v for k, v in RENAME_MAP.items() if k in out.columns and v not in out.columns}
    if rename:
        out = out.rename(columns=rename)
    return out

def _append_sorted(out_cols: list[str], candidates: Iterable[str]) -> list[str]:
    cands = [c for c in candidates if c not in out_cols]
    out_cols.extend(sorted(cands))
    return out_cols

def _unique_scalar_from_aliases(df: pd.DataFrame, aliases: list[str]) -> float | None:
    for col in aliases:
        v = _unique_scalar(df, col)
        if v is not None:
            return v
    return None

def _infer_params_for_b(df: pd.DataFrame) -> tuple[float | None, float | None, float | None, float | None]:
    N = _unique_scalar(df, "param_N")
    gamma = _unique_scalar(df, "meta_gamma")
    delta_ms = _unique_scalar(df, "param_delta_ms")
    delta_app_ms = _unique_scalar(df, "param_delta_app_ms")
    return N, gamma, delta_ms, delta_app_ms

def _infer_clean_params_for_b(df: pd.DataFrame) -> tuple[float | None, float | None, float | None, float | None]:
    N = _unique_scalar_from_aliases(df, ["N", "param_N"])
    gamma = _unique_scalar_from_aliases(df, ["gamma", "meta_gamma"])
    delta_ms = _unique_scalar_from_aliases(df, ["delta_ms", "param_delta_ms"])
    delta_app_ms = _unique_scalar_from_aliases(df, ["Delta_app_ms", "delta_app_ms", "param_delta_app_ms"])
    return N, gamma, delta_ms, delta_app_ms

def _coalesce_col(out: pd.DataFrame, target: str, aliases: list[str]) -> None:
    if target in out.columns:
        return
    for c in aliases:
        if c in out.columns:
            out[target] = out[c]
            return

def _ensure_timing_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    _coalesce_col(out, "max_dur_ms", ["param_max_dur_ms", "param_d_ms", "d_ms", "param_d"])
    _coalesce_col(out, "tm_ms", ["param_tm_ms", "param_TM_ms", "TM_ms", "tm_ms"])

    td_ok = ("td_ms" in out.columns) and not pd.to_numeric(out["td_ms"], errors="coerce").dropna().empty
    if not td_ok and ("max_dur_ms" in out.columns) and ("tm_ms" in out.columns):
        d = pd.to_numeric(out["max_dur_ms"], errors="coerce")
        tm = pd.to_numeric(out["tm_ms"], errors="coerce")
        out["td_ms"] = 2.0 * d + tm

    return out

def _ensure_S0_and_norm(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    if "S0" not in out.columns:
        out["S0"] = pd.NA

    have_keys = all(c in out.columns for c in ["stat","roi","direction","b_step","value"])
    if have_keys:
        s0 = (
            out.loc[out["b_step"] == 0]
            .groupby(["stat","roi","direction"], as_index=False)["value"]
            .mean()
            .rename(columns={"value": "S0"})
        )
        out = out.drop(columns=["S0"]).merge(s0, on=["stat","roi","direction"], how="left")

    if "value_norm" not in out.columns:
        out["value_norm"] = pd.NA

    if "value" in out.columns and "S0" in out.columns:
        v = pd.to_numeric(out["value"], errors="coerce")
        s0v = pd.to_numeric(out["S0"], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            out["value_norm"] = v / s0v
        if "b_step" in out.columns:
            mask0 = (pd.to_numeric(out["b_step"], errors="coerce") == 0) & s0v.notna()
            out.loc[mask0, "value_norm"] = 1.0

    return out

def _ensure_bvalue_derivatives(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    gradient_axis_kind = None
    if "gradient_axis_kind" in out.columns:
        kinds = pd.Series(out["gradient_axis_kind"]).dropna().astype(str).str.lower().unique().tolist()
        if len(kinds) == 1:
            gradient_axis_kind = kinds[0]

    if "bvalue" not in out.columns:
        raise KeyError("Missing 'bvalue' column (after renames).")

    if "bvalue_orig" not in out.columns:
        out["bvalue_orig"] = out["bvalue"]

    for c in ["bvalue_g","bvalue_g_lin_max","bvalue_thorsten"]:
        if c not in out.columns:
            out[c] = pd.NA

    if gradient_axis_kind == "g":
        return out

    N, gamma, delta_ms, delta_app_ms = _infer_clean_params_for_b(out)
    if N is None or gamma is None or delta_ms is None or delta_app_ms is None:
        return out

    if "g" in out.columns and out["bvalue_g"].isna().any():
        g = pd.to_numeric(out["g"], errors="coerce").to_numpy(float)
        out["bvalue_g"] = b_from_g(g, N=N, gamma=gamma, delta_ms=delta_ms, delta_app_ms=delta_app_ms, g_type="g")

    if "g_lin_max" in out.columns and out["bvalue_g_lin_max"].isna().any():
        g = pd.to_numeric(out["g_lin_max"], errors="coerce").to_numpy(float)
        out["bvalue_g_lin_max"] = b_from_g(g, N=N, gamma=gamma, delta_ms=delta_ms, delta_app_ms=delta_app_ms, g_type="g_lin_max")

    if "g_thorsten" in out.columns and out["bvalue_thorsten"].isna().any():
        g = pd.to_numeric(out["g_thorsten"], errors="coerce").to_numpy(float)
        out["bvalue_thorsten"] = b_from_g(g, N=N, gamma=gamma, delta_ms=delta_ms, delta_app_ms=delta_app_ms, g_type="g_thorsten")

    return out

def finalize_clean_signal_long(df: pd.DataFrame) -> pd.DataFrame:
    out = _apply_renames(df)
    raise_on_unrecognized_column_names(out.columns, context="finalize_clean_signal_long")
    required = {"stat", "roi", "direction", "b_step", "bvalue", "value"}
    missing = required - set(out.columns)
    if missing:
        raise ValueError(f"Not a valid clean signal-long table, missing: {sorted(missing)}")

    if "source_file" not in out.columns:
        out["source_file"] = ""

    out = _ensure_S0_and_norm(out)
    out = _ensure_bvalue_derivatives(out)

    for c in CLEAN_SIGNAL_LONG_COLUMNS:
        if c not in out.columns:
            out[c] = pd.NA

    return out[CLEAN_SIGNAL_LONG_COLUMNS].copy()

def finalize_clean_dproj_long(df: pd.DataFrame) -> pd.DataFrame:
    out = _apply_renames(df)
    raise_on_unrecognized_column_names(out.columns, context="finalize_clean_dproj_long")
    required = {"roi", "direction", "b_step", "bvalue", "D_proj"}
    missing = required - set(out.columns)
    if missing:
        raise ValueError(f"Not a valid clean dproj-long table, missing: {sorted(missing)}")

    if "source_file" not in out.columns:
        out["source_file"] = ""
    out = _ensure_bvalue_derivatives(out)

    for c in CLEAN_DPROJ_LONG_COLUMNS:
        if c not in out.columns:
            out[c] = pd.NA

    return out[CLEAN_DPROJ_LONG_COLUMNS].copy()

def ensure_signal_long_schema(df: pd.DataFrame) -> pd.DataFrame:
    out = _apply_renames(df)
    raise_on_unrecognized_column_names(out.columns, context="ensure_signal_long_schema")
    required = {"stat","roi","direction","b_step","bvalue","value"}
    missing = required - set(out.columns)
    if missing:
        raise ValueError(f"Not a valid signal-long table, missing: {sorted(missing)}")

    if "source_file" not in out.columns:
        out["source_file"] = ""

    out = _ensure_timing_columns(out)
    out = _ensure_S0_and_norm(out)
    out = _ensure_bvalue_derivatives(out)

    for c in SIGNAL_LONG_PREFIX:
        if c not in out.columns:
            out[c] = pd.NA

    cols: list[str] = [c for c in SIGNAL_LONG_PREFIX if c in out.columns]
    cols = _append_sorted(cols, [c for c in out.columns if c.startswith("param_")])
    cols = _append_sorted(cols, [c for c in out.columns if c.startswith("meta_")])
    cols = _append_sorted(cols, out.columns)

    return out[cols]

def ensure_dproj_long_schema(df: pd.DataFrame) -> pd.DataFrame:
    out = _apply_renames(df)
    raise_on_unrecognized_column_names(out.columns, context="ensure_dproj_long_schema")
    required = {"roi","direction","b_step","bvalue","D_proj"}
    missing = required - set(out.columns)
    if missing:
        raise ValueError(f"Not a valid dproj-long table, missing: {sorted(missing)}")

    if "source_file" not in out.columns:
        out["source_file"] = ""

    out = _ensure_timing_columns(out)
    out = _ensure_bvalue_derivatives(out)

    for c in DPROJ_LONG_PREFIX:
        if c not in out.columns:
            out[c] = pd.NA

    cols: list[str] = [c for c in DPROJ_LONG_PREFIX if c in out.columns]
    cols = _append_sorted(cols, [c for c in out.columns if c.startswith("param_")])
    cols = _append_sorted(cols, [c for c in out.columns if c.startswith("meta_")])
    cols = _append_sorted(cols, out.columns)

    return out[cols]
