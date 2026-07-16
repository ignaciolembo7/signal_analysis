from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


VALID_G_TYPES = {"g", "g_max", "g_lin_max", "g_thorsten"}
VALID_AXIS_BASES = {
    "bvalue",
    "bvalue_g",
    "bvalue_g_lin_max",
    "bvalue_thorsten",
    *VALID_G_TYPES,
}

_AXIS_TO_GRADIENT_BASE = {
    "bvalue": "g",
    "g": "g",
    "bvalue_g": "g",
    "g_max": "g_max",
    "g_lin_max": "g_lin_max",
    "bvalue_g_lin_max": "g_lin_max",
    "g_thorsten": "g_thorsten",
    "bvalue_thorsten": "g_thorsten",
}

_AXIS_TO_BVALUE_COLUMN = {
    "bvalue": "bvalue",
    "bvalue_g": "bvalue_g",
    "bvalue_g_lin_max": "bvalue_g_lin_max",
    "bvalue_thorsten": "bvalue_thorsten",
}


@dataclass(frozen=True)
class AxisBundle:
    requested_axis: str
    axis_base: str
    gradient_base: str
    side: int | None
    gradient_raw: np.ndarray
    gradient_corr: np.ndarray
    bvalue_raw: np.ndarray | None
    bvalue_corr: np.ndarray | None
    axis_raw: np.ndarray
    axis_corr: np.ndarray


def b_from_g(
    g_mTm: np.ndarray,
    *,
    N: float,
    gamma: float,
    delta_ms: float,
    delta_app_ms: float,
    g_type: str,
) -> np.ndarray:
    """
    Compute b = N * gamma^2 * delta^2 * delta_app * g^2 / 1e9.

    Expected units:
      gamma: 1/(ms*mT)
      delta_ms, delta_app_ms: ms
      g_mTm: mT/m
      b: s/mm^2
    """
    g = np.asarray(g_mTm, dtype=float)
    if g_type not in VALID_G_TYPES:
        raise ValueError(f"b_from_g: unrecognized g_type {g_type!r}. Allowed values: {sorted(VALID_G_TYPES)}.")
    return N * (gamma**2) * (delta_ms**2) * delta_app_ms * (g**2) / 1e9


def normalize_axis_base(axis: str) -> str:
    base = str(axis).strip()
    if base.endswith("_1") or base.endswith("_2"):
        base = base[:-2]
    if base not in VALID_AXIS_BASES:
        raise ValueError(f"Unrecognized axis {axis!r}. Allowed values: {sorted(VALID_AXIS_BASES)}.")
    return base


def split_axis_side(axis: str) -> tuple[str, int | None]:
    raw = str(axis).strip()
    if raw.endswith("_1") or raw.endswith("_2"):
        return normalize_axis_base(raw[:-2]), int(raw[-1])
    return normalize_axis_base(raw), None


def axis_uses_bvalue(axis: str) -> bool:
    return normalize_axis_base(axis).startswith("bvalue")


def gradient_base_for_axis(axis: str) -> str:
    base = normalize_axis_base(axis)
    return _AXIS_TO_GRADIENT_BASE[base]


def axes_share_gradient_family(axis_a: str, axis_b: str) -> bool:
    return gradient_base_for_axis(axis_a) == gradient_base_for_axis(axis_b)


def default_plot_axis_for_fit(axis: str, *, side: int | None = None) -> str:
    base = normalize_axis_base(axis)
    if side is None:
        return base
    return f"{base}_{int(side)}"


def coerce_positive_correction_factor(value: float | None) -> float:
    if value is None:
        return 1.0
    out = float(value)
    if not np.isfinite(out) or out <= 0.0:
        return 1.0
    return out


def gradient_column_name(axis: str, *, side: int | None = None) -> str:
    base = gradient_base_for_axis(axis)
    return base if side is None else f"{base}_{int(side)}"


def bvalue_column_name(axis: str, *, side: int | None = None) -> str | None:
    base = normalize_axis_base(axis)
    col = _AXIS_TO_BVALUE_COLUMN.get(base)
    if col is None:
        return None
    return col if side is None else f"{col}_{int(side)}"


def extract_gradient_array(
    df: pd.DataFrame,
    *,
    axis: str,
    side: int | None = None,
) -> np.ndarray:
    col = gradient_column_name(axis, side=side)
    if col not in df.columns:
        raise KeyError(f"Missing required gradient column {col!r}. Columns={list(df.columns)}")
    arr = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
    return arr


def _extract_stored_bvalue_array(
    df: pd.DataFrame,
    *,
    axis: str,
    side: int | None = None,
) -> np.ndarray | None:
    col = bvalue_column_name(axis, side=side)
    if col is None or col not in df.columns:
        return None
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def bvalue_from_gradient(
    gradient_values: np.ndarray,
    *,
    axis: str,
    N: float | None,
    gamma: float,
    delta_ms: float | None,
    Delta_app_ms: float | None,
) -> np.ndarray:
    if N is None or delta_ms is None or Delta_app_ms is None:
        raise ValueError(
            "Computing a b-value axis from corrected gradients requires N, delta_ms, and Delta_app_ms."
        )
    return b_from_g(
        np.asarray(gradient_values, dtype=float),
        N=float(N),
        gamma=float(gamma),
        delta_ms=float(delta_ms),
        delta_app_ms=float(Delta_app_ms),
        g_type=gradient_base_for_axis(axis),
    )


def axis_from_gradient(
    gradient_values: np.ndarray,
    *,
    axis: str,
    N: float | None,
    gamma: float,
    delta_ms: float | None,
    Delta_app_ms: float | None,
) -> np.ndarray:
    base = normalize_axis_base(axis)
    gradient_values = np.asarray(gradient_values, dtype=float)
    if axis_uses_bvalue(base):
        return bvalue_from_gradient(
            gradient_values,
            axis=base,
            N=N,
            gamma=gamma,
            delta_ms=delta_ms,
            Delta_app_ms=Delta_app_ms,
        )
    return gradient_values


def build_axis_bundle(
    df: pd.DataFrame,
    *,
    axis: str,
    correction_factor: float = 1.0,
    gamma: float = 267.5221900,
    N: float | None = None,
    delta_ms: float | None = None,
    Delta_app_ms: float | None = None,
    side: int | None = None,
) -> AxisBundle:
    axis_base, axis_side = split_axis_side(axis)
    resolved_side = axis_side if side is None else int(side)
    if axis_side is not None and side is not None and int(axis_side) != int(side):
        raise ValueError(f"Axis {axis!r} and explicit side={side} disagree.")

    f_corr = coerce_positive_correction_factor(correction_factor)
    gradient_raw = extract_gradient_array(df, axis=axis_base, side=resolved_side)
    gradient_corr = gradient_raw * float(f_corr)

    bvalue_raw = _extract_stored_bvalue_array(df, axis=axis_base, side=resolved_side)
    bvalue_corr: np.ndarray | None = None

    if axis_uses_bvalue(axis_base):
        if N is None or delta_ms is None or Delta_app_ms is None:
            if bvalue_raw is None:
                raise ValueError(
                    f"Axis {axis!r} needs N, delta_ms, and Delta_app_ms to derive a corrected b-value."
                )
            if not np.isclose(float(f_corr), 1.0):
                raise ValueError(
                    f"Axis {axis!r} requested corrected b-values, but sequence timing parameters are missing."
                )
            bvalue_corr = np.asarray(bvalue_raw, dtype=float)
        else:
            if bvalue_raw is None:
                bvalue_raw = bvalue_from_gradient(
                    gradient_raw,
                    axis=axis_base,
                    N=N,
                    gamma=gamma,
                    delta_ms=delta_ms,
                    Delta_app_ms=Delta_app_ms,
                )
            bvalue_corr = bvalue_from_gradient(
                gradient_corr,
                axis=axis_base,
                N=N,
                gamma=gamma,
                delta_ms=delta_ms,
                Delta_app_ms=Delta_app_ms,
            )

    axis_raw = np.asarray(bvalue_raw, dtype=float) if axis_uses_bvalue(axis_base) else np.asarray(gradient_raw, dtype=float)
    axis_corr = np.asarray(bvalue_corr, dtype=float) if axis_uses_bvalue(axis_base) else np.asarray(gradient_corr, dtype=float)

    return AxisBundle(
        requested_axis=str(axis),
        axis_base=axis_base,
        gradient_base=gradient_base_for_axis(axis_base),
        side=resolved_side,
        gradient_raw=np.asarray(gradient_raw, dtype=float),
        gradient_corr=np.asarray(gradient_corr, dtype=float),
        bvalue_raw=None if bvalue_raw is None else np.asarray(bvalue_raw, dtype=float),
        bvalue_corr=None if bvalue_corr is None else np.asarray(bvalue_corr, dtype=float),
        axis_raw=axis_raw,
        axis_corr=axis_corr,
    )


__all__ = [
    "AxisBundle",
    "VALID_AXIS_BASES",
    "VALID_G_TYPES",
    "axes_share_gradient_family",
    "axis_from_gradient",
    "axis_uses_bvalue",
    "b_from_g",
    "build_axis_bundle",
    "bvalue_column_name",
    "bvalue_from_gradient",
    "coerce_positive_correction_factor",
    "default_plot_axis_for_fit",
    "extract_gradient_array",
    "gradient_base_for_axis",
    "gradient_column_name",
    "normalize_axis_base",
    "split_axis_side",
]
