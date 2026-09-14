from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize_scalar

from data_processing.io import sanitize_output_token
from models.model_fitting import M_ogse_free
from monoexp_fitting.fit_signal import FitOutputs
from plotting.core import render_xy_plot


def _ogse_free_model(
    td_ms: float,
    gradient_mTm: np.ndarray,
    n_value: float,
    m0: float,
    d0_m2_ms: float,
) -> np.ndarray:
    """Evaluate the free OGSE signal using x = td / N."""
    x_ms = float(td_ms) / float(n_value)
    return M_ogse_free(
        float(td_ms),
        np.asarray(gradient_mTm, dtype=float),
        float(n_value),
        x_ms,
        float(m0),
        float(d0_m2_ms),
    )


def _fit_ogse_free_curve(
    *,
    td_ms: float,
    gradient_mTm: np.ndarray,
    n_value: float,
    signal: np.ndarray,
    free_m0: bool,
    fixed_m0: float,
    d0_initial_m2_ms: float,
) -> dict[str, object]:
    gradient_mTm = np.asarray(gradient_mTm, dtype=float)
    signal = np.asarray(signal, dtype=float)
    valid = np.isfinite(gradient_mTm) & np.isfinite(signal) & (signal > 0.0)
    gradient_fit = gradient_mTm[valid]
    signal_fit = signal[valid]
    n_fit = int(len(signal_fit))

    empty = {
        "ok": False,
        "n_fit": n_fit,
        "M0": np.nan,
        "D0_m2_ms": np.nan,
        "D0_mm2_s": np.nan,
        "rmse": np.nan,
    }
    if n_fit < 3:
        return {**empty, "msg": "Too few valid points."}

    seed = float(d0_initial_m2_ms)
    if not np.isfinite(seed) or seed <= 0.0:
        seed = 2.3e-12
    lower = max(seed / 100.0, 1e-15)
    upper = min(seed * 100.0, 1e-9)

    try:
        if free_m0:
            def model(_unused: np.ndarray, m0: float, d0: float) -> np.ndarray:
                return _ogse_free_model(td_ms, gradient_fit, n_value, m0, d0)

            optimum, _ = curve_fit(
                model,
                np.zeros_like(signal_fit),
                signal_fit,
                p0=[float(fixed_m0), seed],
                bounds=([0.0, lower], [5.0, upper]),
                maxfev=400000,
            )
            m0 = float(optimum[0])
            d0 = float(optimum[1])
        else:
            def loss(log_d0: float) -> float:
                candidate = float(np.exp(log_d0))
                prediction = _ogse_free_model(
                    td_ms, gradient_fit, n_value, float(fixed_m0), candidate
                )
                if prediction.shape != signal_fit.shape or not np.all(np.isfinite(prediction)):
                    return np.inf
                return float(np.sum((signal_fit - prediction) ** 2))

            grid = np.linspace(float(np.log(lower)), float(np.log(upper)), 96)
            losses = np.asarray([loss(value) for value in grid], dtype=float)
            best_index = int(np.nanargmin(losses))
            best_log = float(grid[best_index])
            refine_lower = float(grid[max(0, best_index - 1)])
            refine_upper = float(grid[min(len(grid) - 1, best_index + 1)])
            if refine_upper > refine_lower:
                optimum = minimize_scalar(
                    loss,
                    bounds=(refine_lower, refine_upper),
                    method="bounded",
                    options={"xatol": 1e-8},
                )
                if optimum.success and np.isfinite(optimum.fun) and optimum.fun <= losses[best_index]:
                    best_log = float(optimum.x)
            d0 = float(np.exp(best_log))
            m0 = float(fixed_m0)

        prediction = _ogse_free_model(td_ms, gradient_fit, n_value, m0, d0)
        fit_rmse = float(np.sqrt(np.mean((signal_fit - prediction) ** 2)))
        return {
            "ok": True,
            "n_fit": n_fit,
            "M0": m0,
            "D0_m2_ms": d0,
            "D0_mm2_s": d0 * 1e9,
            "rmse": fit_rmse,
            "msg": "",
        }
    except Exception as exc:
        return {**empty, "msg": str(exc)}


def _group_correction_factor(group: pd.DataFrame, column: str | None) -> tuple[float, str]:
    if column is None:
        return 1.0, ""
    if column not in group.columns:
        return np.nan, f"Missing required gradient-correction column {column!r}."
    values = pd.to_numeric(group[column], errors="coerce").to_numpy(float)
    if values.size == 0 or not np.all(np.isfinite(values) & (values > 0.0)):
        return np.nan, f"Missing {column} for this curve."
    return float(np.mean(values)), ""


def _unique_float(group: pd.DataFrame, column: str) -> float | None:
    if column not in group.columns:
        return None
    values = pd.to_numeric(group[column], errors="coerce").dropna().unique()
    return float(values[0]) if len(values) == 1 else None


def _unique_value(group: pd.DataFrame, column: str) -> object:
    if column not in group.columns:
        return np.nan
    values = pd.Series(group[column]).dropna().unique()
    return values[0] if len(values) == 1 else np.nan


def fit_ogse_free_signal(
    df: pd.DataFrame,
    *,
    directions: Sequence[str] | None = None,
    rois: Sequence[str] | None = None,
    ycol: str = "value_norm",
    g_axis: str = "g",
    free_M0: bool = False,
    fix_M0: float = 1.0,
    D0_init_mm2_s: float = 0.0023,
    stat_keep: str = "avg",
    correction_factor_col: str | None = None,
    outdir_plots: Path | None = None,
) -> FitOutputs:
    work = df.copy()
    if "stat" in work.columns and stat_keep and str(stat_keep).upper() != "ALL":
        work = work[work["stat"].astype(str) == str(stat_keep)].copy()
    if directions is not None:
        wanted = {str(value) for value in directions}
        work = work[work["direction"].astype(str).isin(wanted)].copy()
    if rois is not None:
        wanted = {str(value) for value in rois}
        work = work[work["roi"].astype(str).isin(wanted)].copy()
    if work.empty:
        return FitOutputs(fit_params=pd.DataFrame(), fit_points=pd.DataFrame())
    for required in ("b_step", "direction", "roi", "td_ms", "N", g_axis, ycol):
        if required not in work.columns:
            raise KeyError(f"Missing required column {required!r}.")

    if outdir_plots is not None:
        outdir_plots.mkdir(parents=True, exist_ok=True)

    group_cols = [
        column
        for column in ("subj", "sheet", "td_ms", "N", "Hz", "roi", "direction")
        if column in work.columns
    ]
    metadata_cols = [
        "subj", "sheet", "type", "sequence", "group", "protocol", "source_file",
        "td_ms", "N", "Hz", "TE", "TR", "max_dur_ms", "tm_ms", "delta_ms",
        "Delta_app_ms", "bmax",
    ]
    parameter_rows: list[dict[str, object]] = []
    point_rows: list[dict[str, object]] = []

    for _, group in work.groupby(group_cols, sort=False, dropna=False):
        curve = group.sort_values("b_step", kind="stable").copy()
        td_ms = _unique_float(curve, "td_ms")
        n_value = _unique_float(curve, "N")
        if td_ms is None or n_value is None:
            raise ValueError("Each OGSE free-model curve must have one td_ms and one N value.")

        direction = str(curve["direction"].iloc[0])
        roi = str(curve["roi"].iloc[0])
        factor, correction_message = _group_correction_factor(curve, correction_factor_col)
        gradient_raw = pd.to_numeric(curve[g_axis], errors="coerce").to_numpy(float)
        gradient_used = gradient_raw * factor if np.isfinite(factor) else gradient_raw
        signal = pd.to_numeric(curve[ycol], errors="coerce").to_numpy(float)

        if correction_message:
            result = {
                "ok": False, "n_fit": 0, "M0": np.nan, "D0_m2_ms": np.nan,
                "D0_mm2_s": np.nan, "rmse": np.nan, "msg": correction_message,
            }
        else:
            result = _fit_ogse_free_curve(
                td_ms=td_ms,
                gradient_mTm=gradient_used,
                n_value=n_value,
                signal=signal,
                free_m0=bool(free_M0),
                fixed_m0=float(fix_M0),
                d0_initial_m2_ms=float(D0_init_mm2_s) * 1e-9,
            )

        valid = np.isfinite(gradient_used) & np.isfinite(signal) & (signal > 0.0)
        metadata = {column: _unique_value(curve, column) for column in metadata_cols}
        row = {
            **metadata,
            "roi": roi,
            "direction": direction,
            "stat": str(stat_keep),
            "fit_kind": "ogse_free_signal",
            "model": "ogse_free",
            "ycol": str(ycol),
            "b_axis": str(g_axis),
            "g_col": str(g_axis),
            "x_model_ms": td_ms / n_value,
            "f_corr": float(factor),
            "g_corr_scale": float(factor),
            "b_corr_scale": float(factor) ** 2.0,
            "n_points": int(len(curve)),
            "n_fit": int(result.get("n_fit", 0)),
            "ok": bool(result.get("ok", False)),
            "msg": str(result.get("msg", "")),
            "M0": result.get("M0", np.nan),
            "D0_m2_ms": result.get("D0_m2_ms", np.nan),
            "D0_mm2_s": result.get("D0_mm2_s", np.nan),
            "rmse": result.get("rmse", np.nan),
            "method": "curve_fit" if free_M0 else "bounded_logD",
        }
        parameter_rows.append(row)

        for index, (_, source_row) in enumerate(curve.iterrows()):
            point_rows.append(
                {
                    **metadata,
                    "roi": roi,
                    "direction": direction,
                    "stat": str(stat_keep),
                    "fit_kind": "ogse_free_signal",
                    "model": "ogse_free",
                    "ycol": str(ycol),
                    "b_axis": str(g_axis),
                    "x_model_ms": td_ms / n_value,
                    "f_corr": float(factor),
                    "b_step": source_row.get("b_step", np.nan),
                    "g_raw": float(gradient_raw[index]) if np.isfinite(gradient_raw[index]) else np.nan,
                    "g_used": float(gradient_used[index]) if np.isfinite(gradient_used[index]) else np.nan,
                    "y": float(signal[index]) if np.isfinite(signal[index]) else np.nan,
                    "used_for_fit": bool(valid[index]),
                }
            )

        if outdir_plots is not None and bool(result.get("ok", False)) and valid.any():
            g_fine = np.linspace(0.0, float(np.nanmax(gradient_used[valid])), 250)
            fitted = _ogse_free_model(
                td_ms, g_fine, n_value, float(result["M0"]), float(result["D0_m2_ms"])
            )
            tokens = [
                sanitize_output_token(metadata.get("sheet", "sheet")),
                f"N{sanitize_output_token(n_value)}",
                f"td{sanitize_output_token(td_ms)}",
                sanitize_output_token(roi),
                sanitize_output_token(direction),
            ]
            out_png = outdir_plots / (".".join(tokens) + ".ogse_free.png")
            render_xy_plot(
                x=gradient_used[valid],
                y=signal[valid],
                out_png=out_png,
                title=(
                    f"OGSE free signal fit | ROI={roi} | direction={direction} | "
                    f"td_ms={td_ms:g} | N={n_value:g} | f_corr={factor:.4g}"
                ),
                xlabel=f"{g_axis} [mT/m]",
                ylabel=str(ycol),
                connect_data=False,
                fit_x=g_fine,
                fit_y=np.asarray(fitted, dtype=float),
                fit_label=(
                    f"ogse_free, M0={float(result['M0']):.4g}, "
                    f"D0={float(result['D0_mm2_s']):.6g} mm2/s"
                ),
            )

    return FitOutputs(
        fit_params=pd.DataFrame(parameter_rows),
        fit_points=pd.DataFrame(point_rows),
    )


__all__ = ["fit_ogse_free_signal"]
