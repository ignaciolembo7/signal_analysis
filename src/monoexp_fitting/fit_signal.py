from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from data_processing.io import sanitize_output_token, write_table_outputs
from fitting.b_from_g import b_from_g, gradient_base_for_axis, normalize_axis_base
from fitting.core import CurveFitParameter, chi2, fit_curve_fit_parameters, r2_score, rmse, rmse_log
from plotting.core import compact_float, render_xy_plot


AUTO_FIT_MIN_POINTS = 3
AUTO_FIT_MAX_POINTS = 9
AUTO_FIT_REL_TOL = 0.05
AUTO_FIT_ERR_FLOOR = 5e-3
AUTO_FIT_ABS_TOL = 1e-6

BVALUE_AXIS_COLUMNS = {
    "bvalue": "bvalue",
    "bvalue_g": "bvalue_g",
    "bvalue_g_lin_max": "bvalue_g_lin_max",
    "bvalue_thorsten": "bvalue_thorsten",
    "g": "bvalue_g",
    "g_lin_max": "bvalue_g_lin_max",
    "g_thorsten": "bvalue_thorsten",
}


@dataclass(frozen=True)
class FitOutputs:
    fit_params: pd.DataFrame
    fit_points: pd.DataFrame


def monoexp(b: np.ndarray, M0: float, D0_mm2_s: float) -> np.ndarray:
    return float(M0) * np.exp(-np.asarray(b, dtype=float) * float(D0_mm2_s))


def _unique_float(df: pd.DataFrame, col: str) -> float | None:
    if col not in df.columns:
        return None
    vals = pd.to_numeric(df[col], errors="coerce").dropna().unique()
    if len(vals) != 1:
        return None
    return float(vals[0])


def _unique_str(df: pd.DataFrame, col: str) -> str | None:
    if col not in df.columns:
        return None
    vals = pd.Series(df[col]).dropna().astype(str).str.strip().unique()
    vals = [v for v in vals if v]
    if len(vals) != 1:
        return None
    return str(vals[0])


def _coalesce_float(df: pd.DataFrame, explicit: float | None, cols: Sequence[str]) -> float | None:
    if explicit is not None:
        return float(explicit)
    for col in cols:
        value = _unique_float(df, col)
        if value is not None:
            return value
    return None


def resolve_b_axis(
    df: pd.DataFrame,
    *,
    axis: str,
    gamma: float = 267.5221900,
    N: float | None = None,
    delta_ms: float | None = None,
    Delta_app_ms: float | None = None,
    correction_factor: float = 1.0,
) -> tuple[np.ndarray, str]:
    """Return the monoexponential b axis in s/mm2."""
    base = normalize_axis_base(axis)
    b_col = BVALUE_AXIS_COLUMNS.get(base)
    if b_col is not None and b_col in df.columns:
        b = pd.to_numeric(df[b_col], errors="coerce").to_numpy(dtype=float)
    else:
        g_col = gradient_base_for_axis(base)
        if g_col not in df.columns:
            raise KeyError(f"Missing b-value/gradient column for axis {axis!r}.")
        n_val = _coalesce_float(df, N, ["N"])
        delta_val = _coalesce_float(df, delta_ms, ["delta_ms"])
        delta_app_val = _coalesce_float(df, Delta_app_ms, ["Delta_app_ms", "delta_app_ms"])
        if n_val is None or delta_val is None or delta_app_val is None:
            raise ValueError(
                f"Cannot derive b-values for axis {axis!r}; N, delta_ms, and Delta_app_ms are required."
            )
        g = pd.to_numeric(df[g_col], errors="coerce").to_numpy(dtype=float)
        b = b_from_g(
            g,
            N=float(n_val),
            gamma=float(gamma),
            delta_ms=float(delta_val),
            delta_app_ms=float(delta_app_val),
            g_type=g_col,
        )
        b_col = f"derived_from_{g_col}"

    f_corr = float(correction_factor)
    if np.isfinite(f_corr) and f_corr > 0.0 and not np.isclose(f_corr, 1.0):
        b = b * (f_corr**2.0)
    return np.asarray(b, dtype=float), str(b_col)


def fit_monoexp_prefix(
    b: np.ndarray,
    y: np.ndarray,
    *,
    fit_points: int,
    free_M0: bool = False,
    fix_M0: float = 1.0,
    D0_init: float = 0.0023,
    min_valid_points: int = AUTO_FIT_MIN_POINTS,
) -> dict:
    b = np.asarray(b, dtype=float)
    y = np.asarray(y, dtype=float)
    k = min(int(fit_points), len(b))
    prefix_mask = np.zeros(len(b), dtype=bool)
    if k > 0:
        prefix_mask[:k] = True
    valid_mask = np.isfinite(b) & np.isfinite(y) & (y > 0.0) & (b >= 0.0)
    fit_mask = prefix_mask & valid_mask
    b_fit = b[fit_mask].copy()
    y_fit = y[fit_mask].copy()

    result = {
        "ok": False,
        "fit_points": int(k),
        "n_fit": int(len(b_fit)),
        "fit_mask": fit_mask,
        "msg": "",
    }
    if k <= 0:
        result["msg"] = "fit_points must be > 0."
        return result
    if len(b_fit) < int(min_valid_points):
        result["msg"] = "Too few valid points."
        return result

    D0_seed = float(D0_init) if np.isfinite(D0_init) and float(D0_init) > 0.0 else 0.0023
    try:
        method_label = "curve_fit(M0,D0)" if free_M0 else "curve_fit(D0) M0_fixed"

        def model_from_params(params: dict[str, float]) -> np.ndarray:
            return monoexp(b_fit, float(params["M0"]), float(params["D0_mm2_s"]))

        fit = fit_curve_fit_parameters(
            model_from_params,
            y_fit,
            parameters=[
                CurveFitParameter("M0", 1.0 if free_M0 else float(fix_M0), 0.0, 100.0, bool(free_M0)),
                CurveFitParameter("D0_mm2_s", D0_seed, D0_seed / 10.0, 2.0 * D0_seed, True),
            ],
            x=b_fit,
            maxfev=40000,
            method=method_label,
        )
    except Exception as exc:
        result["msg"] = f"curve_fit failed: {exc}"
        return result

    M0_hat = float(fit.values["M0"])
    D0_hat = float(fit.values["D0_mm2_s"])
    yhat = monoexp(b_fit, M0_hat, D0_hat)
    D0_err = float(fit.errors.get("D0_err_mm2_s", np.nan))
    M0_err = float(fit.errors.get("M0_err", np.nan)) if free_M0 else np.nan

    result.update(
        ok=True,
        msg="",
        M0=M0_hat,
        M0_err=M0_err if np.isfinite(M0_err) else np.nan,
        D0_mm2_s=D0_hat,
        D0_err_mm2_s=D0_err if np.isfinite(D0_err) else np.nan,
        D0_m2_ms=D0_hat * 1e-9,
        D0_err_m2_ms=D0_err * 1e-9 if np.isfinite(D0_err) else np.nan,
        rmse=rmse(y_fit, yhat),
        chi2=chi2(y_fit, yhat),
        r2=r2_score(y_fit, yhat),
        rmse_log=rmse_log(y_fit, yhat),
        method=str(fit.method),
    )
    return result


def select_monoexp_fit_result(
    b: np.ndarray,
    y: np.ndarray,
    *,
    fit_points: int | None = None,
    auto_fit_points: bool = False,
    free_M0: bool = False,
    fix_M0: float = 1.0,
    D0_init: float = 0.0023,
    auto_fit_min_points: int = AUTO_FIT_MIN_POINTS,
    auto_fit_max_points: int | None = AUTO_FIT_MAX_POINTS,
    auto_fit_rel_tol: float = AUTO_FIT_REL_TOL,
    auto_fit_err_floor: float = AUTO_FIT_ERR_FLOOR,
) -> dict:
    b = np.asarray(b, dtype=float)
    y = np.asarray(y, dtype=float)

    if auto_fit_points:
        k_min = max(1, int(auto_fit_min_points))
        k_max = len(b) if auto_fit_max_points is None else min(int(auto_fit_max_points), len(b))
        if k_max < k_min:
            return {
                "ok": False,
                "fit_points": np.nan,
                "n_fit": 0,
                "fit_strategy": "auto",
                "auto_fit_metric": "rmse_log",
                "auto_fit_score": np.nan,
                "fit_mask": np.zeros(len(b), dtype=bool),
                "msg": f"Invalid auto_fit_points range: min_k={k_min}, max_k={k_max}.",
            }

        fallback: dict | None = None
        accepted: list[dict] = []
        rejected: list[tuple[int, float]] = []
        tested_until = k_min - 1
        score_limit = max(float(auto_fit_rel_tol), float(auto_fit_err_floor)) + AUTO_FIT_ABS_TOL
        for k in range(k_min, k_max + 1):
            tested_until = k
            cand = fit_monoexp_prefix(
                b,
                y,
                fit_points=k,
                free_M0=free_M0,
                fix_M0=fix_M0,
                D0_init=D0_init,
                min_valid_points=k_min,
            )
            if not cand["ok"]:
                continue
            if fallback is None:
                fallback = cand
            score = float(cand["rmse_log"])
            if np.isfinite(score) and score <= score_limit:
                accepted.append(cand)
            else:
                rejected.append((k, score))

        chosen = accepted[-1] if accepted else fallback
        if chosen is not None:
            selected = dict(chosen)
            selected["fit_strategy"] = "auto"
            selected["auto_fit_metric"] = "rmse_log"
            selected["auto_fit_score"] = float(selected["rmse_log"])
            rejection_summary = ""
            if rejected:
                first_k, first_score = rejected[0]
                rejection_summary = (
                    f" First rejected prefix: k={first_k}, rmse_log={first_score:.4g}."
                )
            fallback_summary = ""
            if not accepted:
                fallback_summary = " No prefix met the threshold; retained the smallest valid prefix."
            selected["msg"] = (
                f"Auto fit_points selected {int(selected['fit_points'])} "
                f"after testing every prefix k={k_min}..{tested_until} with "
                f"rmse_log <= {score_limit:.4g}.{rejection_summary}{fallback_summary}"
            )
            selected["method"] = (
                f"{selected['method']} | "
                f"auto_fit_points_largest_prefix(rmse_log_limit={score_limit:.4g}, "
                f"min_k={k_min}, max_k={k_max})"
            )
            return selected

        valid_total = int(np.sum(np.isfinite(b) & np.isfinite(y) & (y > 0.0) & (b >= 0.0)))
        return {
            "ok": False,
            "fit_points": np.nan,
            "n_fit": valid_total,
            "fit_strategy": "auto",
            "auto_fit_metric": "rmse_log",
            "auto_fit_score": np.nan,
            "fit_mask": np.zeros(len(b), dtype=bool),
            "msg": f"No valid candidate was found for auto_fit_points in the range k={k_min}..{k_max}.",
        }

    selected_fit_points = int(fit_points if fit_points is not None else 6)
    selected = fit_monoexp_prefix(
        b,
        y,
        fit_points=selected_fit_points,
        free_M0=free_M0,
        fix_M0=fix_M0,
        D0_init=D0_init,
        min_valid_points=AUTO_FIT_MIN_POINTS,
    )
    selected["fit_strategy"] = "fixed"
    selected["auto_fit_metric"] = np.nan
    selected["auto_fit_score"] = np.nan
    return selected


def _fit_label(row: Mapping[str, object]) -> str:
    return (
        f"monoexp, M0={compact_float(row.get('M0'))}, "
        f"D0={compact_float(row.get('D0_mm2_s'))} mm2/s"
    )


def plot_monoexp_fit(
    *,
    b: np.ndarray,
    y: np.ndarray,
    fit_row: Mapping[str, object],
    fit_mask: np.ndarray,
    out_png: Path,
    ycol: str,
    b_label: str,
) -> None:
    b = np.asarray(b, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(b) & np.isfinite(y)
    b_plot = b[valid]
    y_plot = y[valid]
    order = np.argsort(b_plot)
    b_plot = b_plot[order]
    y_plot = y_plot[order]

    b_dense = np.linspace(0.0, float(np.nanmax(b_plot)) if b_plot.size else 1.0, 250)
    fit_y = None
    if bool(fit_row.get("ok", False)):
        fit_y = monoexp(b_dense, float(fit_row["M0"]), float(fit_row["D0_mm2_s"]))

    text_lines = [
        f"M0={compact_float(fit_row.get('M0'))}",
        f"D0={compact_float(fit_row.get('D0_mm2_s'))} mm2/s",
        f"fit_points={compact_float(fit_row.get('fit_points'), digits=0)}",
        f"rmse_log={compact_float(fit_row.get('rmse_log'), digits=4)}",
    ]
    title = (
        f"Monoexp signal fit | roi={fit_row.get('roi', 'NA')} | "
        f"dir={fit_row.get('direction', 'NA')} | td={compact_float(fit_row.get('td_ms'))} | "
        f"N={compact_float(fit_row.get('N'), digits=0)}"
    )
    render_xy_plot(
        x=b_plot,
        y=y_plot,
        out_png=out_png,
        title=title,
        xlabel=f"{b_label} (s/mm2)",
        ylabel=ycol,
        data_label="data",
        connect_data=False,
        fit_x=b_dense if fit_y is not None else None,
        fit_y=fit_y,
        fit_label=_fit_label(fit_row),
        highlight_x=b[fit_mask],
        highlight_y=y[fit_mask],
        highlight_label="used for fit",
        text_lines=text_lines,
        text_position="inside_upper_right",
        legend_position="lower_left",
        yscale="log",
    )


def _metadata_from_group(group: pd.DataFrame, columns: Sequence[str]) -> dict[str, object]:
    out: dict[str, object] = {}
    for col in columns:
        if col not in group.columns:
            out[col] = np.nan
            continue
        vals = pd.Series(group[col]).dropna().unique()
        out[col] = vals[0] if len(vals) == 1 else np.nan
    return out


def _correction_factor_from_group(group: pd.DataFrame, col: str) -> tuple[float, str]:
    if col not in group.columns:
        return np.nan, f"Missing required gradient-correction column {col!r}."
    values = pd.to_numeric(group[col], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(values) & (values > 0.0)
    if not np.all(valid):
        return np.nan, f"Missing {col} for this curve."
    return float(np.mean(values)), ""


def fit_signal_monoexp(
    df: pd.DataFrame,
    *,
    directions: Sequence[str] | None = None,
    rois: Sequence[str] | None = None,
    ycol: str = "value_norm",
    b_axis: str = "bvalue_g",
    fit_points: int | None = None,
    auto_fit_points: bool = True,
    auto_fit_min_points: int = AUTO_FIT_MIN_POINTS,
    auto_fit_max_points: int | None = AUTO_FIT_MAX_POINTS,
    auto_fit_rel_tol: float = AUTO_FIT_REL_TOL,
    auto_fit_err_floor: float = AUTO_FIT_ERR_FLOOR,
    free_M0: bool = False,
    fix_M0: float = 1.0,
    D0_init: float = 0.0023,
    gamma: float = 267.5221900,
    stat_keep: str = "avg",
    f_by_direction: Mapping[str, float] | None = None,
    correction_factor_col: str | None = None,
    outdir_plots: Path | None = None,
    plot_prefix: str = "",
) -> FitOutputs:
    dfa = df.copy()
    if "stat" in dfa.columns and stat_keep and str(stat_keep).upper() != "ALL":
        dfa = dfa[dfa["stat"].astype(str) == str(stat_keep)].copy()
    if directions is not None:
        wanted = {str(v) for v in directions}
        dfa = dfa[dfa["direction"].astype(str).isin(wanted)].copy()
    if rois is not None:
        wanted = {str(v) for v in rois}
        dfa = dfa[dfa["roi"].astype(str).isin(wanted)].copy()
    if dfa.empty:
        return FitOutputs(fit_params=pd.DataFrame(), fit_points=pd.DataFrame())
    if ycol not in dfa.columns:
        raise KeyError(f"Missing y column {ycol!r}.")
    if "b_step" not in dfa.columns:
        raise KeyError("Missing required column 'b_step'.")

    if outdir_plots is not None:
        outdir_plots.mkdir(parents=True, exist_ok=True)

    params_rows: list[dict] = []
    point_rows: list[dict] = []
    group_cols = [c for c in ["subj", "sheet", "td_ms", "N", "Hz", "roi", "direction"] if c in dfa.columns]
    meta_cols = [
        "subj", "sheet", "type", "sequence", "group", "protocol", "source_file",
        "td_ms", "N", "Hz", "TE", "TR", "max_dur_ms", "tm_ms", "delta_ms", "Delta_app_ms", "bmax",
    ]

    for key, group in dfa.groupby(group_cols, sort=False, dropna=False):
        del key
        g = group.sort_values("b_step", kind="stable").copy()
        direction = str(g["direction"].iloc[0]) if "direction" in g.columns else "NA"
        roi = str(g["roi"].iloc[0]) if "roi" in g.columns else "NA"
        f_corr = float(f_by_direction.get(direction, 1.0)) if f_by_direction else 1.0
        correction_msg = ""
        if correction_factor_col:
            f_corr, correction_msg = _correction_factor_from_group(g, correction_factor_col)
        y = pd.to_numeric(g[ycol], errors="coerce").to_numpy(dtype=float)
        b, b_col = resolve_b_axis(
            g,
            axis=b_axis,
            gamma=gamma,
            correction_factor=f_corr if np.isfinite(f_corr) else 1.0,
        )
        if correction_msg:
            result = {
                "ok": False,
                "fit_points": np.nan,
                "n_fit": 0,
                "fit_mask": np.zeros(len(g), dtype=bool),
                "fit_strategy": "auto" if auto_fit_points else "fixed",
                "auto_fit_metric": "rmse_log" if auto_fit_points else np.nan,
                "auto_fit_score": np.nan,
                "msg": correction_msg,
            }
        else:
            result = select_monoexp_fit_result(
                b,
                y,
                fit_points=fit_points,
                auto_fit_points=auto_fit_points,
                free_M0=free_M0,
                fix_M0=fix_M0,
                D0_init=D0_init,
                auto_fit_min_points=auto_fit_min_points,
                auto_fit_max_points=auto_fit_max_points,
                auto_fit_rel_tol=auto_fit_rel_tol,
                auto_fit_err_floor=auto_fit_err_floor,
            )
        fit_mask = np.asarray(result.get("fit_mask", np.zeros(len(g), dtype=bool)), dtype=bool)
        meta = _metadata_from_group(g, meta_cols)
        fit_points_value = result.get("fit_points", np.nan)
        fit_points_out = int(fit_points_value) if np.isfinite(fit_points_value) else np.nan
        row = {
            **meta,
            "roi": roi,
            "direction": direction,
            "stat": str(stat_keep),
            "fit_kind": "monoexp",
            "model": "monoexp",
            "ycol": str(ycol),
            "b_axis": str(b_axis),
            "b_col": str(b_col),
            "fit_points": fit_points_out,
            "fit_strategy": str(result.get("fit_strategy", "auto" if auto_fit_points else "fixed")),
            "auto_fit_metric": result.get("auto_fit_metric", np.nan),
            "auto_fit_score": result.get("auto_fit_score", np.nan),
            "f_corr": float(f_corr),
            "b_corr_scale": float(f_corr) ** 2.0,
            "n_points": int(len(g)),
            "n_fit": int(result.get("n_fit", 0)),
            "ok": bool(result.get("ok", False)),
            "msg": str(result.get("msg", "")),
        }
        for name in [
            "M0", "M0_err", "D0_mm2_s", "D0_err_mm2_s", "D0_m2_ms", "D0_err_m2_ms",
            "rmse", "rmse_log", "chi2", "r2", "method",
        ]:
            row[name] = result.get(name, np.nan)
        params_rows.append(row)

        for idx, (_, src_row) in enumerate(g.iterrows()):
            point_rows.append(
                {
                    **{col: meta.get(col, np.nan) for col in meta_cols},
                    "roi": roi,
                    "direction": direction,
                    "stat": str(stat_keep),
                    "fit_kind": "monoexp",
                    "model": "monoexp",
                    "ycol": str(ycol),
                    "b_axis": str(b_axis),
                    "b_col": str(b_col),
                    "fit_points": fit_points_out,
                    "fit_strategy": row["fit_strategy"],
                    "f_corr": float(f_corr),
                    "b_step": src_row.get("b_step", np.nan),
                    "bvalue_used": float(b[idx]) if np.isfinite(b[idx]) else np.nan,
                    "y": float(y[idx]) if np.isfinite(y[idx]) else np.nan,
                    "used_for_fit": bool(fit_mask[idx]) if idx < fit_mask.size else False,
                }
            )

        if outdir_plots is not None:
            pieces = [
                plot_prefix,
                sanitize_output_token(meta.get("sheet", "sheet")),
                f"N{sanitize_output_token(meta.get('N', 'NA'))}",
                f"td{sanitize_output_token(meta.get('td_ms', 'NA'))}",
                sanitize_output_token(roi),
                sanitize_output_token(direction),
            ]
            out_png = outdir_plots / (".".join([p for p in pieces if p]) + ".monoexp.png")
            plot_monoexp_fit(
                b=b,
                y=y,
                fit_row=row,
                fit_mask=fit_mask,
                out_png=out_png,
                ycol=ycol,
                b_label=str(b_axis),
            )

    return FitOutputs(fit_params=pd.DataFrame(params_rows), fit_points=pd.DataFrame(point_rows))


def write_fit_outputs(
    outputs: FitOutputs,
    out_root: Path,
    *,
    b_axis: str,
    ycol: str,
    model: str = "monoexp",
) -> tuple[Path | None, Path | None]:
    out_root.mkdir(parents=True, exist_ok=True)
    axis_tag = sanitize_output_token(b_axis)
    y_tag = sanitize_output_token(ycol)
    model_tag = sanitize_output_token(model)
    params_path = out_root / f"signal_fit_params.{model_tag}.{axis_tag}.{y_tag}.parquet"
    points_path = out_root / f"signal_fit_points.{model_tag}.{axis_tag}.{y_tag}.parquet"
    written_params: Path | None = None
    written_points: Path | None = None
    if not outputs.fit_params.empty:
        written_params = write_table_outputs(outputs.fit_params, params_path, xlsx_path=params_path.with_suffix(".xlsx"))
    if not outputs.fit_points.empty:
        written_points = write_table_outputs(outputs.fit_points, points_path, xlsx_path=points_path.with_suffix(".xlsx"))
    return written_params, written_points


__all__ = [
    "AUTO_FIT_ABS_TOL",
    "AUTO_FIT_ERR_FLOOR",
    "AUTO_FIT_MAX_POINTS",
    "AUTO_FIT_MIN_POINTS",
    "AUTO_FIT_REL_TOL",
    "FitOutputs",
    "fit_monoexp_prefix",
    "fit_signal_monoexp",
    "monoexp",
    "resolve_b_axis",
    "select_monoexp_fit_result",
    "write_fit_outputs",
]
