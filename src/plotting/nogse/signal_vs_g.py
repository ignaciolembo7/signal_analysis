from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from data_processing.result_signals import X_SWEEP_NAME
from plotting.core import (
    XYSeries,
    compact_float,
    ensure_dir,
    render_multi_series_plot,
    render_xy_plot,
    sanitize_token,
)


SIGNAL_KEYS = ["subj", "roi", "N", "td_ms", "direction", "b_step"]


def _signal_keys(df: pd.DataFrame) -> list[str]:
    """Return curve identity columns, keeping acquisition types independent."""
    keys = list(SIGNAL_KEYS)
    if "session_id" in df.columns and df["session_id"].notna().any():
        keys.insert(1, "session_id")
    if "type" in df.columns:
        keys.append("type")
    return keys


def x_sweep_mask(df: pd.DataFrame) -> pd.Series:
    """Rows of NOGSE x sweeps (autox, sequence version v2.00.3): G fixed, x changes from volume to volume.

    Same rule as ``data_processing.result_signals.is_x_sweep_params``: type ``AUTOX`` or a protocol
    name with ``<D>dir<X>x``.
    """
    mask = pd.Series(False, index=df.index)
    if "type" in df.columns:
        mask |= df["type"].astype(str).str.strip().str.upper().eq("AUTOX")
    if "protocol" in df.columns:
        mask |= df["protocol"].astype(str).str.contains(X_SWEEP_NAME)
    return mask


def analysis_id_from_path(path: Path) -> str:
    stem = path.stem
    if stem.endswith(".long"):
        stem = stem[: -len(".long")]
    return stem


def split_all_or_values(values: Sequence[str] | None) -> list[str] | None:
    if values is None:
        return None
    out = [str(v) for v in values]
    if len(out) == 1 and out[0].upper() == "ALL":
        return None
    return out


def unique_scalar(df: pd.DataFrame, col: str):
    if col not in df.columns:
        return None
    values = pd.Series(df[col]).dropna().unique().tolist()
    if len(values) != 1:
        return None
    return values[0]


def prepare_signal_series(
    avg_df: pd.DataFrame,
    std_df: pd.DataFrame | None,
    *,
    xcol: str,
    ycol: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    data = avg_df.copy()
    data[xcol] = pd.to_numeric(data[xcol], errors="coerce")
    data[ycol] = pd.to_numeric(data[ycol], errors="coerce")
    data = data.dropna(subset=[xcol, ycol]).sort_values(xcol)
    if data.empty:
        return np.array([]), np.array([]), None

    sigma = None
    if std_df is not None and not std_df.empty:
        err = std_df.copy()
        err[xcol] = pd.to_numeric(err[xcol], errors="coerce")
        err["value"] = pd.to_numeric(err["value"], errors="coerce")
        err = err.dropna(subset=[xcol, "value"]).sort_values(xcol)
        if not err.empty:
            sigma = err["value"].to_numpy(dtype=float)
            if ycol == "value_norm" and "S0" in data.columns:
                s0 = pd.to_numeric(data["S0"], errors="coerce").to_numpy(dtype=float)
                with np.errstate(divide="ignore", invalid="ignore"):
                    sigma = sigma / s0
            if sigma.shape[0] != data.shape[0]:
                sigma = None

    return data[xcol].to_numpy(dtype=float), data[ycol].to_numpy(dtype=float), sigma


def _prepare_avg_std(
    df: pd.DataFrame,
    *,
    xcol: str,
    ycol: str,
    stat: str,
    extra_keys: Sequence[str] = (),
) -> pd.DataFrame:
    avg = df[df["stat"].astype(str) == str(stat)].copy()
    std = df[df["stat"].astype(str) == "std"].copy()

    if avg.empty:
        raise ValueError(f"No rows found for stat={stat!r}.")

    signal_keys = _signal_keys(df) + list(extra_keys)
    avg_cols = signal_keys + [xcol, ycol]
    if ycol == "value_norm" and "S0" in avg.columns:
        avg_cols.append("S0")
    avg = avg[avg_cols].rename(columns={xcol: "x_raw", ycol: "y_mean"})

    if std.empty:
        avg["y_std"] = np.nan
        merged = avg
    else:
        std = std[signal_keys + ["value"]].rename(columns={"value": "y_std"})
        merged = avg.merge(std, on=signal_keys, how="left")

    merged["y_mean"] = pd.to_numeric(merged["y_mean"], errors="coerce")
    merged["y_std"] = pd.to_numeric(merged["y_std"], errors="coerce")
    if ycol == "value_norm" and "S0" in merged.columns:
        s0 = pd.to_numeric(merged["S0"], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            merged["y_std"] = merged["y_std"] / s0
    return merged.sort_values(["subj", "roi", "N", "direction", "td_ms", "b_step"], kind="stable")


def _finite_sigma_or_none(values: pd.Series) -> np.ndarray | None:
    sigma = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return sigma if np.isfinite(sigma).any() else None


def _format_value(value: object, *, digits: int = 3) -> str:
    try:
        return compact_float(float(value), digits=digits)
    except (TypeError, ValueError):
        return str(value)


def _sort_key(value: object) -> tuple[int, float | str]:
    try:
        return (0, float(value))
    except (TypeError, ValueError):
        return (1, str(value))


def _plot_token(parts: dict[str, object]) -> str:
    return "_".join(f"{sanitize_token(key)}-{sanitize_token(_format_value(value))}" for key, value in parts.items())


TITLE_KEYS = {"session_id": "ses"}


def _plot_title(parts: dict[str, object], *, stat: str) -> str:
    title_parts = [f"{TITLE_KEYS.get(key, key)}={_format_value(value)}" for key, value in parts.items()]
    title_parts.append(f"stat={stat}")
    return " | ".join(title_parts)


def _series_label(column: str, value: object) -> str:
    digits = 0 if column == "N" else 3
    return f"{column}={_format_value(value, digits=digits)}"


def _make_series(group: pd.DataFrame, *, label: str) -> XYSeries | None:
    curve = group.dropna(subset=["x", "y_mean"]).sort_values(["x", "b_step"], kind="stable")
    if curve.empty:
        return None
    return XYSeries(
        x=curve["x"].to_numpy(float),
        y=curve["y_mean"].to_numpy(float),
        sigma=_finite_sigma_or_none(curve["y_std"]),
        label=label,
    )


def _plot_series_by(
    *,
    merged: pd.DataFrame,
    out_dir: Path,
    group_cols: list[str],
    series_col: str,
    file_prefix: str,
    xlabel: str,
    ylabel: str,
    stat: str,
) -> list[Path]:
    out_paths: list[Path] = []
    for keys, group in merged.groupby(group_cols, sort=False):
        key_tuple = keys if isinstance(keys, tuple) else (keys,)
        group_values = dict(zip(group_cols, key_tuple))
        series: list[XYSeries] = []

        for value in sorted(group[series_col].drop_duplicates().tolist(), key=_sort_key):
            curve = _make_series(group[group[series_col].eq(value)], label=_series_label(series_col, value))
            if curve is not None:
                series.append(curve)

        if not series:
            continue

        out_png = out_dir / f"{file_prefix}_{_plot_token(group_values)}.png"
        render_multi_series_plot(
            series=series,
            out_png=out_png,
            title=_plot_title(group_values, stat=stat),
            xlabel=xlabel,
            ylabel=ylabel,
            legend_title=series_col,
        )
        out_paths.append(out_png)

    return out_paths


def fit_parameter_fragments(fit_row: dict[str, object]) -> list[str]:
    fragments: list[str] = []
    if "tc_ms" in fit_row:
        fragments.append(f"tc_ms={compact_float(fit_row.get('tc_ms'))}")
    if "alpha" in fit_row:
        fragments.append(f"alpha={compact_float(fit_row.get('alpha'))}")
    if "M0" in fit_row:
        fragments.append(f"M0={compact_float(fit_row.get('M0'))}")
    if "D0_m2_ms" in fit_row:
        fragments.append(f"D0={compact_float(fit_row.get('D0_m2_ms'))} m2/ms")
    elif "D0_mm2_s" in fit_row:
        fragments.append(f"D0={compact_float(fit_row.get('D0_mm2_s'))} mm2/s")
    return fragments


def build_fit_label(fit_row: dict[str, object]) -> str:
    model = str(fit_row.get("model", "fit"))
    fragments = fit_parameter_fragments(fit_row)
    return ", ".join([model] + fragments) if fragments else model


def build_title(
    *,
    analysis_id: str | None,
    roi: str,
    direction: str,
    tn: object | None,
    signal_type: str | None,
    model: str | None,
) -> str:
    parts: list[str] = []
    if analysis_id:
        parts.append(str(analysis_id))
    parts.extend([f"ROI={roi}", f"direction={direction}"])
    if tn is not None:
        parts.append(f"TN={compact_float(tn)}")
    if signal_type:
        parts.append(f"type={signal_type}")
    if model:
        parts.append(f"model={model}")
    return " | ".join(parts)


def plot_nogse_signal_group(
    *,
    avg_df: pd.DataFrame,
    std_df: pd.DataFrame | None,
    xcol: str,
    ycol: str,
    out_png: Path,
    analysis_id: str,
    roi: str,
    direction: str,
    signal_type: str | None,
    fit_row: dict[str, object] | None = None,
    fit_curve: np.ndarray | None = None,
    x_data: np.ndarray | None = None,
    y_data: np.ndarray | None = None,
    fit_points: int | None = None,
    data_label: str = "signal",
    connect_data: bool = True,
) -> None:
    if x_data is None or y_data is None:
        x_vals, y_vals, sigma = prepare_signal_series(avg_df, std_df, xcol=xcol, ycol=ycol)
    else:
        x_vals = np.asarray(x_data, dtype=float)
        y_vals = np.asarray(y_data, dtype=float)
        _x, _y, sigma = prepare_signal_series(avg_df, std_df, xcol=xcol, ycol=ycol)

    if x_vals.size == 0 or y_vals.size == 0:
        return

    model = str(fit_row.get("model")) if fit_row is not None and fit_row.get("model") is not None else None
    tn = fit_row.get("TN") if fit_row is not None else unique_scalar(avg_df, "TN")
    title = build_title(
        analysis_id="",
        roi=str(roi),
        direction=str(direction),
        tn=tn,
        signal_type=None if signal_type is None else str(signal_type),
        model=model,
    )

    fit_x = None
    fit_y = None
    fit_label = None
    if fit_curve is not None and fit_curve.size:
        fit_x = np.asarray(fit_curve[:, 0], dtype=float)
        fit_y = np.asarray(fit_curve[:, 1], dtype=float)
        if fit_row is not None:
            fit_label = build_fit_label(fit_row)

    highlight_x = None
    highlight_y = None
    highlight_label = None
    if fit_points is not None and fit_points > 0:
        k = min(int(fit_points), len(x_vals))
        highlight_x = np.asarray(x_vals[:k], dtype=float)
        highlight_y = np.asarray(y_vals[:k], dtype=float)
        highlight_label = f"fit first {k}"

    render_xy_plot(
        x=x_vals,
        y=y_vals,
        sigma=sigma,
        out_png=out_png,
        title=title,
        xlabel=xcol,
        ylabel=ycol,
        data_label=data_label,
        connect_data=connect_data,
        fit_x=fit_x,
        fit_y=fit_y,
        fit_label=fit_label,
        highlight_x=highlight_x,
        highlight_y=highlight_y,
        highlight_label=highlight_label,
    )


def plot_nogse_signal_table(
    df: pd.DataFrame,
    *,
    out_root: Path,
    analysis_id: str,
    xcol: str,
    ycol: str,
    stat: str,
    rois: list[str] | None,
    directions: list[str] | None,
) -> list[Path]:
    work = df.copy()
    if rois is not None:
        work = work[work["roi"].astype(str).isin(rois)].copy()
    if directions is not None:
        work = work[work["direction"].astype(str).isin(directions)].copy()

    out_dir = out_root / analysis_id if analysis_id else out_root
    ensure_dir(out_dir)

    # x sweeps are plotted against x (one curve per G or per direction), every other
    # NOGSE acquisition against ``xcol`` (one curve per td_ms or per direction).
    sweep = x_sweep_mask(work)
    out_paths: list[Path] = []
    if (~sweep).any():
        out_paths.extend(_plot_g_curves(work.loc[~sweep], out_dir=out_dir, xcol=xcol, ycol=ycol, stat=stat))
    if sweep.any():
        out_paths.extend(_plot_x_sweep_curves(work.loc[sweep], out_dir=out_dir, ycol=ycol, stat=stat))
    return out_paths


def _curve_groups(merged: pd.DataFrame) -> tuple[list[str], list[str]]:
    type_group = []
    if "type" in merged.columns and merged["type"].notna().any():
        # A Hahn acquisition and a CPMG acquisition are different experiments,
        # even when all their other curve coordinates coincide.
        merged["type"] = merged["type"].fillna("unspecified").astype(str)
        type_group = ["type"]
    session_group = ["session_id"] if "session_id" in merged.columns and merged["session_id"].notna().any() else []
    return session_group, type_group


def _plot_x_sweep_curves(df: pd.DataFrame, *, out_dir: Path, ycol: str, stat: str) -> list[Path]:
    # The G series of one sweep share session, type, td_ms and b_step: G keeps their avg/std apart.
    merged = _prepare_avg_std(df, xcol="x", ycol=ycol, stat=stat, extra_keys=["G"])
    merged["x"] = pd.to_numeric(merged["x_raw"], errors="coerce")
    merged["G"] = pd.to_numeric(merged["G"], errors="coerce")
    merged = merged.dropna(subset=["x", "G"])  # the b0 volumes have no x
    session_group, type_group = _curve_groups(merged)
    base_cols = ["subj", *session_group, *type_group, "roi", "N", "td_ms"]

    xlabel = "NOGSE period x [ms]"
    ylabel = f"NOGSE signal ({ycol})"
    out_paths = _plot_series_by(
        merged=merged,
        out_dir=out_dir,
        group_cols=[*base_cols, "direction"],
        series_col="G",
        file_prefix="NOGSE_vs_x_by-G",
        xlabel=xlabel,
        ylabel=ylabel,
        stat=stat,
    )
    out_paths.extend(
        _plot_series_by(
            merged=merged,
            out_dir=out_dir,
            group_cols=[*base_cols, "G"],
            series_col="direction",
            file_prefix="NOGSE_vs_x_by-direction",
            xlabel=xlabel,
            ylabel=ylabel,
            stat=stat,
        )
    )
    return out_paths


def _plot_g_curves(df: pd.DataFrame, *, out_dir: Path, xcol: str, ycol: str, stat: str) -> list[Path]:
    merged = _prepare_avg_std(df, xcol=xcol, ycol=ycol, stat=stat)
    merged["x"] = pd.to_numeric(merged["x_raw"], errors="coerce")
    session_group, type_group = _curve_groups(merged)

    xlabel = f"Modulation gradient strength G [mT/m] ({xcol})"
    ylabel = f"NOGSE signal ({ycol})"

    out_paths = _plot_series_by(
        merged=merged,
        out_dir=out_dir,
        group_cols=["subj", *session_group, *type_group, "roi", "N", "direction"],
        series_col="td_ms",
        file_prefix="NOGSE_vs_G_by-td",
        xlabel=xlabel,
        ylabel=ylabel,
        stat=stat,
    )
    out_paths.extend(
        _plot_series_by(
            merged=merged,
            out_dir=out_dir,
            group_cols=["subj", *session_group, *type_group, "roi", "N", "td_ms"],
            series_col="direction",
            file_prefix="NOGSE_vs_G_by-direction",
            xlabel=xlabel,
            ylabel=ylabel,
            stat=stat,
        )
    )

    return out_paths
