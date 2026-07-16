from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotting.core import (
    XYSeries,
    compact_float,
    ensure_dir,
    render_multi_series_plot,
    render_xy_plot,
    sanitize_token,
)


SIGNAL_KEYS = ["subj", "roi", "N", "td_ms", "direction", "b_step"]


def load_long_parquet(path: str | Path) -> pd.DataFrame:
    return pd.read_parquet(path)


def compute_gradient_axis(values: pd.Series, *, xcol: str) -> np.ndarray:
    del xcol
    return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)


def _prepare_avg_std(df: pd.DataFrame, *, xcol: str, ycol: str, stat: str) -> pd.DataFrame:
    avg = df[df["stat"].astype(str) == str(stat)].copy()
    std = df[df["stat"].astype(str) == "std"].copy()

    if avg.empty:
        raise ValueError(f"No rows with stat={stat!r} were found in the input table.")

    avg = avg[SIGNAL_KEYS + [xcol, ycol]].rename(columns={ycol: "y_mean", xcol: "x_raw"})
    if std.empty:
        avg["y_std"] = np.nan
        merged = avg
    else:
        std = std[SIGNAL_KEYS + ["value"]].rename(columns={"value": "y_std"})
        merged = avg.merge(std, on=SIGNAL_KEYS, how="left")

    merged["y_mean"] = pd.to_numeric(merged["y_mean"], errors="coerce")
    merged["y_std"] = pd.to_numeric(merged["y_std"], errors="coerce")
    return merged.sort_values(["subj", "roi", "N", "direction", "td_ms", "b_step"], kind="stable")


def _finite_sigma_or_none(values: pd.Series) -> np.ndarray | None:
    sigma = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return sigma if np.isfinite(sigma).any() else None


def _format_group_value(value: object, *, digits: int = 3) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    return compact_float(numeric, digits=digits)


def _series_sort_key(value: object) -> tuple[int, float | str]:
    try:
        return (0, float(value))
    except (TypeError, ValueError):
        return (1, str(value))


def _plot_token(parts: dict[str, object]) -> str:
    return "_".join(f"{sanitize_token(key)}-{sanitize_token(_format_group_value(value))}" for key, value in parts.items())


def _plot_title(parts: dict[str, object], *, stat: str) -> str:
    title_parts = [f"{key}={_format_group_value(value)}" for key, value in parts.items()]
    title_parts.append(f"stat={stat}")
    return " | ".join(title_parts)


def _series_label(column: str, value: object) -> str:
    digits = 0 if column == "N" else 3
    return f"{column}={_format_group_value(value, digits=digits)}"


def _make_series(group: pd.DataFrame, *, label: str) -> XYSeries | None:
    curve = group.dropna(subset=["x", "y_mean"]).sort_values(["x", "b_step"], kind="stable")
    if curve.empty:
        return None
    return XYSeries(
        x=curve["x"].to_numpy(dtype=float),
        y=curve["y_mean"].to_numpy(dtype=float),
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
    ylim: tuple[float, float] | None,
) -> list[Path]:
    out_paths: list[Path] = []
    for keys, group in merged.groupby(group_cols, sort=False):
        key_tuple = keys if isinstance(keys, tuple) else (keys,)
        group_values = dict(zip(group_cols, key_tuple))
        series: list[XYSeries] = []

        values = sorted(group[series_col].drop_duplicates().tolist(), key=_series_sort_key)
        for value in values:
            sub = group[group[series_col].eq(value)].copy()
            curve = _make_series(sub, label=_series_label(series_col, value))
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
            ylim=ylim,
        )
        out_paths.append(out_png)

    return out_paths


def plot_ogse_signal_summary(
    df: pd.DataFrame,
    out_dir: str | Path,
    *,
    xcol: str = "g_thorsten",
    ycol: str = "value_norm",
    stat: str = "avg",
    ylim: tuple[float, float] | None = (0.0, 1.0),
) -> list[Path]:
    out_dir = Path(out_dir)
    ensure_dir(out_dir)

    merged = _prepare_avg_std(df, xcol=xcol, ycol=ycol, stat=stat)
    merged["x"] = compute_gradient_axis(merged["x_raw"], xcol=xcol)

    out_paths: list[Path] = []

    y_limits = None if ycol == "value" else ylim
    xlabel = f"Modulation gradient strength G [mT/m] ({xcol})"
    ylabel = f"OGSE signal ({ycol})"

    # A) For fixed subject, ROI, N and direction: one curve per td_ms.
    out_paths.extend(
        _plot_series_by(
            merged=merged,
            out_dir=out_dir,
            group_cols=["subj", "roi", "N", "direction"],
            series_col="td_ms",
            file_prefix="OGSE_vs_G_by-td",
            xlabel=xlabel,
            ylabel=ylabel,
            stat=stat,
            ylim=y_limits,
        )
    )

    # B) For fixed subject, ROI, N and td_ms: one curve per direction.
    out_paths.extend(
        _plot_series_by(
            merged=merged,
            out_dir=out_dir,
            group_cols=["subj", "roi", "N", "td_ms"],
            series_col="direction",
            file_prefix="OGSE_vs_G_by-direction",
            xlabel=xlabel,
            ylabel=ylabel,
            stat=stat,
            ylim=y_limits,
        )
    )

    return out_paths


def build_ogse_signal_fit_label(fit_row: dict[str, object]) -> str:
    model = str(fit_row.get("model", "unknown"))
    parts = [model]
    if "M0" in fit_row:
        parts.append(f"M0={compact_float(fit_row.get('M0'))}")
    if "D0_mm2_s" in fit_row:
        parts.append(f"D0={compact_float(fit_row.get('D0_mm2_s'))} mm2/s")
    elif "D0_m2_ms" in fit_row:
        parts.append(f"D0={compact_float(fit_row.get('D0_m2_ms'))} m2/ms")
    return ", ".join(parts)


def build_ogse_signal_fit_title(*, roi: str, direction: str, td_txt: str, n_txt: str, model: str) -> str:
    return (
        f"OGSE signal fit | ROI={roi} | direction={direction} | "
        f"model={model} | td_ms={td_txt} | N={n_txt}"
    )


def plot_ogse_signal_fit(
    *,
    x: np.ndarray,
    y: np.ndarray,
    fit_row: dict[str, object],
    out_png: Path,
    ycol: str,
    x_label: str,
    fit_points: int,
    fit_x: np.ndarray | None,
    fit_y: np.ndarray | None,
) -> None:
    roi = str(fit_row.get("roi", "roi"))
    direction = str(fit_row.get("direction", "direction"))
    td_txt = compact_float(fit_row.get("td_ms"))
    n_txt = compact_float(fit_row.get("N"), digits=0)
    model = str(fit_row.get("model", "unknown"))
    k = max(0, int(fit_points))
    fit_label = build_ogse_signal_fit_label(fit_row)
    render_xy_plot(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        out_png=out_png,
        title=build_ogse_signal_fit_title(
            roi=roi,
            direction=direction,
            td_txt=td_txt,
            n_txt=n_txt,
            model=model,
        ),
        xlabel=str(x_label),
        ylabel=ycol,
        data_label="data",
        connect_data=False,
        fit_x=None if fit_x is None else np.asarray(fit_x, dtype=float),
        fit_y=None if fit_y is None else np.asarray(fit_y, dtype=float),
        fit_label=fit_label,
        highlight_x=np.asarray(x[:k], dtype=float) if k > 0 else None,
        highlight_y=np.asarray(y[:k], dtype=float) if k > 0 else None,
        highlight_label=f"fit first {k}",
        yscale="log",
    )
