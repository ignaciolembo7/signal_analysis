from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotting.core import XYSeries, distinct_colors, ensure_dir, render_multi_series_plot

VALID_YCOLS = {"value", "value_norm"}


def canonical_xcol(xcol: str) -> str:
    return str(xcol)


def canonical_ycol(ycol: str) -> str:
    y = str(ycol)
    if y not in VALID_YCOLS:
        raise ValueError(f"Unrecognized y column {y!r}. Allowed values: {sorted(VALID_YCOLS)}.")
    return y


def plot_x_series(df: pd.DataFrame, xcol: str) -> np.ndarray:
    return pd.to_numeric(df[xcol], errors="coerce").to_numpy(dtype=float)


def require_columns(df: pd.DataFrame, cols: list[str], *, label: str) -> None:
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise KeyError(f"{label}: missing required columns {missing}. Columns={sorted(df.columns)}")


def filter_stat(df: pd.DataFrame, stat: str | None) -> pd.DataFrame:
    if "stat" not in df.columns or stat is None or str(stat).upper() == "ALL":
        return df.copy()
    out = df[df["stat"].astype(str) == str(stat)].copy()
    if out.empty:
        available = sorted(df["stat"].dropna().astype(str).unique().tolist())
        raise ValueError(f"No rows remain after filtering stat={stat!r}. Available stats: {available}")
    return out


def unique_int_or_none(df: pd.DataFrame, col: str) -> int | None:
    if col not in df.columns:
        return None
    vals = pd.to_numeric(df[col], errors="coerce").dropna().unique()
    if len(vals) != 1:
        return None
    return int(vals[0])


def selected_rois(df: pd.DataFrame, requested: list[str] | None) -> list[str]:
    available = sorted(df["roi"].dropna().astype(str).unique().tolist())
    if not requested:
        return available
    chosen = [roi for roi in requested if roi in available]
    missing = [roi for roi in requested if roi not in available]
    if missing:
        raise ValueError(f"ROIs not found: {missing}. Available ROIs: {available}")
    return chosen


def contrast_ylabel(*, ycol: str, n1: int | None, n2: int | None, family_label: str) -> str:
    if n1 is not None and n2 is not None:
        return rf"{family_label} $\Delta M_{{N{n1}-N{n2}}}$"
    return ycol


def prep_contrast_curve(df: pd.DataFrame, *, xcol: str, ycol: str) -> tuple[np.ndarray, np.ndarray]:
    x = plot_x_series(df, xcol)
    y = pd.to_numeric(df[ycol], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        return np.array([]), np.array([])
    x_valid = x[mask]
    y_valid = y[mask]
    order = np.argsort(x_valid)
    return x_valid[order], y_valid[order]


def plot_contrast_summary(
    df: pd.DataFrame,
    *,
    out_root: str | Path,
    exp_id: str,
    xcol: str,
    ycol: str,
    directions: list[str],
    rois_requested: list[str] | None,
    family_label: str,
) -> list[Path]:
    outdir = Path(out_root) / exp_id
    ensure_dir(outdir)

    xcol = canonical_xcol(xcol)
    ycol = canonical_ycol(ycol)
    require_columns(df, ["roi", "direction", xcol, ycol], label="plot_contrast_summary")

    rois = selected_rois(df, None)
    selected_subset = selected_rois(df, rois_requested)
    n1 = unique_int_or_none(df, "N_1")
    n2 = unique_int_or_none(df, "N_2")
    ylabel = contrast_ylabel(ycol=ycol, n1=n1, n2=n2, family_label=family_label)
    title_suffix = f"N{n1}-N{n2}" if (n1 is not None and n2 is not None) else exp_id

    outputs: list[Path] = []

    for direction in directions:
        direction_df = df[df["direction"].astype(str) == str(direction)].copy()
        if direction_df.empty:
            continue

        colors = distinct_colors(len(rois))
        series_all: list[XYSeries] = []
        for roi, color in zip(rois, colors):
            curve_df = direction_df[direction_df["roi"].astype(str) == str(roi)].copy()
            x_vals, y_vals = prep_contrast_curve(curve_df, xcol=xcol, ycol=ycol)
            if x_vals.size == 0:
                continue
            series_all.append(XYSeries(x=x_vals, y=y_vals, label=str(roi), color=color))

        if series_all:
            out_png = outdir / f"{exp_id}.{ycol}_vs_{xcol}.direction-{direction}.allROIs.png"
            render_multi_series_plot(
                series=series_all,
                out_png=out_png,
                title=f"direction={direction} | {title_suffix}",
                xlabel=f"Gradient strength ({xcol}) [mT/m]",
                ylabel=ylabel,
                legend_title="ROI",
            )
            outputs.append(out_png)

        if rois_requested:
            colors_subset = distinct_colors(len(selected_subset))
            series_subset: list[XYSeries] = []
            for roi, color in zip(selected_subset, colors_subset):
                curve_df = direction_df[direction_df["roi"].astype(str) == str(roi)].copy()
                x_vals, y_vals = prep_contrast_curve(curve_df, xcol=xcol, ycol=ycol)
                if x_vals.size == 0:
                    continue
                series_subset.append(XYSeries(x=x_vals, y=y_vals, label=str(roi), color=color))

            if series_subset:
                roi_tag = "-".join(selected_subset)
                out_png = outdir / f"{exp_id}.{ycol}_vs_{xcol}.direction-{direction}.subsetROIs-{roi_tag}.png"
                render_multi_series_plot(
                    series=series_subset,
                    out_png=out_png,
                    title=f"direction={direction} | {title_suffix}",
                    xlabel=f"Gradient strength ({xcol}) [mT/m]",
                    ylabel=ylabel,
                    legend_title="ROI subset",
                )
                outputs.append(out_png)

    for roi in rois:
        roi_df = df[(df["roi"].astype(str) == str(roi)) & (df["direction"].astype(str).isin(directions))].copy()
        if roi_df.empty:
            continue

        colors = distinct_colors(len(directions))
        series: list[XYSeries] = []
        for direction, color in zip(directions, colors):
            curve_df = roi_df[roi_df["direction"].astype(str) == str(direction)].copy()
            x_vals, y_vals = prep_contrast_curve(curve_df, xcol=xcol, ycol=ycol)
            if x_vals.size == 0:
                continue
            series.append(XYSeries(x=x_vals, y=y_vals, label=str(direction), color=color))

        if not series:
            continue

        out_png = outdir / f"{exp_id}.{ycol}_vs_{xcol}.ROI-{roi}.directions.png"
        render_multi_series_plot(
            series=series,
            out_png=out_png,
            title=f"ROI={roi} | {title_suffix}",
            xlabel=f"Gradient strength ({xcol}) [mT/m]",
            ylabel=ylabel,
            legend_title="direction",
        )
        outputs.append(out_png)

    return outputs


def contrast_family_prefix(family_label: str | None) -> str | None:
    if not family_label:
        return None
    first_token = str(family_label).strip().split(maxsplit=1)[0].lower()
    if first_token in {"nogse", "ogse"}:
        return first_token
    return None


def build_display_model_name(model: object, *, family_label: str | None) -> str:
    model_name = str(model or "fit").strip() or "fit"
    family_prefix = contrast_family_prefix(family_label)
    if family_prefix is None or model_name == "fit":
        return model_name
    if model_name.lower().startswith(f"{family_prefix} "):
        return model_name
    return f"{family_prefix} {model_name}"
