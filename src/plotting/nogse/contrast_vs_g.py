from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from plotting.contrast_vs_g import (
    build_display_model_name,
    canonical_xcol,
    canonical_ycol,
    filter_stat,
    plot_contrast_summary,
)
from plotting.core import compact_float, render_xy_plot

FAMILY_LABEL = "NOGSE contrast"


def fit_parameter_fragments(fit_row: dict[str, object]) -> list[str]:
    model = str(fit_row.get("model", ""))
    fragments: list[str] = []
    if model == "tort" and "alpha" in fit_row:
        fragments.append(f"alpha={compact_float(fit_row.get('alpha'))}")
    if model == "rest" and "tc_ms" in fit_row:
        fragments.append(f"tc_ms={compact_float(fit_row.get('tc_ms'))}")
    if model == "nogse_free_grad_offset" and "g0_mTm" in fit_row:
        fragments.append(f"g0={compact_float(fit_row.get('g0_mTm'))} mT/m")

    if "M0" in fit_row:
        fragments.append(f"M0={compact_float(fit_row.get('M0'))}")

    if "D0_m2_ms" in fit_row:
        fragments.append(f"D0={compact_float(fit_row.get('D0_m2_ms'))} m2/ms")
    elif "D0_mm2_s" in fit_row:
        fragments.append(f"D0={compact_float(fit_row.get('D0_mm2_s'))} mm2/s")

    return fragments


def build_contrast_fit_label(fit_row: dict[str, object], *, family_label: str | None = None) -> str:
    model = build_display_model_name(fit_row.get("model", "fit"), family_label=family_label)
    fragments = fit_parameter_fragments(fit_row)
    return ", ".join([model] + fragments) if fragments else model


def plot_nogse_contrast_summary(
    df: pd.DataFrame,
    *,
    out_root: str | Path,
    exp_id: str,
    xcol: str,
    ycol: str,
    directions: list[str],
    rois_requested: list[str] | None,
    stat: str,
) -> list[Path]:
    return plot_contrast_summary(
        df,
        out_root=out_root,
        exp_id=exp_id,
        xcol=xcol,
        ycol=ycol,
        directions=directions,
        rois_requested=rois_requested,
        family_label=FAMILY_LABEL,
    )


def plot_nogse_contrast_fit(
    *,
    x: np.ndarray,
    y: np.ndarray,
    fit_x: np.ndarray | None,
    fit_y: np.ndarray | None,
    fit_row: dict[str, object],
    out_png: Path,
    x_label: str,
    ycol: str,
) -> None:
    roi = str(fit_row.get("roi", "roi"))
    direction = str(fit_row.get("direction", "direction"))
    td_txt = compact_float(fit_row.get("td_ms"))
    n1_txt = compact_float(fit_row.get("N_1"), digits=0)
    n2_txt = compact_float(fit_row.get("N_2"), digits=0)
    display_model = build_display_model_name(fit_row.get("model", "fit"), family_label=FAMILY_LABEL)
    fit_label = build_contrast_fit_label(fit_row, family_label=FAMILY_LABEL)

    render_xy_plot(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        out_png=out_png,
        title=(
            f"ROI={roi} | direction={direction} | "
            f"model={display_model} | td_ms={td_txt} | N{n1_txt}-N{n2_txt}"
        ),
        xlabel=str(x_label),
        ylabel=ycol,
        data_label="data",
        connect_data=False,
        fit_x=None if fit_x is None else np.asarray(fit_x, dtype=float),
        fit_y=None if fit_y is None else np.asarray(fit_y, dtype=float),
        fit_label=fit_label,
    )


__all__ = [
    "FAMILY_LABEL",
    "build_contrast_fit_label",
    "canonical_xcol",
    "canonical_ycol",
    "filter_stat",
    "plot_nogse_contrast_fit",
    "plot_nogse_contrast_summary",
]
