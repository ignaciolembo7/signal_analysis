from __future__ import annotations

from pathlib import Path
import re
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tools.brain_labels import canonical_sheet_name, infer_subj_label


def _sanitize_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(value).strip())


def _join_sorted_unique(values: pd.Series) -> str:
    items = sorted({str(x).strip() for x in values.dropna() if str(x).strip()})
    return "|".join(items)


def _join_limited_sorted_unique(values: pd.Series, *, max_items: int = 6) -> str:
    items = sorted({str(x).strip() for x in values.dropna() if str(x).strip()})
    if len(items) <= max_items:
        return "|".join(items)
    return "|".join(items[:max_items] + [f"...(+{len(items) - max_items})"])


def _rms_numeric(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(np.square(arr))))


def _valid_text(value: object) -> str | None:
    token = str(value).strip()
    if not token or token.lower() in {"nan", "none", "<na>"}:
        return None
    return token


def _infer_sheet(path: Path, df: pd.DataFrame) -> str:
    if "source_file" in df.columns:
        vals = pd.Series(df["source_file"]).dropna().astype(str).unique().tolist()
        if len(vals) == 1:
            sheet = canonical_sheet_name(vals[0])
            if sheet:
                return sheet

    parent = path.parent.parent.name
    if parent:
        return str(parent)

    return canonical_sheet_name(path.parent.name) or path.parent.name


def discover_fit_param_files(root: str | Path, pattern: str = "**/fit_params*.parquet") -> list[Path]:
    base = Path(root)
    if base.is_file():
        return [base]
    files = sorted(base.glob(pattern))
    if not files:
        raise FileNotFoundError(f"Could not find fit_params with pattern={pattern!r} under {base}")
    return files


def load_monoexp_fit_measurements(
    fits_root: str | Path,
    *,
    pattern: str = "**/fit_params*.parquet",
    subjs: Sequence[str] | None = None,
    rois: Sequence[str] | None = None,
    directions: Sequence[str] | None = None,
    Ns: Sequence[float | str] | None = None,
    stat: str = "avg",
    ycol: str = "value_norm",
) -> pd.DataFrame:
    subj_set = {str(x) for x in subjs} if subjs is not None else None
    roi_set = {str(x) for x in rois} if rois is not None else None
    direction_set = {str(x) for x in directions} if directions is not None else None
    n_set = {float(x) for x in Ns} if Ns is not None else None

    frames: list[pd.DataFrame] = []
    for path in discover_fit_param_files(fits_root, pattern=pattern):
        df = pd.read_parquet(path)
        if df.empty:
            continue

        out = df.copy()
        if "fit_kind" in out.columns:
            out = out[out["fit_kind"].astype(str) == "monoexp"].copy()
        if "model" in out.columns:
            out = out[out["model"].astype(str) == "monoexp"].copy()
        if "ok" in out.columns:
            out = out[out["ok"].astype(bool)].copy()
        if "stat" in out.columns:
            out = out[out["stat"].astype(str) == str(stat)].copy()
        if "ycol" in out.columns:
            out = out[out["ycol"].astype(str) == str(ycol)].copy()
        if out.empty:
            continue

        out["roi"] = out["roi"].astype(str)
        out["direction"] = out["direction"].astype(str)
        out["td_ms"] = pd.to_numeric(out.get("td_ms", np.nan), errors="coerce")
        out["Delta_app_ms"] = pd.to_numeric(out.get("Delta_app_ms", np.nan), errors="coerce")
        out["N"] = pd.to_numeric(out.get("N", np.nan), errors="coerce")
        out["D0_mm2_s"] = pd.to_numeric(out.get("D0_mm2_s", np.nan), errors="coerce")
        out["D0_err_mm2_s"] = pd.to_numeric(out.get("D0_err_mm2_s", np.nan), errors="coerce")

        sheet = _infer_sheet(path, out)
        out["sheet"] = sheet

        if "source_file" in out.columns:
            source_vals = pd.Series(out["source_file"]).dropna().astype(str).unique().tolist()
            source_name = source_vals[0] if len(source_vals) == 1 else path.name
        else:
            source_name = path.name
            out["source_file"] = source_name

        if "subj" in out.columns:
            subj_series = out["subj"].map(_valid_text)
            out["subj"] = subj_series.where(subj_series.notna(), infer_subj_label(sheet, source_name=source_name))
        else:
            out["subj"] = infer_subj_label(sheet, source_name=source_name)

        out["fit_file"] = str(path)

        if subj_set is not None:
            out = out[out["subj"].isin(subj_set)]
        if roi_set is not None:
            out = out[out["roi"].isin(roi_set)]
        if direction_set is not None:
            out = out[out["direction"].isin(direction_set)]
        if n_set is not None:
            out = out[out["N"].isin(n_set)]

        out = out.dropna(subset=["D0_mm2_s"])
        if out.empty:
            continue
        frames.append(out)

    if not frames:
        raise ValueError("No monoexp fits remained after filtering.")

    cols = [
        "subj",
        "sheet",
        "roi",
        "direction",
        "N",
        "td_ms",
        "Delta_app_ms",
        "D0_mm2_s",
        "D0_err_mm2_s",
        "n_fit",
        "source_file",
        "fit_file",
    ]
    df_all = pd.concat(frames, ignore_index=True)
    for col in cols:
        if col not in df_all.columns:
            df_all[col] = np.nan
    return df_all[cols].copy()


def aggregate_monoexp_by_x(df: pd.DataFrame, *, xcol: str) -> pd.DataFrame:
    if xcol not in {"td_ms", "Delta_app_ms"}:
        raise ValueError(f"Invalid xcol: {xcol}")

    work = df.copy()
    work[xcol] = pd.to_numeric(work[xcol], errors="coerce")
    work["N"] = pd.to_numeric(work.get("N", np.nan), errors="coerce")
    work = work.dropna(subset=[xcol, "D0_mm2_s"])
    if work.empty:
        raise ValueError(f"No valid data remained for xcol={xcol}.")

    other_x = "Delta_app_ms" if xcol == "td_ms" else "td_ms"
    grouped = (
        work.groupby(["subj", "sheet", "roi", "direction", "N", xcol], as_index=False)
        .agg(
            D_mean_mm2_s=("D0_mm2_s", "mean"),
            D_std_mm2_s=("D0_err_mm2_s", _rms_numeric),
            n_measurements=("D0_mm2_s", "count"),
            other_x_mean=(other_x, "mean"),
            source_files=("source_file", _join_sorted_unique),
        )
    )
    grouped["D_std_mm2_s"] = grouped["D_std_mm2_s"].fillna(0.0)
    return grouped.sort_values(["subj", "sheet", "roi", "direction", "N", xcol], kind="stable").reset_index(drop=True)


def _plot_group_curves(
    df: pd.DataFrame,
    *,
    group_key: tuple[str, str],
    curve_col: str,
    xcol: str,
    out_path: Path,
    title: str,
    xlabel: str,
) -> Path:
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.get_cmap("tab10")

    curve_vals = sorted(df[curve_col].astype(str).unique().tolist())
    for idx, curve_val in enumerate(curve_vals):
        sub = df[df[curve_col].astype(str) == str(curve_val)].sort_values(xcol, kind="stable")
        x = pd.to_numeric(sub[xcol], errors="coerce").to_numpy(float)
        y = pd.to_numeric(sub["D_mean_mm2_s"], errors="coerce").to_numpy(float)
        yerr = pd.to_numeric(sub["D_std_mm2_s"], errors="coerce").to_numpy(float)
        color = colors(idx / max(1, len(curve_vals) - 1))
        ax.plot(x, y, marker="o", linewidth=2, label=str(curve_val), color=color)
        if np.isfinite(yerr).any():
            low = y - np.nan_to_num(yerr, nan=0.0)
            high = y + np.nan_to_num(yerr, nan=0.0)
            ax.fill_between(x, low, high, color=color, alpha=0.2)

    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(r"D [mm$^2$/s]", fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend(title=curve_col, fontsize=9)
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


def _compare_n_curve_score(df: pd.DataFrame, *, xcol: str) -> int:
    work = df.copy()
    work["N"] = pd.to_numeric(work.get("N", np.nan), errors="coerce")
    work[xcol] = pd.to_numeric(work.get(xcol, np.nan), errors="coerce")
    work = work.dropna(subset=["N", xcol])
    if work.empty:
        return 0

    n_x_counts = work.groupby("N")[xcol].nunique()
    if len(n_x_counts) < 2:
        return 0
    return int(n_x_counts.max())


def _collapse_compare_n_scope(df: pd.DataFrame, *, xcol: str) -> pd.DataFrame:
    group_cols = ["subj", "roi", "direction", "N", xcol]
    missing = [col for col in group_cols if col not in df.columns]
    if missing:
        raise KeyError(f"Missing columns for compare_N collapse: {missing}")

    work = df.copy()
    work["D_mean_mm2_s"] = pd.to_numeric(work.get("D_mean_mm2_s", np.nan), errors="coerce")
    work["D_std_mm2_s"] = pd.to_numeric(work.get("D_std_mm2_s", np.nan), errors="coerce")
    work["n_measurements"] = pd.to_numeric(work.get("n_measurements", 1), errors="coerce").fillna(1)

    if "source_files" not in work.columns:
        work["source_files"] = ""
    if "sheet" not in work.columns:
        work["sheet"] = ""
    if "other_x_mean" not in work.columns:
        work["other_x_mean"] = np.nan

    collapsed = (
        work.groupby(group_cols, as_index=False)
        .agg(
            sheet=("sheet", _join_limited_sorted_unique),
            D_mean_mm2_s=("D_mean_mm2_s", "mean"),
            D_std_mm2_s=("D_std_mm2_s", _rms_numeric),
            n_measurements=("n_measurements", "sum"),
            other_x_mean=("other_x_mean", "mean"),
            source_files=("source_files", _join_limited_sorted_unique),
        )
    )
    collapsed["D_std_mm2_s"] = collapsed["D_std_mm2_s"].fillna(0.0)
    return collapsed.sort_values(group_cols, kind="stable").reset_index(drop=True)


def plot_compare_roi(df_avg: pd.DataFrame, *, xcol: str, out_dir: str | Path) -> list[Path]:
    out_dir = Path(out_dir)
    xlabel = "td_ms [ms]" if xcol == "td_ms" else r"$\Delta_{app}$ [ms]"
    outputs: list[Path] = []
    for (subj, Nval, roi), sub in df_avg.groupby(["subj", "N", "roi"], sort=True):
        out_path = out_dir / f"{xcol}" / "compare_roi" / f"{_sanitize_token(subj)}__N={_sanitize_token(Nval)}__{_sanitize_token(roi)}.png"
        outputs.append(
            _plot_group_curves(
                sub,
                group_key=(subj, roi),
                curve_col="direction",
                xcol=xcol,
                out_path=out_path,
                title=f"{subj} | N={Nval:g} | {roi} | D vs {xcol}",
                xlabel=xlabel,
            )
        )
    return outputs


def plot_compare_direction(df_avg: pd.DataFrame, *, xcol: str, out_dir: str | Path) -> list[Path]:
    out_dir = Path(out_dir)
    xlabel = "td_ms [ms]" if xcol == "td_ms" else r"$\Delta_{app}$ [ms]"
    outputs: list[Path] = []
    for (subj, Nval, direction), sub in df_avg.groupby(["subj", "N", "direction"], sort=True):
        out_path = out_dir / f"{xcol}" / "compare_direction" / f"{_sanitize_token(subj)}__N={_sanitize_token(Nval)}__dir={_sanitize_token(direction)}.png"
        outputs.append(
            _plot_group_curves(
                sub,
                group_key=(subj, direction),
                curve_col="roi",
                xcol=xcol,
                out_path=out_path,
                title=f"{subj} | N={Nval:g} | dir={direction} | D vs {xcol}",
                xlabel=xlabel,
            )
        )
    return outputs


plot_by_roi = plot_compare_roi
plot_by_direction = plot_compare_direction


def plot_compare_N_within_sheet(df_avg: pd.DataFrame, *, xcol: str, out_dir: str | Path) -> list[Path]:
    out_dir = Path(out_dir)
    xlabel = "td_ms [ms]" if xcol == "td_ms" else r"$\Delta_{app}$ [ms]"
    outputs: list[Path] = []
    for (subj, roi, direction), subj_sub in df_avg.groupby(["subj", "roi", "direction"], sort=True):
        collapsed = _collapse_compare_n_scope(subj_sub, xcol=xcol)
        subject_score = _compare_n_curve_score(collapsed, xcol=xcol)
        if subject_score < 2:
            continue

        plotted_sheet_scope = False
        for sheet, sheet_sub in subj_sub.groupby("sheet", sort=True):
            if _compare_n_curve_score(sheet_sub, xcol=xcol) < subject_score:
                continue
            plotted_sheet_scope = True
            out_path = (
                out_dir
                / f"{xcol}"
                / "compare_N"
                / f"{_sanitize_token(subj)}__{_sanitize_token(sheet)}__{_sanitize_token(roi)}__dir={_sanitize_token(direction)}.png"
            )
            outputs.append(
                _plot_group_curves(
                    sheet_sub,
                    group_key=(subj, direction),
                    curve_col="N",
                    xcol=xcol,
                    out_path=out_path,
                    title=f"{subj} | {sheet} | {roi} | dir={direction} | compare N",
                    xlabel=xlabel,
                )
            )

        if plotted_sheet_scope:
            continue

        out_path = (
            out_dir
            / f"{xcol}"
            / "compare_N"
            / f"{_sanitize_token(subj)}__all-sheets__{_sanitize_token(roi)}__dir={_sanitize_token(direction)}.png"
        )
        outputs.append(
            _plot_group_curves(
                collapsed,
                group_key=(subj, direction),
                curve_col="N",
                xcol=xcol,
                out_path=out_path,
                title=f"{subj} | {roi} | dir={direction} | compare N",
                xlabel=xlabel,
            )
        )
    return outputs
