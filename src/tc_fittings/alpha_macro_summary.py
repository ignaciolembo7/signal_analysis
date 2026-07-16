from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from data_processing.io import write_xlsx_csv_outputs
from tools.brain_labels import infer_subj_label


DEFAULT_CC_REGIONS = ["PostCC", "MidPostCC", "CentralCC", "MidAntCC", "AntCC"]
DEFAULT_LATERAL_VENTRICLES = ["Left-Lateral-Ventricle", "Right-Lateral-Ventricle"]
DEFAULT_DIRECTION_ALIASES = {"x": "long", "y": "tra", "z": "tra"}
DEFAULT_SUBJ_COLORS = {"BRAIN": "#377eb8", "LUDG": "#ff7f00", "MBBL": "#4daf4a"}
REQUIRED_DPROJ_COLUMNS = {"roi", "direction", "bvalue", "D_proj"}


def _as_str(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip()


def _sanitize_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(value).strip())


def _merge_pipe_lists(values: Iterable[object]) -> str:
    items: set[str] = set()
    for value in values:
        for token in str(value).split("|"):
            token = token.strip()
            if token:
                items.add(token)
    return "|".join(sorted(items))


def _select_bvalue_and_bstep(
    bvalues: Sequence[float],
    *,
    selected_bstep: int | None,
    group_label: str,
) -> tuple[float, int]:
    if len(bvalues) == 0:
        raise ValueError(f"No bvalues available for {group_label}.")
    if selected_bstep is not None and selected_bstep < 1:
        raise ValueError("selected_bstep must be >= 1.")

    ordered = np.sort(np.asarray(bvalues, dtype=float))
    if selected_bstep is None:
        return float(ordered[-1]), int(len(ordered))
    if selected_bstep > len(ordered):
        raise ValueError(
            "selected_bstep="
            f"{selected_bstep} fuera de rango para {group_label}. "
            f"Solo hay {len(ordered)} bsteps disponibles: {list(map(float, ordered))}"
        )
    return float(ordered[selected_bstep - 1]), int(selected_bstep)


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=0)
    raise ValueError(f"Unsupported format for {path}")


def discover_dproj_files(root: str | Path, pattern: str = "**/*.Dproj.long.parquet") -> list[Path]:
    base = Path(root)
    if base.is_file():
        return [base]
    files = sorted(base.glob(pattern))
    if not files:
        raise FileNotFoundError(f"Could not find files with pattern={pattern!r} under {base}")
    return files


def parse_direction_aliases(items: Sequence[str] | None) -> dict[str, str]:
    aliases = dict(DEFAULT_DIRECTION_ALIASES)
    if not items:
        return aliases
    for raw in items:
        if "=" not in raw:
            raise ValueError(f"Invalid alias {raw!r}. Use source=target format, for example x=long.")
        source, target = raw.split("=", 1)
        source = source.strip()
        target = target.strip()
        if not source or not target:
            raise ValueError(f"Invalid alias {raw!r}. Use source=target format.")
        aliases[source] = target
    return aliases


def load_dproj_measurements(
    dproj_root: str | Path,
    *,
    pattern: str = "**/*.Dproj.long.parquet",
    subjs: Sequence[str] | None = None,
    rois: Sequence[str] | None = None,
    directions: Sequence[str] | None = None,
    N: float | None = None,
    Hz: float | None = None,
    bvalue_decimals: int = 1,
) -> pd.DataFrame:
    if N is not None and Hz is not None:
        raise ValueError("Choose only one of N and Hz for filtering.")

    frames: list[pd.DataFrame] = []
    for path in discover_dproj_files(dproj_root, pattern=pattern):
        df = _read_table(path)
        missing = REQUIRED_DPROJ_COLUMNS - set(df.columns)
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")

        out = df.copy()
        out["roi"] = _as_str(out["roi"])
        out["direction"] = _as_str(out["direction"])
        out["bvalue"] = pd.to_numeric(out["bvalue"], errors="coerce")
        out["D_proj"] = pd.to_numeric(out["D_proj"], errors="coerce")
        out["Delta_app_ms"] = pd.to_numeric(out.get("Delta_app_ms", np.nan), errors="coerce")
        out["delta_ms"] = pd.to_numeric(out.get("delta_ms", np.nan), errors="coerce")
        out["td_ms"] = pd.to_numeric(out.get("td_ms", np.nan), errors="coerce")
        out["N"] = pd.to_numeric(out.get("N", np.nan), errors="coerce")
        out["Hz"] = pd.to_numeric(out.get("Hz", np.nan), errors="coerce")
        out["sheet"] = _as_str(out["sheet"]) if "sheet" in out.columns else path.stem
        if "subj" in out.columns:
            subj = _as_str(out["subj"])
            invalid = subj.str.lower().isin({"", "nan", "none", "<na>"})
            out["subj"] = subj.mask(invalid, out["sheet"].map(lambda x: infer_subj_label(x, source_name=path.name)))
        else:
            out["subj"] = out["sheet"].map(lambda x: infer_subj_label(x, source_name=path.name))
        out["source_file"] = path.name

        frames.append(out)

    if not frames:
        raise FileNotFoundError(f"No D_proj files found in {dproj_root} with pattern={pattern!r}")

    df_all = pd.concat(frames, ignore_index=True)
    return load_dproj_measurements_from_table(
        df_all,
        subjs=subjs,
        rois=rois,
        directions=directions,
        N=N,
        Hz=Hz,
        bvalue_decimals=bvalue_decimals,
    )


def load_dproj_measurements_from_table(
    df: pd.DataFrame,
    *,
    subjs: Sequence[str] | None = None,
    rois: Sequence[str] | None = None,
    directions: Sequence[str] | None = None,
    N: float | None = None,
    Hz: float | None = None,
    bvalue_decimals: int = 1,
) -> pd.DataFrame:
    if N is not None and Hz is not None:
        raise ValueError("Choose only one of N and Hz for filtering.")

    missing = REQUIRED_DPROJ_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"D_proj table is missing required columns: {sorted(missing)}")

    subj_set = {str(x) for x in subjs} if subjs is not None else None
    roi_set = {str(x) for x in rois} if rois is not None else None
    direction_set = {str(x) for x in directions} if directions is not None else None

    out = df.copy()
    out["roi"] = _as_str(out["roi"])
    out["direction"] = _as_str(out["direction"])
    out["bvalue"] = pd.to_numeric(out["bvalue"], errors="coerce")
    out["D_proj"] = pd.to_numeric(out["D_proj"], errors="coerce")
    out["Delta_app_ms"] = pd.to_numeric(out.get("Delta_app_ms", np.nan), errors="coerce")
    out["delta_ms"] = pd.to_numeric(out.get("delta_ms", np.nan), errors="coerce")
    out["td_ms"] = pd.to_numeric(out.get("td_ms", np.nan), errors="coerce")
    out["N"] = pd.to_numeric(out.get("N", np.nan), errors="coerce")
    out["Hz"] = pd.to_numeric(out.get("Hz", np.nan), errors="coerce")
    out["sheet"] = _as_str(out["sheet"]) if "sheet" in out.columns else ""
    if "subj" in out.columns:
        subj = _as_str(out["subj"])
        invalid = subj.str.lower().isin({"", "nan", "none", "<na>"})
        out["subj"] = subj.mask(invalid, out["sheet"].map(lambda x: infer_subj_label(x, source_name="")))
    else:
        out["subj"] = out["sheet"].map(lambda x: infer_subj_label(x, source_name=""))
    if "source_file" not in out.columns:
        out["source_file"] = ""

    if subj_set is not None:
        out = out[out["subj"].isin(subj_set)]
    if roi_set is not None:
        out = out[out["roi"].isin(roi_set)]
    if direction_set is not None:
        out = out[out["direction"].isin(direction_set)]
    if N is not None:
        out = out[np.isclose(out["N"], float(N), atol=1e-6, equal_nan=False)]
    if Hz is not None:
        out = out[np.isclose(out["Hz"], float(Hz), atol=1e-6, equal_nan=False)]

    out = out.dropna(subset=["bvalue", "D_proj", "Delta_app_ms"])
    if out.empty:
        selector = f"N={N}" if N is not None else (f"Hz={Hz}" if Hz is not None else "no N/Hz filter")
        raise ValueError(f"No data remained after filtering ({selector}).")

    df_all = out
    df_all["bvalue"] = df_all["bvalue"].round(int(bvalue_decimals))

    group_cols = ["subj", "roi", "direction", "bvalue", "Delta_app_ms"]
    df_avg = (
        df_all.groupby(group_cols, as_index=False)
        .agg(
            D_mean_mm2_s=("D_proj", "mean"),
            D_std_mm2_s=("D_proj", "std"),
            n_measurements=("D_proj", "count"),
            Hz=("Hz", "mean"),
            N=("N", "mean"),
            delta_ms=("delta_ms", "mean"),
            td_ms=("td_ms", "mean"),
            sheet_count=("sheet", "nunique"),
            sheets=("sheet", lambda s: "|".join(sorted({str(x) for x in s.dropna()}))),
            source_files=("source_file", lambda s: "|".join(sorted({str(x) for x in s.dropna()}))),
        )
    )
    df_avg["D_std_mm2_s"] = df_avg["D_std_mm2_s"].fillna(0.0)
    return df_avg.sort_values(["subj", "roi", "direction", "bvalue", "Delta_app_ms"], kind="stable").reset_index(drop=True)


def plot_d_vs_delta_curves(
    df_avg: pd.DataFrame,
    *,
    out_dir: str | Path,
    selected_bstep: int | None = None,
    selected_bstep_by_group: Mapping[tuple[str, str, str], int] | None = None,
    reference_D0: float | None = None,
    reference_D0_error: float | None = None,
) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs: list[Path] = []
    for (subj, roi, direction), sub in df_avg.groupby(["subj", "roi", "direction"], sort=True):
        sub = sub.sort_values(["bvalue", "Delta_app_ms"], kind="stable")
        bvalues = sorted(sub["bvalue"].dropna().unique().tolist())
        if not bvalues:
            continue

        fig, ax = plt.subplots(figsize=(8, 6))
        cmap = plt.get_cmap("viridis")
        group_key = (str(subj), str(roi), str(direction))
        group_selected_bstep = selected_bstep
        if selected_bstep_by_group is not None:
            roi_norm_key = (str(subj), str(roi).replace("_norm", "").lower(), str(direction))
            if group_key in selected_bstep_by_group:
                group_selected_bstep = int(selected_bstep_by_group[group_key])
            elif roi_norm_key in selected_bstep_by_group:
                group_selected_bstep = int(selected_bstep_by_group[roi_norm_key])
        selected_bvalue, chosen_bstep = _select_bvalue_and_bstep(
            bvalues,
            selected_bstep=group_selected_bstep,
            group_label=f"subj={subj}, roi={roi}, direction={direction}",
        )
        selected = sub[np.isclose(sub["bvalue"], selected_bvalue, atol=1e-9)].sort_values("Delta_app_ms")

        for idx, bvalue in enumerate(bvalues):
            curve = sub[np.isclose(sub["bvalue"], bvalue, atol=1e-9)].sort_values("Delta_app_ms")
            if curve.empty:
                continue
            color = cmap(idx / max(1, len(bvalues) - 1))
            ax.plot(
                curve["Delta_app_ms"],
                curve["D_mean_mm2_s"],
                marker="o",
                linewidth=2,
                color=color,
                label=f"b-value = {bvalue:g} s/mm$^2$",
            )

        if not selected.empty:
            d0_mean = float(np.nanmean(selected["D_mean_mm2_s"].to_numpy(float)))
            d0_std = float(np.nanstd(selected["D_mean_mm2_s"].to_numpy(float)))
            label = f"D (bstep={chosen_bstep}, b={selected_bvalue:g}) = {d0_mean:.6f}"
            if reference_D0 is not None and np.isfinite(reference_D0) and reference_D0 > 0:
                alpha = d0_mean / float(reference_D0)
                if reference_D0_error is not None and np.isfinite(reference_D0_error):
                    rel_d = d0_std / d0_mean if d0_mean > 0 else np.nan
                    rel_ref = float(reference_D0_error) / float(reference_D0)
                    alpha_err = np.sqrt(rel_d**2 + rel_ref**2) * alpha if np.isfinite(rel_d) else np.nan
                    label = f"{label}\n$\\alpha$ = ({alpha:.4f} $\\pm$ {alpha_err:.4f})"
                else:
                    label = f"{label}\n$\\alpha$ = ({alpha:.4f})"
            ax.axhline(d0_mean, color="black", linestyle="--", linewidth=1.1, alpha=0.7, label=label)

        if reference_D0 is not None and np.isfinite(reference_D0):
            title = f"{subj} | {roi} | dir={direction} | $D_0$ = {float(reference_D0):.4f} mm$^2$/s"
        else:
            title = f"{subj} | {roi} | dir={direction}"
        ax.set_title(title, fontsize=15)
        ax.set_xlabel(r"$\Delta_{app}$ [ms]", fontsize=13)
        ax.set_ylabel(r"D [mm$^2$/s]", fontsize=13)
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.legend(fontsize=9, title="b-value")
        plt.tight_layout()

        out_path = out_dir / f"D_vs_delta_app_{_sanitize_token(subj)}_{_sanitize_token(roi)}_dir={_sanitize_token(direction)}.png"
        fig.savefig(out_path, dpi=300)
        plt.close(fig)
        outputs.append(out_path)

    return outputs


def _expand_direction_alias_rows(df_summary: pd.DataFrame, aliases: dict[str, str]) -> pd.DataFrame:
    if df_summary.empty:
        return df_summary.copy()

    raw = df_summary.copy()
    raw["direction_kind"] = "raw"
    raw["direction_components"] = raw["direction"]

    derived_rows: list[dict[str, object]] = []
    for (subj, roi, alias), sub in (
        raw.assign(direction_alias=raw["direction"].map(lambda x: aliases.get(str(x), str(x))))
        .groupby(["subj", "roi", "direction_alias"], sort=True)
    ):
        directions = sorted(sub["direction"].astype(str).unique().tolist())
        if len(directions) == 1 and directions[0] == alias:
            continue
        derived_rows.append(
            {
                "subj": subj,
                "sheet": _merge_pipe_lists(sub["sheet"].dropna()),
                "roi": roi,
                "direction": alias,
                "direction_kind": "derived",
                "direction_components": "|".join(directions),
                "alpha_macro": float(sub["alpha_macro"].mean()),
                "alpha_macro_error": float(sub["alpha_macro_error"].mean()),
                "D0_mean_mm2_s": float(sub["D0_mean_mm2_s"].mean()),
                "D0_std_mm2_s": float(sub["D0_std_mm2_s"].mean()),
                "selected_bstep": int(round(sub["selected_bstep"].mean())) if "selected_bstep" in sub.columns else np.nan,
                "selected_bvalue": float(sub["selected_bvalue"].mean()),
                "n_delta_app": int(round(sub["n_delta_app"].mean())),
                "delta_app_min_ms": float(sub["delta_app_min_ms"].min()),
                "delta_app_max_ms": float(sub["delta_app_max_ms"].max()),
                "reference_D0_mm2_s": float(sub["reference_D0_mm2_s"].iloc[0]),
                "reference_D0_error_mm2_s": float(sub["reference_D0_error_mm2_s"].iloc[0]),
                "model": "fixed_bvalue_mean",
            }
        )

    if not derived_rows:
        return raw
    derived = pd.DataFrame(derived_rows)
    return pd.concat([raw, derived], ignore_index=True, sort=False)


def compute_alpha_macro_summary(
    df_avg: pd.DataFrame,
    *,
    reference_D0: float = 0.0032,
    reference_D0_error: float = 0.0000283512,
    selected_bstep: int | None = None,
    roi_selected_bsteps: dict[str, int] | None = None,
    direction_aliases: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if reference_D0 <= 0:
        raise ValueError("reference_D0 must be > 0.")
    if df_avg.empty:
        raise ValueError("df_avg is empty.")

    rows: list[dict[str, object]] = []
    for (subj, roi, direction), sub in df_avg.groupby(["subj", "roi", "direction"], sort=True):
        bvalues = np.sort(sub["bvalue"].dropna().unique())
        if len(bvalues) == 0:
            continue
        roi_bstep = None
        if roi_selected_bsteps is not None:
            roi_bstep = roi_selected_bsteps.get(str(roi))
        selected_bvalue, chosen_bstep = _select_bvalue_and_bstep(
            bvalues,
            selected_bstep=roi_bstep if roi_bstep is not None else selected_bstep,
            group_label=f"subj={subj}, roi={roi}, direction={direction}",
        )
        chosen = sub[np.isclose(sub["bvalue"], selected_bvalue, atol=1e-9)].sort_values("Delta_app_ms")
        if chosen.empty:
            continue

        dvals = chosen["D_mean_mm2_s"].to_numpy(float)
        d0_mean = float(np.nanmean(dvals))
        d0_std = float(np.nanstd(dvals))
        alpha = float(d0_mean / reference_D0)

        rel_d = d0_std / d0_mean if d0_mean > 0 else np.nan
        rel_ref = float(reference_D0_error) / float(reference_D0) if reference_D0_error is not None else 0.0
        alpha_error = float(np.sqrt(rel_d**2 + rel_ref**2) * alpha) if np.isfinite(rel_d) else np.nan

        rows.append(
            {
                "subj": str(subj),
                "sheet": _merge_pipe_lists(chosen["sheets"].astype(str).tolist()),
                "roi": str(roi),
                "direction": str(direction),
                "alpha_macro": alpha,
                "alpha_macro_error": alpha_error,
                "D0_mean_mm2_s": d0_mean,
                "D0_std_mm2_s": d0_std,
                "selected_bstep": chosen_bstep,
                "selected_bvalue": selected_bvalue,
                "n_delta_app": int(chosen["Delta_app_ms"].nunique()),
                "delta_app_min_ms": float(chosen["Delta_app_ms"].min()),
                "delta_app_max_ms": float(chosen["Delta_app_ms"].max()),
                "reference_D0_mm2_s": float(reference_D0),
                "reference_D0_error_mm2_s": float(reference_D0_error),
                "model": "fixed_bvalue_mean",
            }
        )

    df_summary = pd.DataFrame(rows)
    if df_summary.empty:
        raise ValueError("Could not build alpha_macro: there were no valid groups.")

    aliases = direction_aliases or dict(DEFAULT_DIRECTION_ALIASES)
    df_summary = _expand_direction_alias_rows(df_summary, aliases)
    df_summary["region"] = df_summary["roi"]
    df_summary["direccion"] = df_summary["direction"]
    return df_avg.copy(), df_summary.sort_values(["subj", "roi", "direction"], kind="stable").reset_index(drop=True)


def plot_alpha_macro_vs_roi(
    df_summary: pd.DataFrame,
    *,
    out_png: str | Path,
    roi_order: Sequence[str],
    directions: Sequence[str],
    subjs: Sequence[str] | None = None,
    title_prefix: str = r"$\alpha_{macro}$",
) -> Path:
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    sub = df_summary[
        (df_summary["roi"].isin([str(x) for x in roi_order]))
        & (df_summary["direction"].isin([str(x) for x in directions]))
    ].copy()
    if "direction_kind" in sub.columns and (sub["direction_kind"] == "derived").any():
        sub = sub[sub["direction_kind"] == "derived"].copy()
    if sub.empty:
        raise ValueError(f"No data available for plotting {out_png.name}")

    subj_list = list(subjs) if subjs is not None else sorted(sub["subj"].astype(str).unique().tolist())
    fig, axes = plt.subplots(1, len(directions), figsize=(8 * len(directions), 6), sharey=True)
    if len(directions) == 1:
        axes = [axes]

    x = np.arange(len(roi_order))
    for ax, direction in zip(axes, directions):
        ddir = sub[sub["direction"] == str(direction)]
        for subj in subj_list:
            db = ddir[ddir["subj"] == str(subj)]
            yvals: list[float] = []
            evals: list[float] = []
            for roi in roi_order:
                row = db[db["roi"] == str(roi)]
                if row.empty:
                    yvals.append(np.nan)
                    evals.append(np.nan)
                    continue
                yvals.append(float(row["alpha_macro"].iloc[0]))
                evals.append(float(row["alpha_macro_error"].iloc[0]))

            y = np.asarray(yvals, dtype=float)
            yerr = np.asarray(evals, dtype=float)
            color = DEFAULT_SUBJ_COLORS.get(str(subj), None)
            ax.plot(x, y, "o-", linewidth=2, markersize=7, color=color, label=str(subj))
            if np.isfinite(yerr).any():
                low = y - np.nan_to_num(yerr, nan=0.0)
                high = y + np.nan_to_num(yerr, nan=0.0)
                ax.fill_between(x, low, high, alpha=0.2, color=color)

        ax.set_title(f"{title_prefix} | dir={direction}", fontsize=15)
        ax.set_xticks(x)
        ax.set_xticklabels(list(roi_order), rotation=45, ha="right")
        ax.set_xlabel("ROI", fontsize=12)
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.legend()

    axes[0].set_ylabel(r"$\alpha_{macro}$", fontsize=12)
    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    return out_png


def write_alpha_macro_outputs(
    df_avg: pd.DataFrame,
    df_summary: pd.DataFrame,
    *,
    out_summary_xlsx: Path,
    out_avg_xlsx: Path | None = None,
) -> None:
    write_xlsx_csv_outputs(df_summary, out_summary_xlsx, csv_path=out_summary_xlsx.with_suffix(".csv"))
    if out_avg_xlsx is not None:
        write_xlsx_csv_outputs(df_avg, out_avg_xlsx, csv_path=out_avg_xlsx.with_suffix(".csv"))


def ensure_list_or_none(values: Iterable[str] | None) -> list[str] | None:
    if values is None:
        return None
    out = [str(x) for x in values]
    return out or None
