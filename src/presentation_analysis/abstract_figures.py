"""Build ESMRMB-style comparison figures from signal-analysis master tables.

The functions in this module intentionally keep discovery, fitting, and plotting
separate.  A notebook can therefore select any available extraction variant
without duplicating the scientific calculations in individual cells.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.optimize import least_squares

from models.model_fitting import M_ogse_rest_offset, OGSE_contrast_vs_g_rest


CC_ROIS = ("AntCC", "MidAntCC", "CentralCC", "MidPostCC", "PostCC")
DEFAULT_DIRECTIONS = ("long", "tra")
GAMMA_RAD_MS_MT = 267.5221900

_KEY_COLUMNS = (
    "row_kind",
    "subj",
    "sheet",
    "roi",
    "direction",
    "stat",
    "source_file",
    "b_step",
)


def _ordered_unique(values: Iterable[object]) -> list[str]:
    return list(dict.fromkeys(str(value) for value in values))


def _safe_tag(value: object) -> str:
    return (
        str(value)
        .strip()
        .replace(" ", "_")
        .replace("/", "-")
        .replace("\\", "-")
        .replace(":", "-")
    )


def _require_columns(table: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise KeyError(f"{label} is missing required columns: {missing}")


def discover_analysis_runs(
    project_root: str | Path,
    *,
    dataset: str = "brains",
    sequence: str = "ogse",
    dwi_level: str | None = None,
) -> pd.DataFrame:
    """Discover master tables and parse their extraction tags.

    Expected layout::

        Data-BIDS/derivatives/signal_analysis/
          <dwi-level>--<variant>/<dataset>/<sequence>_experiments/master.long.parquet
    """

    project_root = Path(project_root).resolve()
    analysis_root = project_root / "Data-BIDS" / "derivatives" / "signal_analysis"
    rows: list[dict[str, object]] = []
    if not analysis_root.exists():
        raise FileNotFoundError(f"Signal-analysis root does not exist: {analysis_root}")

    pattern = f"*/{dataset}/{sequence}_experiments/master.long.parquet"
    for master_path in sorted(analysis_root.glob(pattern)):
        analysis_tag = master_path.parents[2].name
        if "--" in analysis_tag:
            parsed_dwi, variant = analysis_tag.split("--", 1)
        else:
            parsed_dwi, variant = analysis_tag, "default"
        if dwi_level is not None and parsed_dwi != str(dwi_level):
            continue
        experiment_root = master_path.parent
        rows.append(
            {
                "analysis_tag": analysis_tag,
                "dwi_level": parsed_dwi,
                "variant": variant,
                "dataset": dataset,
                "sequence": sequence,
                "experiment_root": experiment_root,
                "master_path": master_path,
                "alpha_path": experiment_root
                / "alpha_macro"
                / "master"
                / "summary_alpha_values.xlsx",
            }
        )

    runs = pd.DataFrame(rows)
    if runs.empty:
        scope = f"dataset={dataset!r}, sequence={sequence!r}, dwi_level={dwi_level!r}"
        raise FileNotFoundError(f"No signal-analysis masters were found for {scope} under {analysis_root}")
    return runs.sort_values(["dwi_level", "variant"]).reset_index(drop=True)


def _select_runs(runs: pd.DataFrame, variants: Sequence[str] | None) -> pd.DataFrame:
    if variants is None:
        return runs.copy()
    requested = [str(variant) for variant in variants]
    selected = runs[runs["variant"].isin(requested)].copy()
    missing = [variant for variant in requested if variant not in set(selected["variant"])]
    if missing:
        available = sorted(runs["variant"].astype(str).unique())
        raise ValueError(f"Unknown variants {missing}. Available variants: {available}")
    rank = {variant: index for index, variant in enumerate(requested)}
    return selected.sort_values("variant", key=lambda col: col.map(rank)).reset_index(drop=True)


def load_masters(runs: pd.DataFrame, variants: Sequence[str] | None = None) -> pd.DataFrame:
    """Load selected master tables and attach their variant provenance."""

    selected = _select_runs(runs, variants)
    frames: list[pd.DataFrame] = []
    for row in selected.itertuples(index=False):
        frame = pd.read_parquet(row.master_path)
        frame.insert(0, "variant", str(row.variant))
        frame.insert(1, "analysis_tag", str(row.analysis_tag))
        frame.insert(2, "dataset", str(row.dataset))
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def load_alpha_summaries(runs: pd.DataFrame, variants: Sequence[str] | None = None) -> pd.DataFrame:
    """Load the alpha-macro summary associated with each selected master."""

    selected = _select_runs(runs, variants)
    frames: list[pd.DataFrame] = []
    for row in selected.itertuples(index=False):
        alpha_path = Path(row.alpha_path)
        if not alpha_path.exists():
            raise FileNotFoundError(f"Alpha-macro summary does not exist: {alpha_path}")
        frame = pd.read_excel(alpha_path)
        rename = {}
        if "region" in frame.columns and "roi" not in frame.columns:
            rename["region"] = "roi"
        if "direccion" in frame.columns and "direction" not in frame.columns:
            rename["direccion"] = "direction"
        if "alpha" in frame.columns and "alpha_macro" not in frame.columns:
            rename["alpha"] = "alpha_macro"
        if "alpha_error" in frame.columns and "alpha_macro_error" not in frame.columns:
            rename["alpha_error"] = "alpha_macro_error"
        frame = frame.rename(columns=rename)
        _require_columns(frame, ["subj", "roi", "direction", "alpha_macro"], str(alpha_path))
        frame.insert(0, "variant", str(row.variant))
        frame.insert(1, "dataset", str(row.dataset))
        frame["roi"] = frame["roi"].astype(str).str.replace("_norm", "", regex=False)
        frame["direction"] = frame["direction"].astype(str).str.strip()
        frame["alpha_macro"] = pd.to_numeric(frame["alpha_macro"], errors="coerce")
        if "alpha_macro_error" not in frame.columns:
            frame["alpha_macro_error"] = np.nan
        frame["alpha_macro_error"] = pd.to_numeric(frame["alpha_macro_error"], errors="coerce")
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def master_qc_summary(masters: pd.DataFrame) -> pd.DataFrame:
    """Return one compact quality-control row per variant."""

    _require_columns(masters, ["variant", *_KEY_COLUMNS], "masters")
    rows: list[dict[str, object]] = []
    for variant, frame in masters.groupby("variant", sort=False):
        counts = frame["row_kind"].value_counts()
        rows.append(
            {
                "variant": variant,
                "rows": len(frame),
                "signal_rows": int(counts.get("signal", 0)),
                "rotated_rows": int(counts.get("signal_rotated", 0)),
                "subjects": frame["subj"].nunique(),
                "sheets": frame["sheet"].nunique(),
                "sources": frame["source_file"].nunique(),
                "exact_duplicates": int(frame.duplicated().sum()),
                "logical_duplicates": int(frame.duplicated(list(_KEY_COLUMNS)).sum()),
            }
        )
    return pd.DataFrame(rows)


def compare_with_reference(
    masters: pd.DataFrame,
    *,
    reference_variant: str,
    row_kind: str = "signal_rotated",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare logical row coverage with a reference variant.

    Returns a summary table and the missing logical keys for every variant.
    """

    frame = masters[masters["row_kind"].eq(row_kind)].copy()
    variants = _ordered_unique(frame["variant"])
    if reference_variant not in variants:
        raise ValueError(f"Reference variant {reference_variant!r} is not loaded")
    compare_keys = [column for column in _KEY_COLUMNS if column != "row_kind"]
    reference = frame[frame["variant"].eq(reference_variant)][compare_keys].drop_duplicates()
    summaries: list[dict[str, object]] = []
    missing_frames: list[pd.DataFrame] = []
    for variant in variants:
        candidate = frame[frame["variant"].eq(variant)][compare_keys].drop_duplicates()
        merged = reference.merge(candidate, on=compare_keys, how="left", indicator=True)
        missing = merged[merged["_merge"].eq("left_only")].drop(columns="_merge")
        missing.insert(0, "variant", variant)
        missing_frames.append(missing)
        reverse = candidate.merge(reference, on=compare_keys, how="left", indicator=True)
        summaries.append(
            {
                "variant": variant,
                "reference": reference_variant,
                "reference_rows": len(reference),
                "matched_rows": len(reference) - len(missing),
                "missing_rows": len(missing),
                "extra_rows": int(reverse["_merge"].eq("left_only").sum()),
                "coverage_percent": 100.0 * (len(reference) - len(missing)) / len(reference),
            }
        )
    return pd.DataFrame(summaries), pd.concat(missing_frames, ignore_index=True)


def _unique_by_pair(frame: pd.DataFrame, pair_keys: Sequence[str], label: str) -> pd.DataFrame:
    duplicated = frame.duplicated(list(pair_keys), keep=False)
    if duplicated.any():
        example = frame.loc[duplicated, list(pair_keys) + ["source_file"]].head(10)
        raise ValueError(f"{label} has ambiguous N-pair rows. Examples:\n{example.to_string(index=False)}")
    return frame


def build_contrast_table(
    masters: pd.DataFrame,
    *,
    n_high: int = 8,
    n_low: int = 4,
    directions: Sequence[str] = DEFAULT_DIRECTIONS,
    rois: Sequence[str] | None = CC_ROIS,
    value_column: str = "value_norm",
    gradient_column: str = "g_thorsten",
    apply_gradient_correction: bool = True,
) -> pd.DataFrame:
    """Build a direct b-step-paired contrast for diagnostic auditing only.

    This table has two different gradients per row and must not be used for
    peak-gradient or transformed-axis analyses. Use
    :func:`build_resampled_contrasts` for the ESMRMB analysis.
    """

    required = [
        "variant",
        "row_kind",
        "subj",
        "sheet",
        "roi",
        "direction",
        "stat",
        "td_ms",
        "b_step",
        "N",
        "source_file",
        value_column,
        gradient_column,
        "grad_correction_factor",
    ]
    _require_columns(masters, required, "masters")
    frame = masters[
        masters["row_kind"].eq("signal_rotated")
        & masters["direction"].astype(str).isin([str(value) for value in directions])
        & pd.to_numeric(masters["N"], errors="coerce").isin([int(n_high), int(n_low)])
        & masters["stat"].astype(str).isin(["avg", "std"])
    ].copy()
    if rois is not None:
        frame = frame[frame["roi"].astype(str).isin([str(roi) for roi in rois])].copy()
    if frame.empty:
        raise ValueError("No rotated rows match the requested contrast configuration")

    pair_keys = ["variant", "subj", "sheet", "roi", "direction", "td_ms", "b_step", "stat"]
    keep = pair_keys + [
        "N",
        "Hz",
        "sequence",
        "source_file",
        value_column,
        gradient_column,
        "grad_correction_factor",
        "delta_ms",
        "Delta_app_ms",
    ]
    high = frame[pd.to_numeric(frame["N"], errors="coerce").eq(int(n_high))][keep]
    low = frame[pd.to_numeric(frame["N"], errors="coerce").eq(int(n_low))][keep]
    high = _unique_by_pair(high, pair_keys, f"N={n_high}")
    low = _unique_by_pair(low, pair_keys, f"N={n_low}")
    paired = high.merge(low, on=pair_keys, how="inner", suffixes=("_high", "_low"), validate="one_to_one")
    if paired.empty:
        raise ValueError(f"No pointwise N={n_high} versus N={n_low} pairs were found")

    mean_rows = paired[paired["stat"].eq("avg")].copy()
    mean_rows["contrast"] = (
        pd.to_numeric(mean_rows[f"{value_column}_high"], errors="coerce")
        - pd.to_numeric(mean_rows[f"{value_column}_low"], errors="coerce")
    )

    std_rows = paired[paired["stat"].eq("std")].copy()
    std_keys = [key for key in pair_keys if key != "stat"]
    if not std_rows.empty:
        std_rows["contrast_std"] = np.sqrt(
            pd.to_numeric(std_rows[f"{value_column}_high"], errors="coerce") ** 2
            + pd.to_numeric(std_rows[f"{value_column}_low"], errors="coerce") ** 2
        )
        mean_rows = mean_rows.merge(std_rows[std_keys + ["contrast_std"]], on=std_keys, how="left")
    else:
        mean_rows["contrast_std"] = np.nan

    mean_rows["N_high"] = int(n_high)
    mean_rows["N_low"] = int(n_low)
    mean_rows["g_high_raw"] = pd.to_numeric(mean_rows[f"{gradient_column}_high"], errors="coerce")
    mean_rows["g_low_raw"] = pd.to_numeric(mean_rows[f"{gradient_column}_low"], errors="coerce")
    if apply_gradient_correction:
        corr_high = pd.to_numeric(mean_rows["grad_correction_factor_high"], errors="coerce")
        corr_low = pd.to_numeric(mean_rows["grad_correction_factor_low"], errors="coerce")
        if corr_high.isna().any() or corr_low.isna().any():
            raise ValueError("Gradient correction was requested, but one or more N-pair rows lack a factor")
    else:
        corr_high = pd.Series(1.0, index=mean_rows.index)
        corr_low = pd.Series(1.0, index=mean_rows.index)
    mean_rows["correction_high"] = corr_high
    mean_rows["correction_low"] = corr_low
    mean_rows["g_high_corr"] = mean_rows["g_high_raw"] * corr_high
    mean_rows["g_low_corr"] = mean_rows["g_low_raw"] * corr_low
    return mean_rows.drop(columns="stat").sort_values(std_keys).reset_index(drop=True)


def _r2(y: np.ndarray, yhat: np.ndarray) -> float:
    residual = float(np.sum((y - yhat) ** 2))
    total = float(np.sum((y - np.mean(y)) ** 2))
    return float(1.0 - residual / total) if total > 0 else np.nan


def build_resampled_contrasts(
    masters: pd.DataFrame,
    *,
    d0_m2_ms: float,
    n_high: int = 8,
    n_low: int = 4,
    directions: Sequence[str] = DEFAULT_DIRECTIONS,
    rois: Sequence[str] | None = CC_ROIS,
    value_column: str = "value_norm",
    gradient_column: str = "g_thorsten",
    apply_gradient_correction: bool = True,
    tc_bounds_ms: tuple[float, float] = (0.1, 10000.0),
    rician_c_bounds: tuple[float, float] = (0.0, 2.0),
    grid_size: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit paired OGSE signals and resample their difference on one gradient grid.

    The two branches are not subtracted by acquisition index. Each branch has
    its own restricted-model correlation time, while the normalized Rician
    floor ``C`` is shared by the N-high/N-low pair. Both fitted signals are
    then evaluated on the overlap of their corrected gradient ranges.
    """

    required = [
        "variant",
        "row_kind",
        "subj",
        "sheet",
        "roi",
        "direction",
        "stat",
        "td_ms",
        "b_step",
        "N",
        "source_file",
        value_column,
        gradient_column,
        "grad_correction_factor",
    ]
    _require_columns(masters, required, "masters")
    frame = masters[
        masters["row_kind"].eq("signal_rotated")
        & masters["direction"].astype(str).isin([str(value) for value in directions])
        & pd.to_numeric(masters["N"], errors="coerce").isin([int(n_high), int(n_low)])
        & masters["stat"].astype(str).eq("avg")
    ].copy()
    if rois is not None:
        frame = frame[frame["roi"].astype(str).isin([str(roi) for roi in rois])].copy()
    if frame.empty:
        raise ValueError("No rotated mean-signal rows match the requested resampled-contrast configuration")
    if int(grid_size) < 32:
        raise ValueError("grid_size must be at least 32")

    group_columns = ["variant", "subj", "sheet", "roi", "direction", "td_ms"]
    fit_rows: list[dict[str, object]] = []
    contrast_frames: list[pd.DataFrame] = []
    for key, group in frame.groupby(group_columns, sort=False, dropna=False):
        metadata = dict(zip(group_columns, key))
        high = group[pd.to_numeric(group["N"], errors="coerce").eq(int(n_high))].copy()
        low = group[pd.to_numeric(group["N"], errors="coerce").eq(int(n_low))].copy()
        fit_row: dict[str, object] = {
            **metadata,
            "N_high": int(n_high),
            "N_low": int(n_low),
            "D0_m2_ms": float(d0_m2_ms),
            "M0": 1.0,
            "n_points_high": 0,
            "n_points_low": 0,
            "ok": False,
            "message": "",
        }
        if high.empty or low.empty:
            fit_row["message"] = "Missing N-high or N-low signal branch"
            fit_rows.append(fit_row)
            continue

        def prepare_branch(branch: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            g_raw = pd.to_numeric(branch[gradient_column], errors="coerce").to_numpy(float)
            signal = pd.to_numeric(branch[value_column], errors="coerce").to_numpy(float)
            if apply_gradient_correction:
                correction = pd.to_numeric(branch["grad_correction_factor"], errors="coerce").to_numpy(float)
            else:
                correction = np.ones(len(branch), dtype=float)
            finite = np.isfinite(g_raw) & np.isfinite(signal) & np.isfinite(correction)
            order = np.argsort(g_raw[finite] * correction[finite])
            return g_raw[finite][order], (g_raw[finite] * correction[finite])[order], signal[finite][order]

        g_high_raw, g_high_corr, signal_high = prepare_branch(high)
        g_low_raw, g_low_corr, signal_low = prepare_branch(low)
        fit_row["n_points_high"] = int(len(signal_high))
        fit_row["n_points_low"] = int(len(signal_low))
        if len(signal_high) < 3 or len(signal_low) < 3:
            fit_row["message"] = "Fewer than three finite points in one signal branch"
            fit_rows.append(fit_row)
            continue
        if np.nanmax(g_high_corr) <= 0 or np.nanmax(g_low_corr) <= 0:
            fit_row["message"] = "A signal branch has no positive corrected gradient"
            fit_rows.append(fit_row)
            continue

        td_ms = float(metadata["td_ms"])

        def predict(g_corr: np.ndarray, n_value: int, tc_ms: float, c_value: float) -> np.ndarray:
            return np.asarray(
                M_ogse_rest_offset(
                    td_ms,
                    g_corr,
                    int(n_value),
                    td_ms / int(n_value),
                    tc_ms,
                    1.0,
                    float(d0_m2_ms),
                    c_value,
                ),
                dtype=float,
            )

        def residual(parameters: np.ndarray) -> np.ndarray:
            tc_high_ms = float(np.exp(parameters[0]))
            tc_low_ms = float(np.exp(parameters[1]))
            c_value = float(parameters[2])
            return np.concatenate(
                [
                    predict(g_high_corr, n_high, tc_high_ms, c_value) - signal_high,
                    predict(g_low_corr, n_low, tc_low_ms, c_value) - signal_low,
                ]
            )

        observed_floor = float(np.nanmedian(np.r_[signal_high[-2:], signal_low[-2:]]))
        observed_floor = float(np.clip(observed_floor, 0.0, 0.89))
        c_initial = observed_floor / np.sqrt(max(1e-12, 1.0 - observed_floor**2))
        c_initial = float(np.clip(c_initial, rician_c_bounds[0] + 1e-8, rician_c_bounds[1] - 1e-8))
        lower = np.array([np.log(tc_bounds_ms[0]), np.log(tc_bounds_ms[0]), rician_c_bounds[0]])
        upper = np.array([np.log(tc_bounds_ms[1]), np.log(tc_bounds_ms[1]), rician_c_bounds[1]])
        try:
            result = least_squares(
                residual,
                x0=np.array([np.log(5.0), np.log(5.0), c_initial]),
                bounds=(lower, upper),
                max_nfev=200000,
            )
            tc_high_ms = float(np.exp(result.x[0]))
            tc_low_ms = float(np.exp(result.x[1]))
            c_value = float(result.x[2])
            errors = _standard_errors(result, ("exp", "exp", "identity"))
            predicted_high = predict(g_high_corr, n_high, tc_high_ms, c_value)
            predicted_low = predict(g_low_corr, n_low, tc_low_ms, c_value)
            observed = np.concatenate([signal_high, signal_low])
            predicted = np.concatenate([predicted_high, predicted_low])

            # Restrict the common grid to the gradient interval supported by
            # both acquisitions; this avoids extrapolation-driven peaks.
            g_max_corr = float(min(np.nanmax(g_high_corr), np.nanmax(g_low_corr)))
            g_resampled = np.linspace(0.0, g_max_corr, int(grid_size))
            signal_high_fit = predict(g_resampled, n_high, tc_high_ms, c_value)
            signal_low_fit = predict(g_resampled, n_low, tc_low_ms, c_value)
            contrast = signal_high_fit - signal_low_fit
            peak_index = int(np.nanargmax(contrast))
            g_peak_corr = float(g_resampled[peak_index])
            lcf_peak_um, tc_peak_ms = _derived_peak_axes(td_ms, g_peak_corr, float(d0_m2_ms))

            fit_row.update(
                {
                    "ok": bool(result.success),
                    "message": str(result.message),
                    "tc_high_ms": tc_high_ms,
                    "tc_high_error_ms": float(errors[0]),
                    "tc_low_ms": tc_low_ms,
                    "tc_low_error_ms": float(errors[1]),
                    "C": c_value,
                    "C_error": float(errors[2]),
                    "rician_floor": c_value / np.sqrt(1.0 + c_value**2),
                    "r2_high": _r2(signal_high, predicted_high),
                    "r2_low": _r2(signal_low, predicted_low),
                    "r2": _r2(observed, predicted),
                    "rmse": float(np.sqrt(np.mean((observed - predicted) ** 2))),
                    "g_max_corr_mTm": g_max_corr,
                    "g_peak_corr_mTm": g_peak_corr,
                    "signal_peak": float(contrast[peak_index]),
                    "lcf_peak_um": lcf_peak_um,
                    "tc_peak_ms": tc_peak_ms,
                    "source_file_high": str(high["source_file"].iloc[0]),
                    "source_file_low": str(low["source_file"].iloc[0]),
                }
            )
            lcf_axes = [
                _derived_peak_axes(td_ms, float(g_value), float(d0_m2_ms))
                for g_value in g_resampled
            ]
            contrast_frames.append(
                pd.DataFrame(
                    {
                        **{column: value for column, value in metadata.items()},
                        "grid_index": np.arange(len(g_resampled), dtype=int),
                        "N_high": int(n_high),
                        "N_low": int(n_low),
                        "g_resampled_corr": g_resampled,
                        "g_high_corr": g_resampled,
                        "g_low_corr": g_resampled,
                        "signal_high_fit": signal_high_fit,
                        "signal_low_fit": signal_low_fit,
                        "contrast": contrast,
                        "contrast_std": np.nan,
                        "lcf_um": [value[0] for value in lcf_axes],
                        "tc_axis_ms": [value[1] for value in lcf_axes],
                    }
                )
            )
        except Exception as error:  # preserve failed groups in the audit table
            fit_row["message"] = f"{type(error).__name__}: {error}"
        fit_rows.append(fit_row)

    fit_summary = pd.DataFrame(fit_rows)
    if not contrast_frames:
        raise ValueError("No resampled contrasts could be constructed")
    resampled = pd.concat(contrast_frames, ignore_index=True)
    return resampled, fit_summary


def _standard_errors(result, transform: Sequence[str]) -> np.ndarray:
    errors = np.full(len(transform), np.nan, dtype=float)
    jacobian = np.asarray(result.jac, dtype=float)
    dof = jacobian.shape[0] - jacobian.shape[1]
    if dof <= 0:
        return errors
    try:
        covariance = np.linalg.pinv(jacobian.T @ jacobian) * (2.0 * result.cost / dof)
        errors = np.sqrt(np.clip(np.diag(covariance), 0.0, np.inf))
    except np.linalg.LinAlgError:
        return errors
    for index, kind in enumerate(transform):
        if kind == "exp":
            errors[index] *= float(np.exp(result.x[index]))
    return errors


def _derived_peak_axes(td_ms: float, g_raw_mtm: float, d0_m2_ms: float) -> tuple[float, float]:
    if not np.isfinite(g_raw_mtm) or g_raw_mtm <= 0:
        return np.nan, np.nan
    length_g = ((2.0 ** 1.5) * d0_m2_ms / (GAMMA_RAD_MS_MT * g_raw_mtm)) ** (1.0 / 3.0)
    length_d = np.sqrt(2.0 * d0_m2_ms * td_ms)
    dimensionless_d = length_d / length_g
    dimensionless_cf = ((3.0 / 2.0) ** 0.25) * dimensionless_d ** (-0.5)
    length_cf_m = dimensionless_cf * length_g
    tc_peak_ms = length_cf_m**2 / (2.0 * d0_m2_ms)
    return float(length_cf_m * 1e6), float(tc_peak_ms)


def fit_rest_contrasts(
    contrast: pd.DataFrame,
    *,
    d0_m2_ms: float,
    tc_bounds_ms: tuple[float, float] = (0.1, 1000.0),
    m0_bounds: tuple[float, float] = (0.0, 5.0),
    peak_grid_size: int = 1000,
) -> pd.DataFrame:
    """Fit a two-gradient diagnostic contrast; not used by the ESMRMB notebook."""

    group_columns = ["variant", "subj", "sheet", "roi", "direction", "td_ms"]
    _require_columns(
        contrast,
        group_columns + ["contrast", "g_high_corr", "g_low_corr", "g_high_raw", "correction_high"],
        "contrast",
    )
    rows: list[dict[str, object]] = []
    for key, group in contrast.groupby(group_columns, sort=False, dropna=False):
        metadata = dict(zip(group_columns, key))
        group = group.sort_values("b_step")
        y = pd.to_numeric(group["contrast"], errors="coerce").to_numpy(float)
        g_high = pd.to_numeric(group["g_high_corr"], errors="coerce").to_numpy(float)
        g_low = pd.to_numeric(group["g_low_corr"], errors="coerce").to_numpy(float)
        finite = np.isfinite(y) & np.isfinite(g_high) & np.isfinite(g_low)
        y, g_high, g_low = y[finite], g_high[finite], g_low[finite]
        n_high = int(group["N_high"].iloc[0])
        n_low = int(group["N_low"].iloc[0])
        td_ms = float(metadata["td_ms"])
        fit_row: dict[str, object] = {
            **metadata,
            "N_high": n_high,
            "N_low": n_low,
            "D0_m2_ms": float(d0_m2_ms),
            "n_points": int(len(y)),
            "ok": False,
            "message": "",
        }
        if len(y) < 3:
            fit_row["message"] = "Fewer than three finite contrast points"
            rows.append(fit_row)
            continue

        def residual(parameters: np.ndarray) -> np.ndarray:
            tc_ms = float(np.exp(parameters[0]))
            m0 = float(parameters[1])
            prediction = OGSE_contrast_vs_g_rest(td_ms, g_high, g_low, n_high, n_low, tc_ms, m0, d0_m2_ms)
            return np.asarray(prediction, dtype=float) - y

        try:
            result = least_squares(
                residual,
                x0=np.array([np.log(5.0), 1.0]),
                bounds=(
                    np.array([np.log(tc_bounds_ms[0]), m0_bounds[0]]),
                    np.array([np.log(tc_bounds_ms[1]), m0_bounds[1]]),
                ),
                max_nfev=200000,
            )
            tc_ms = float(np.exp(result.x[0]))
            m0 = float(result.x[1])
            errors = _standard_errors(result, ("exp", "identity"))
            prediction = np.asarray(
                OGSE_contrast_vs_g_rest(td_ms, g_high, g_low, n_high, n_low, tc_ms, m0, d0_m2_ms),
                dtype=float,
            )
            fraction = np.linspace(0.0, 1.0, max(64, int(peak_grid_size)))
            g_high_max = float(np.nanmax(g_high))
            g_low_max = float(np.nanmax(g_low))
            peak_curve = np.asarray(
                OGSE_contrast_vs_g_rest(
                    td_ms,
                    fraction * g_high_max,
                    fraction * g_low_max,
                    n_high,
                    n_low,
                    tc_ms,
                    m0,
                    d0_m2_ms,
                ),
                dtype=float,
            )
            peak_index = int(np.nanargmax(peak_curve))
            peak_fraction = float(fraction[peak_index])
            correction_high = float(pd.to_numeric(group["correction_high"], errors="coerce").median())
            g_peak_corr = peak_fraction * g_high_max
            g_peak_raw = g_peak_corr / correction_high
            lcf_peak_um, tc_peak_ms = _derived_peak_axes(td_ms, g_peak_raw, d0_m2_ms)
            fit_row.update(
                {
                    "ok": bool(result.success),
                    "message": str(result.message),
                    "tc_fit_ms": tc_ms,
                    "tc_fit_error_ms": float(errors[0]),
                    "M0": m0,
                    "M0_error": float(errors[1]),
                    "rmse": float(np.sqrt(np.mean((y - prediction) ** 2))),
                    "r2": _r2(y, prediction),
                    "peak_fraction": peak_fraction,
                    "signal_peak": float(peak_curve[peak_index]),
                    "g_high_max_corr_mTm": g_high_max,
                    "g_low_max_corr_mTm": g_low_max,
                    "correction_high": correction_high,
                    "g_peak_corr_mTm": g_peak_corr,
                    "g_peak_raw_mTm": g_peak_raw,
                    "lcf_peak_um": lcf_peak_um,
                    "tc_peak_ms": tc_peak_ms,
                }
            )
        except Exception as error:  # keep a complete audit table for failed fits
            fit_row["message"] = f"{type(error).__name__}: {error}"
        rows.append(fit_row)
    return pd.DataFrame(rows)


def model_curve_for_fit(fit_row: Mapping[str, object], *, grid_size: int = 500) -> pd.DataFrame:
    """Evaluate one fitted contrast curve and return raw gradient and lcf axes."""

    fraction = np.linspace(0.0, 1.0, max(32, int(grid_size)))
    g_high_corr = fraction * float(fit_row["g_high_max_corr_mTm"])
    g_low_corr = fraction * float(fit_row["g_low_max_corr_mTm"])
    values = OGSE_contrast_vs_g_rest(
        float(fit_row["td_ms"]),
        g_high_corr,
        g_low_corr,
        int(fit_row["N_high"]),
        int(fit_row["N_low"]),
        float(fit_row["tc_fit_ms"]),
        float(fit_row["M0"]),
        float(fit_row["D0_m2_ms"]),
    )
    g_high_raw = g_high_corr / float(fit_row["correction_high"])
    axes = [
        _derived_peak_axes(float(fit_row["td_ms"]), float(g_value), float(fit_row["D0_m2_ms"]))
        for g_value in g_high_raw
    ]
    return pd.DataFrame(
        {
            "fraction": fraction,
            "g_high_raw": g_high_raw,
            "g_high_corr": g_high_corr,
            "lcf_um": [value[0] for value in axes],
            "tc_axis_ms": [value[1] for value in axes],
            "contrast_fit": np.asarray(values, dtype=float),
        }
    )


def aggregate_alpha(alpha: pd.DataFrame) -> pd.DataFrame:
    """Collapse repeated sheets to one alpha estimate per subject, ROI, and direction."""

    keys = ["variant", "subj", "roi", "direction"]
    _require_columns(alpha, keys + ["alpha_macro", "alpha_macro_error"], "alpha")
    rows: list[dict[str, object]] = []
    for key, group in alpha.groupby(keys, sort=False, dropna=False):
        values = pd.to_numeric(group["alpha_macro"], errors="coerce").dropna().to_numpy(float)
        errors = pd.to_numeric(group["alpha_macro_error"], errors="coerce").dropna().to_numpy(float)
        if not len(values):
            continue
        rows.append(
            {
                **dict(zip(keys, key)),
                "alpha_macro": float(np.mean(values)),
                "alpha_macro_error": float(np.sqrt(np.sum(errors**2)) / len(errors)) if len(errors) else np.nan,
                "alpha_n_sheets": int(group["sheet"].nunique()) if "sheet" in group else 1,
            }
        )
    return pd.DataFrame(rows)


def _tc_pseudohuber(td_ms: np.ndarray, c: float, delta: float, alpha_macro: float) -> np.ndarray:
    td_ms = np.asarray(td_ms, dtype=float)
    return c + alpha_macro * delta * (np.sqrt(1.0 + (td_ms / delta) ** 2) - 1.0)


def fit_tc_pseudohuber(
    contrast_fits: pd.DataFrame,
    alpha_aggregated: pd.DataFrame,
    *,
    exclude_td_ms: Sequence[float] = (76.0,),
    td_range_ms: tuple[float, float] = (50.0, 250.0),
    c_bounds_ms: tuple[float, float] = (0.0, 10.0),
    delta_bounds_ms: tuple[float, float] = (1e-6, 10000.0),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit peak correlation time versus diffusion time with fixed alpha macro."""

    keys = ["variant", "subj", "roi", "direction"]
    valid = contrast_fits[contrast_fits["ok"].fillna(False)].copy()
    valid = valid[pd.to_numeric(valid["tc_peak_ms"], errors="coerce").notna()]
    valid = valid[
        pd.to_numeric(valid["td_ms"], errors="coerce").between(td_range_ms[0], td_range_ms[1], inclusive="both")
    ]
    for excluded in exclude_td_ms:
        valid = valid[~np.isclose(pd.to_numeric(valid["td_ms"], errors="coerce"), float(excluded), atol=1e-3)]

    tc_data = (
        valid.groupby(keys + ["td_ms"], as_index=False, dropna=False)
        .agg(tc_peak_ms=("tc_peak_ms", "mean"), tc_peak_std_ms=("tc_peak_ms", "std"), n_curves=("tc_peak_ms", "size"))
        .sort_values(keys + ["td_ms"])
    )
    alpha_lookup = alpha_aggregated.set_index(keys)
    rows: list[dict[str, object]] = []
    for key, group in tc_data.groupby(keys, sort=False, dropna=False):
        metadata = dict(zip(keys, key))
        fit_row: dict[str, object] = {**metadata, "ok": False, "message": ""}
        if key not in alpha_lookup.index:
            fit_row["message"] = "No matching alpha-macro estimate"
            rows.append(fit_row)
            continue
        alpha_row = alpha_lookup.loc[key]
        if isinstance(alpha_row, pd.DataFrame):
            alpha_row = alpha_row.iloc[0]
        alpha_macro = float(alpha_row["alpha_macro"])
        alpha_error = float(alpha_row.get("alpha_macro_error", np.nan))
        x = pd.to_numeric(group["td_ms"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(group["tc_peak_ms"], errors="coerce").to_numpy(float)
        finite = np.isfinite(x) & np.isfinite(y)
        x, y = x[finite], y[finite]
        fit_row.update({"alpha_macro": alpha_macro, "alpha_macro_error": alpha_error, "n_td": int(len(x))})
        if len(x) < 3:
            fit_row["message"] = "Fewer than three diffusion times"
            rows.append(fit_row)
            continue

        def residual(parameters: np.ndarray) -> np.ndarray:
            return _tc_pseudohuber(x, float(parameters[0]), float(parameters[1]), alpha_macro) - y

        try:
            result = least_squares(
                residual,
                x0=np.array(
                    [
                        np.clip(float(np.min(y)), c_bounds_ms[0], c_bounds_ms[1]),
                        np.clip(float(np.median(x)), delta_bounds_ms[0], delta_bounds_ms[1]),
                    ]
                ),
                bounds=(np.array([c_bounds_ms[0], delta_bounds_ms[0]]), np.array([c_bounds_ms[1], delta_bounds_ms[1]])),
                max_nfev=200000,
            )
            c, delta = (float(value) for value in result.x)
            errors = _standard_errors(result, ("identity", "identity"))
            prediction = _tc_pseudohuber(x, c, delta, alpha_macro)
            fit_row.update(
                {
                    "ok": bool(result.success),
                    "message": str(result.message),
                    "c_ms": c,
                    "c_error_ms": float(errors[0]),
                    "delta_ms": delta,
                    "delta_error_ms": float(errors[1]),
                    "r2": _r2(y, prediction),
                    "rmse_ms": float(np.sqrt(np.mean((y - prediction) ** 2))),
                }
            )
        except Exception as error:
            fit_row["message"] = f"{type(error).__name__}: {error}"
        rows.append(fit_row)
    return tc_data, pd.DataFrame(rows)


def merge_alpha_delta(alpha_aggregated: pd.DataFrame, tc_summary: pd.DataFrame) -> pd.DataFrame:
    """Create one subject-level table used by the Figure 3/4 analogues."""

    keys = ["variant", "subj", "roi", "direction"]
    columns = keys + ["delta_ms", "delta_error_ms", "c_ms", "r2", "n_td", "ok"]
    available = [column for column in columns if column in tc_summary.columns]
    return alpha_aggregated.merge(tc_summary[available], on=keys, how="inner", validate="one_to_one")


def _save_figure(fig: plt.Figure, output_path: str | Path | None, *, dpi: int = 300) -> None:
    if output_path is None:
        return
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=int(dpi), bbox_inches="tight")


def _resolve_example(
    table: pd.DataFrame,
    *,
    subject: str | None,
    roi: str | None,
    direction: str | None,
    td_ms: float | None,
) -> tuple[str, str, str, float]:
    subset = table.copy()
    if subject is not None:
        subset = subset[subset["subj"].eq(subject)]
    if roi is not None:
        subset = subset[subset["roi"].eq(roi)]
    if direction is not None:
        subset = subset[subset["direction"].eq(direction)]
    if subset.empty:
        raise ValueError(
            f"No example data exist for subject={subject!r}, roi={roi!r}, direction={direction!r}"
        )

    # Prefer a subject/ROI/direction combination represented in the largest
    # number of variants. This keeps auto-selected examples comparable when a
    # segmentation removes an ROI completely in one acquisition.
    combination_keys = ["subj", "roi", "direction"]
    if "variant" in subset.columns:
        combinations = (
            subset.groupby(combination_keys, as_index=False, dropna=False)["variant"]
            .nunique()
            .rename(columns={"variant": "variant_count"})
            .sort_values(["variant_count", *combination_keys], ascending=[False, True, True, True])
        )
    else:
        combinations = subset[combination_keys].drop_duplicates().sort_values(combination_keys)
    selected = combinations.iloc[0]
    resolved_subject = str(selected["subj"])
    resolved_roi = str(selected["roi"])
    resolved_direction = str(selected["direction"])
    subset = subset[
        subset["subj"].eq(selected["subj"])
        & subset["roi"].eq(selected["roi"])
        & subset["direction"].eq(selected["direction"])
    ]

    td_values = pd.to_numeric(subset["td_ms"], errors="coerce")
    if "variant" in subset.columns:
        td_candidates = (
            subset.assign(_td_ms=td_values)
            .dropna(subset=["_td_ms"])
            .groupby("_td_ms", as_index=False)["variant"]
            .nunique()
            .rename(columns={"variant": "variant_count"})
        )
        maximum_coverage = td_candidates["variant_count"].max()
        available_td = sorted(td_candidates.loc[td_candidates["variant_count"].eq(maximum_coverage), "_td_ms"])
    else:
        available_td = sorted(td_values.dropna().unique())
    if td_ms is None:
        resolved_td = float(available_td[len(available_td) // 2])
    else:
        resolved_td = float(min(available_td, key=lambda value: abs(float(value) - float(td_ms))))
    return resolved_subject, resolved_roi, resolved_direction, resolved_td


def plot_signal_contrast_example(
    masters: pd.DataFrame,
    contrast: pd.DataFrame,
    contrast_fits: pd.DataFrame,
    *,
    variants: Sequence[str],
    subject: str | None,
    sheet: str | None = None,
    roi: str | None = None,
    direction: str | None = None,
    td_ms: float | None = 143.4,
    n_high: int = 8,
    n_low: int = 4,
    output_path: str | Path | None = None,
    dpi: int = 300,
) -> plt.Figure:
    """Plot the Figure 1b analogue for every selected extraction variant."""

    resolver_source = contrast if not contrast.empty else contrast_fits
    resolver_table = (
        resolver_source
        if sheet is None
        else resolver_source[resolver_source["sheet"].astype(str).eq(str(sheet))]
    )
    subject, roi, direction, td_ms = _resolve_example(
        resolver_table,
        subject=subject,
        roi=roi,
        direction=direction,
        td_ms=td_ms,
    )
    fig, axes = plt.subplots(len(variants), 2, figsize=(11, 4.1 * len(variants)), squeeze=False)
    sheet_mask_master = pd.Series(True, index=masters.index) if sheet is None else masters["sheet"].astype(str).eq(str(sheet))
    sheet_mask_contrast = (
        pd.Series(True, index=contrast.index) if sheet is None else contrast["sheet"].astype(str).eq(str(sheet))
    )
    sheet_mask_fit = (
        pd.Series(True, index=contrast_fits.index)
        if sheet is None
        else contrast_fits["sheet"].astype(str).eq(str(sheet))
    )
    for row_index, variant in enumerate(variants):
        signal_ax, contrast_ax = axes[row_index]
        signal = masters[
            masters["variant"].eq(variant)
            & masters["row_kind"].eq("signal_rotated")
            & masters["subj"].eq(subject)
            & masters["roi"].eq(roi)
            & masters["direction"].eq(direction)
            & masters["stat"].eq("avg")
            & sheet_mask_master
            & np.isclose(pd.to_numeric(masters["td_ms"], errors="coerce"), td_ms)
            & pd.to_numeric(masters["N"], errors="coerce").isin([n_high, n_low])
        ].copy()
        fit = contrast_fits[
            contrast_fits["variant"].eq(variant)
            & contrast_fits["subj"].eq(subject)
            & contrast_fits["roi"].eq(roi)
            & contrast_fits["direction"].eq(direction)
            & sheet_mask_fit
            & np.isclose(pd.to_numeric(contrast_fits["td_ms"], errors="coerce"), td_ms)
            & contrast_fits["ok"].fillna(False)
        ]
        for n_value, color in [(n_high, "#1f77b4"), (n_low, "#ff7f0e")]:
            series = signal[pd.to_numeric(signal["N"], errors="coerce").eq(n_value)].sort_values("g_thorsten")
            x = pd.to_numeric(series["g_thorsten"], errors="coerce") * pd.to_numeric(
                series["grad_correction_factor"], errors="coerce"
            )
            signal_ax.plot(x, series["value_norm"], "o", color=color, label=f"N={n_value} data")
            if not fit.empty and len(x):
                fit_row = fit.iloc[0]
                tc_column = "tc_high_ms" if n_value == n_high else "tc_low_ms"
                fit_grid = np.linspace(0.0, float(np.nanmax(x)), 500)
                fit_values = M_ogse_rest_offset(
                    td_ms,
                    fit_grid,
                    n_value,
                    td_ms / n_value,
                    float(fit_row[tc_column]),
                    1.0,
                    float(fit_row["D0_m2_ms"]),
                    float(fit_row["C"]),
                )
                signal_ax.plot(fit_grid, fit_values, "-", color=color, label=f"N={n_value} fit")
        points = contrast[
            contrast["variant"].eq(variant)
            & contrast["subj"].eq(subject)
            & contrast["roi"].eq(roi)
            & contrast["direction"].eq(direction)
            & sheet_mask_contrast
            & np.isclose(pd.to_numeric(contrast["td_ms"], errors="coerce"), td_ms)
        ].sort_values("g_resampled_corr")
        signal_ax.set_ylabel("Normalized signal")
        signal_ax.set_xlabel("Corrected gradient [mT/m]")
        signal_ax.set_title(f"{variant}: shared-C signal fits")
        signal_ax.grid(alpha=0.25)
        if signal.empty:
            signal_ax.text(0.5, 0.5, "No matching signal pair", ha="center", va="center", transform=signal_ax.transAxes)
        else:
            signal_ax.legend()

        contrast_ax.plot(
            points["g_resampled_corr"],
            points["contrast"],
            "-",
            color="#009E73",
            label="fitted/resampled contrast",
        )
        if not fit.empty:
            peak = fit.iloc[0]
            contrast_ax.plot(peak["g_peak_corr_mTm"], peak["signal_peak"], "*", color="black", markersize=9, label="peak")
        contrast_ax.set_ylabel(f"Contrast N={n_high} − N={n_low}")
        contrast_ax.set_xlabel("Common corrected gradient [mT/m]")
        contrast_ax.set_title(f"{variant}: resampled contrast")
        contrast_ax.grid(alpha=0.25)
        if points.empty:
            contrast_ax.text(
                0.5,
                0.5,
                "No successful resampled contrast",
                ha="center",
                va="center",
                transform=contrast_ax.transAxes,
            )
        else:
            contrast_ax.legend()
    sheet_label = f" | {sheet}" if sheet is not None else ""
    fig.suptitle(f"{subject}{sheet_label} | {roi} | {direction} | $T_D$={td_ms:g} ms", fontsize=14)
    fig.tight_layout()
    _save_figure(fig, output_path, dpi=dpi)
    return fig


def export_all_signal_contrast_panels(
    masters: pd.DataFrame,
    contrast: pd.DataFrame,
    contrast_fits: pd.DataFrame,
    *,
    variants: Sequence[str],
    output_dir: str | Path,
    n_high: int = 8,
    n_low: int = 4,
    dpi: int = 180,
    progress_every: int = 25,
    clean_output: bool = True,
) -> pd.DataFrame:
    """Export one Figure-1b-style panel for every measured combination.

    The batch key includes the acquisition sheet so repeated acquisitions do
    not overwrite each other. The output manifest records which selected
    variants produced a successful paired fit for each combination.
    """

    keys = ["subj", "sheet", "roi", "direction", "td_ms"]
    _require_columns(contrast_fits, ["variant", "ok", *keys], "contrast_fits")
    selected_variants = [str(variant) for variant in variants]
    candidates = contrast_fits[contrast_fits["variant"].astype(str).isin(selected_variants)].copy()
    if candidates.empty:
        raise ValueError("No contrast-fit candidates match the selected variants")

    combinations = candidates[keys].drop_duplicates().sort_values(keys, kind="stable").reset_index(drop=True)

    def normalized_key(values: Sequence[object]) -> tuple[str, str, str, str, float]:
        return (
            str(values[0]),
            str(values[1]),
            str(values[2]),
            str(values[3]),
            float(values[4]),
        )

    master_candidates = masters[
        masters["variant"].astype(str).isin(selected_variants)
        & masters["row_kind"].eq("signal_rotated")
        & masters["stat"].eq("avg")
        & pd.to_numeric(masters["N"], errors="coerce").isin([int(n_high), int(n_low)])
    ]
    master_groups = {
        normalized_key(key): group
        for key, group in master_candidates.groupby(keys, sort=False, dropna=False)
    }
    contrast_groups = {
        normalized_key(key): group
        for key, group in contrast[
            contrast["variant"].astype(str).isin(selected_variants)
        ].groupby(keys, sort=False, dropna=False)
    }
    fit_groups = {
        normalized_key(key): group
        for key, group in candidates.groupby(keys, sort=False, dropna=False)
    }
    output_dir = Path(output_dir)
    if clean_output and output_dir.exists():
        for previous_figure in output_dir.rglob("*.png"):
            previous_figure.unlink()
        previous_manifest = output_dir / "manifest.csv"
        if previous_manifest.exists():
            previous_manifest.unlink()
        for metadata_name in ("Thumbs.db", ".DS_Store"):
            for metadata_file in output_dir.rglob(metadata_name):
                metadata_file.unlink()
        # Remove the now-empty ROI/direction directories from the former
        # nested gallery layout, while preserving the gallery root itself.
        previous_directories = sorted(
            (path for path in output_dir.rglob("*") if path.is_dir()),
            key=lambda path: len(path.parts),
            reverse=True,
        )
        for previous_directory in previous_directories:
            try:
                previous_directory.rmdir()
            except OSError:
                pass
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    total = len(combinations)
    for index, combination in combinations.iterrows():
        combination_key = normalized_key([combination[column] for column in keys])
        group = fit_groups[combination_key]
        successful = [
            variant
            for variant in selected_variants
            if bool(group.loc[group["variant"].astype(str).eq(variant), "ok"].fillna(False).any())
        ]
        missing_or_failed = [variant for variant in selected_variants if variant not in successful]
        td_tag = _safe_tag(f"{float(combination['td_ms']):g}")
        figure_path = (
            output_dir
            / f"sub-{_safe_tag(combination['subj'])}"
            / f"sheet-{_safe_tag(combination['sheet'])}"
            / (
                f"roi-{_safe_tag(combination['roi'])}"
                f"__dir-{_safe_tag(combination['direction'])}"
                f"__td-{td_tag}ms.png"
            )
        )
        figure = plot_signal_contrast_example(
            master_groups.get(combination_key, masters.iloc[0:0]),
            contrast_groups.get(combination_key, contrast.iloc[0:0]),
            group,
            variants=selected_variants,
            subject=str(combination["subj"]),
            sheet=str(combination["sheet"]),
            roi=str(combination["roi"]),
            direction=str(combination["direction"]),
            td_ms=float(combination["td_ms"]),
            n_high=n_high,
            n_low=n_low,
            output_path=figure_path,
            dpi=dpi,
        )
        plt.close(figure)
        manifest_rows.append(
            {
                **combination.to_dict(),
                "selected_variants": "|".join(selected_variants),
                "successful_variants": "|".join(successful),
                "missing_or_failed_variants": "|".join(missing_or_failed),
                "all_selected_variants_successful": not missing_or_failed,
                "figure_path": str(figure_path.relative_to(output_dir)),
            }
        )
        completed = index + 1
        if progress_every > 0 and (completed % int(progress_every) == 0 or completed == total):
            print(f"Exported {completed}/{total} signal/contrast panels")

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(output_dir / "manifest.csv", index=False)
    return manifest


def plot_contrast_lcf_grid(
    contrast: pd.DataFrame,
    contrast_fits: pd.DataFrame,
    tc_data: pd.DataFrame,
    tc_summary: pd.DataFrame,
    *,
    variants: Sequence[str],
    subject: str | None,
    roi: str | None = None,
    directions: Sequence[str] = DEFAULT_DIRECTIONS,
    lcf_limits_um: tuple[float, float] | None = (3.0, 12.0),
    output_path: str | Path | None = None,
) -> plt.Figure:
    """Plot Figure 2-style contrast curves with peak-time insets."""

    subject, roi, _, _ = _resolve_example(
        contrast,
        subject=subject,
        roi=roi,
        direction=str(directions[0]),
        td_ms=None,
    )
    fig, axes = plt.subplots(
        len(variants),
        len(directions),
        figsize=(6.4 * len(directions), 4.8 * len(variants)),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    all_td = sorted(pd.to_numeric(contrast["td_ms"], errors="coerce").dropna().unique())
    colors = {td: plt.cm.viridis(index / max(1, len(all_td) - 1)) for index, td in enumerate(all_td)}
    for row_index, variant in enumerate(variants):
        for column_index, direction in enumerate(directions):
            ax = axes[row_index, column_index]
            fits = contrast_fits[
                contrast_fits["variant"].eq(variant)
                & contrast_fits["subj"].eq(subject)
                & contrast_fits["roi"].eq(roi)
                & contrast_fits["direction"].eq(direction)
                & contrast_fits["ok"].fillna(False)
            ].sort_values("td_ms")
            for fit_row in fits.itertuples(index=False):
                td_value = float(fit_row.td_ms)
                color = colors[td_value]
                curve = contrast[
                    contrast["variant"].eq(variant)
                    & contrast["subj"].eq(subject)
                    & contrast["roi"].eq(roi)
                    & contrast["direction"].eq(direction)
                    & np.isclose(pd.to_numeric(contrast["td_ms"], errors="coerce"), td_value)
                ].copy()
                curve = curve[np.isfinite(curve["lcf_um"])].sort_values("lcf_um")
                ax.plot(curve["lcf_um"], curve["contrast"], color=color, label=f"{td_value:g} ms")
            ax.set_title(f"{variant} | {direction}")
            ax.set_xlabel(r"Center filter length $l_{c,f}$ [$\mu$m]")
            ax.set_ylabel("OGSE contrast")
            ax.grid(alpha=0.25)
            if lcf_limits_um is not None:
                ax.set_xlim(*lcf_limits_um)
            if len(fits):
                ax.legend(title=r"$T_D$", fontsize=8, loc="lower left", frameon=False)

            inset = ax.inset_axes([0.58, 0.55, 0.38, 0.38])
            data = tc_data[
                tc_data["variant"].eq(variant)
                & tc_data["subj"].eq(subject)
                & tc_data["roi"].eq(roi)
                & tc_data["direction"].eq(direction)
            ]
            inset.plot(data["td_ms"], data["tc_peak_ms"], "o", color="black", markersize=3)
            summary = tc_summary[
                tc_summary["variant"].eq(variant)
                & tc_summary["subj"].eq(subject)
                & tc_summary["roi"].eq(roi)
                & tc_summary["direction"].eq(direction)
                & tc_summary["ok"].fillna(False)
            ]
            if not summary.empty and len(data):
                row = summary.iloc[0]
                grid = np.linspace(float(data["td_ms"].min()), float(data["td_ms"].max()), 200)
                inset.plot(grid, _tc_pseudohuber(grid, row["c_ms"], row["delta_ms"], row["alpha_macro"]), color="black")
            inset.set_xlabel(r"$T_D$", fontsize=7)
            inset.set_ylabel(r"$\tau_c$", fontsize=7)
            inset.tick_params(labelsize=7)
            inset.grid(alpha=0.2)
    fig.suptitle(f"{subject} | {roi}: NOGSE-like contrast and transition fits", fontsize=15)
    fig.tight_layout()
    _save_figure(fig, output_path)
    return fig


def _subject_markers(subjects: Sequence[str]) -> dict[str, str]:
    markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "*")
    return {subject: markers[index % len(markers)] for index, subject in enumerate(subjects)}


def plot_metric_comparison(
    metrics: pd.DataFrame,
    *,
    variants: Sequence[str],
    rois: Sequence[str] = CC_ROIS,
    directions: Sequence[str] = DEFAULT_DIRECTIONS,
    output_path: str | Path | None = None,
) -> plt.Figure:
    """Plot Figure 3-style alpha and delta summaries with variant bars and subject points."""

    metrics = metrics[metrics["ok"].fillna(False)].copy()
    subjects = sorted(metrics["subj"].unique())
    markers = _subject_markers(subjects)
    variant_colors = {variant: plt.cm.Set2(index % 8) for index, variant in enumerate(variants)}
    figure, axes = plt.subplots(2, len(directions), figsize=(6.3 * len(directions), 8.0), squeeze=False, sharex="col")
    plot_specs = [("delta_ms", r"Transition time $\delta$ [ms]"), ("alpha_macro", r"Macroscopic tortuosity $\alpha$")]
    x_base = np.arange(len(rois), dtype=float)
    width = 0.78 / max(1, len(variants))
    for row_index, (column, label) in enumerate(plot_specs):
        for direction_index, direction in enumerate(directions):
            ax = axes[row_index, direction_index]
            for variant_index, variant in enumerate(variants):
                offset = (variant_index - (len(variants) - 1) / 2.0) * width
                means, sems = [], []
                for roi in rois:
                    values = pd.to_numeric(
                        metrics[
                            metrics["variant"].eq(variant)
                            & metrics["direction"].eq(direction)
                            & metrics["roi"].eq(roi)
                        ][column],
                        errors="coerce",
                    ).dropna()
                    means.append(float(values.mean()) if len(values) else np.nan)
                    sems.append(float(values.sem()) if len(values) > 1 else 0.0)
                positions = x_base + offset
                ax.bar(
                    positions,
                    means,
                    width=width * 0.9,
                    yerr=sems,
                    color=variant_colors[variant],
                    alpha=0.75,
                    capsize=3,
                    label=variant,
                )
                for roi_index, roi in enumerate(rois):
                    values = metrics[
                        metrics["variant"].eq(variant)
                        & metrics["direction"].eq(direction)
                        & metrics["roi"].eq(roi)
                    ][["subj", column]].dropna()
                    for point_index, point in enumerate(values.itertuples(index=False)):
                        jitter = (point_index - (len(values) - 1) / 2.0) * width * 0.06
                        ax.plot(
                            positions[roi_index] + jitter,
                            getattr(point, column),
                            marker=markers[str(point.subj)],
                            markerfacecolor="white",
                            markeredgecolor=variant_colors[variant],
                            linestyle="none",
                            markersize=5,
                        )
            ax.set_title(direction)
            ax.set_ylabel(label)
            ax.set_xticks(x_base)
            ax.set_xticklabels(rois, rotation=35, ha="right")
            ax.grid(axis="y", alpha=0.25)
            if row_index == 0 and direction_index == 0:
                ax.legend(title="Variant", fontsize=8)
    subject_handles = [
        Line2D([], [], marker=markers[subject], linestyle="none", color="black", markerfacecolor="white", label=subject)
        for subject in subjects
    ]
    if subject_handles:
        figure.legend(
            handles=subject_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.955),
            ncol=min(6, len(subject_handles)),
            title="Subjects",
        )
    figure.suptitle("ESMRMB Figure 3 analogue: ROI-variant comparison", y=0.995, fontsize=15)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    _save_figure(figure, output_path)
    return figure


def plot_alpha_delta_scatter(
    metrics: pd.DataFrame,
    *,
    variants: Sequence[str],
    rois: Sequence[str] = CC_ROIS,
    directions: Sequence[str] = DEFAULT_DIRECTIONS,
    output_path: str | Path | None = None,
) -> plt.Figure:
    """Plot Figure 4-style alpha-versus-transition-time panels."""

    metrics = metrics[metrics["ok"].fillna(False)].copy()
    subjects = sorted(metrics["subj"].unique())
    markers = _subject_markers(subjects)
    roi_colors = {roi: plt.cm.tab10(index % 10) for index, roi in enumerate(rois)}
    figure, axes = plt.subplots(
        len(variants),
        len(directions),
        figsize=(6.0 * len(directions), 4.5 * len(variants)),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    for row_index, variant in enumerate(variants):
        for column_index, direction in enumerate(directions):
            ax = axes[row_index, column_index]
            panel = metrics[metrics["variant"].eq(variant) & metrics["direction"].eq(direction)]
            for point in panel.itertuples(index=False):
                ax.scatter(
                    point.alpha_macro,
                    point.delta_ms,
                    color=roi_colors.get(str(point.roi), "gray"),
                    marker=markers[str(point.subj)],
                    s=48,
                )
            ax.set_title(f"{variant} | {direction}")
            ax.set_xlabel(r"Macroscopic tortuosity $\alpha$")
            ax.set_ylabel(r"Transition time $\delta$ [ms]")
            ax.grid(alpha=0.25)
    roi_handles = [Line2D([], [], marker="o", linestyle="none", color=roi_colors[roi], label=roi) for roi in rois]
    subject_handles = [
        Line2D([], [], marker=markers[subject], linestyle="none", color="black", label=subject) for subject in subjects
    ]
    figure.legend(
        handles=roi_handles + subject_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=min(6, len(roi_handles + subject_handles)),
    )
    figure.suptitle("ESMRMB Figure 4 analogue: transition time versus macroscopic tortuosity", y=0.995, fontsize=15)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    _save_figure(figure, output_path)
    return figure


def compare_with_reference_metrics(
    metrics: pd.DataFrame,
    *,
    reference_variant: str,
    value_columns: Sequence[str] = ("alpha_macro", "delta_ms"),
) -> pd.DataFrame:
    """Return paired absolute and relative differences against one reference variant."""

    keys = ["subj", "roi", "direction"]
    reference = metrics[metrics["variant"].eq(reference_variant)][keys + list(value_columns)].copy()
    rows: list[pd.DataFrame] = []
    for variant in _ordered_unique(metrics["variant"]):
        candidate = metrics[metrics["variant"].eq(variant)][keys + list(value_columns)].copy()
        paired = candidate.merge(reference, on=keys, suffixes=("", "_reference"), validate="one_to_one")
        paired.insert(0, "variant", variant)
        paired.insert(1, "reference_variant", reference_variant)
        for column in value_columns:
            paired[f"{column}_difference"] = paired[column] - paired[f"{column}_reference"]
            paired[f"{column}_percent_difference"] = 100.0 * paired[f"{column}_difference"] / paired[f"{column}_reference"]
        rows.append(paired)
    return pd.concat(rows, ignore_index=True)
