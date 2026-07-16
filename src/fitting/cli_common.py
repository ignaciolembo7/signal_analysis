from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from data_processing.master_table import load_master_table, select_signal, filter_master_rows, write_master_table
from data_processing.io import fit_params_output_basename, write_table_outputs
from data_processing.master_table import append_master_rows
from fitting.parameter_modes import (
    FitParameterConfig,
    fixed_parameter_value,
    make_parameter_config,
    normalize_parameter_mode,
)


MASTER_SELECTOR_ARGS = (
    ("analysis_id", "analysis_id"),
    ("subj", "subj"),
    ("sheet", "sheet"),
    ("roi", "roi"),
    ("direction", "direction"),
    ("stat", "stat"),
    ("source_file", "source_file"),
)


PARAM_ALIASES = {
    "m0": "M0",
    "M0": "M0",
    "d0": "D0_m2_ms",
    "D0": "D0_m2_ms",
    "D0_m2_ms": "D0_m2_ms",
    "tc": "tc_ms",
    "tc_ms": "tc_ms",
    "alpha": "alpha",
    "c": "C",
    "C": "C",
    "rn": "RN",
    "RN": "RN",
    "g0": "g0_mTm",
    "g0_mTm": "g0_mTm",
}


@dataclass(frozen=True)
class CommonParameterPlan:
    configs: dict[str, FitParameterConfig]

    def mode(self, name: str, default: str = "free") -> str:
        config = self.configs.get(canonical_param_name(name))
        return default if config is None else config.mode

    def init(self, name: str, default: float) -> float:
        config = self.configs.get(canonical_param_name(name))
        return float(default) if config is None else float(config.init)

    def fixed(self, name: str, default: float | None = None) -> float | None:
        config = self.configs.get(canonical_param_name(name))
        if config is None:
            return default
        if config.mode != "fixed":
            return default
        return fixed_parameter_value(config)

    def bounds(self, name: str, default: tuple[float, float]) -> tuple[float, float]:
        config = self.configs.get(canonical_param_name(name))
        return default if config is None else config.bounds

    def global_params(self) -> list[str]:
        return [name for name, config in self.configs.items() if config.mode in {"global_td", "global_contrast"}]


def canonical_param_name(name: str) -> str:
    key = str(name).strip()
    try:
        return PARAM_ALIASES[key]
    except KeyError as exc:
        raise ValueError(f"Unsupported parameter name {name!r}. Allowed values: {sorted(PARAM_ALIASES)}.") from exc


def split_cli_values(values: Sequence[str] | None) -> list[str] | None:
    if values is None:
        return None
    if isinstance(values, str):
        return str(values).replace(",", " ").split() or None
    out: list[str] = []
    for value in values:
        out.extend(str(value).replace(",", " ").split())
    return out or None


def add_master_source_args(
    parser: argparse.ArgumentParser,
    *,
    default_row_kind: str | None = None,
    include_roi: bool = True,
    include_direction: bool = True,
    include_stat: bool = True,
    include_td_ms: bool = True,
    include_N: bool = True,
    include_Hz: bool = True,
) -> None:
    parser.add_argument("--master-parquet", type=Path, default=None, help="Read input rows from the master table.")
    parser.add_argument("--row-kind", default=default_row_kind, help="Master row_kind selector. Defaults to the fitter's natural input kind.")
    parser.add_argument("--analysis-id", action="append", default=None, help="Master analysis_id selector. Can be repeated or comma-separated.")
    parser.add_argument("--subj", action="append", default=None, help="Master subj selector. Can be repeated or comma-separated.")
    parser.add_argument("--sheet", action="append", default=None, help="Master sheet selector. Can be repeated or comma-separated.")
    if include_roi:
        parser.add_argument("--roi", action="append", default=None, help="Master ROI selector. Can be repeated or comma-separated.")
    if include_direction:
        parser.add_argument("--direction", action="append", default=None, help="Master direction selector. Can be repeated or comma-separated.")
    parser.add_argument("--source-file", action="append", default=None, help="Master source_file selector. Can be repeated or comma-separated.")
    if include_stat:
        parser.add_argument("--stat", default=None, help="Master stat selector.")
    if include_td_ms:
        parser.add_argument("--td_ms", type=float, default=None, help="Master td_ms selector.")
    if include_N:
        parser.add_argument("--N", type=float, default=None, help="Master N selector.")
    if include_Hz:
        parser.add_argument("--Hz", type=float, default=None, help="Master Hz selector.")


def add_parameter_mode_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--param-mode",
        "--param_mode",
        action="append",
        default=None,
        metavar="PARAM=MODE",
        help="Unified parameter mode. MODE is fixed, free, global_td, or global_contrast. Can be repeated.",
    )
    parser.add_argument(
        "--param-init",
        "--param_init",
        action="append",
        default=None,
        metavar="PARAM=VALUE",
        help="Unified parameter initial value. Can be repeated.",
    )
    parser.add_argument(
        "--param-fixed",
        "--param_fixed",
        action="append",
        default=None,
        metavar="PARAM=VALUE",
        help="Unified fixed parameter value. Implies --param-mode PARAM=fixed if no mode is supplied.",
    )
    parser.add_argument(
        "--param-bounds",
        "--param_bounds",
        action="append",
        default=None,
        metavar="PARAM=LOW:HIGH",
        help="Unified parameter bounds. Can be repeated.",
    )


def add_fit_master_output_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--append-fit-params-to-master",
        action="store_true",
        help="Also append fit_params rows to --master-parquet with row_kind='fit_params'.",
    )


def _parse_assignments(values: Sequence[str] | None, *, option: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in values or []:
        for item in str(raw).replace(",", " ").split():
            if "=" not in item:
                raise ValueError(f"{option} expects PARAM=VALUE, received {item!r}.")
            name, value = item.split("=", 1)
            out[canonical_param_name(name)] = value
    return out


def _parse_bounds(values: Sequence[str] | None) -> dict[str, tuple[float, float]]:
    raw = _parse_assignments(values, option="--param-bounds")
    out: dict[str, tuple[float, float]] = {}
    for name, text in raw.items():
        parts = str(text).replace(":", " ").split()
        if len(parts) != 2:
            raise ValueError(f"--param-bounds for {name} expects LOW:HIGH, received {text!r}.")
        lo, hi = float(parts[0]), float(parts[1])
        if not lo < hi:
            raise ValueError(f"--param-bounds for {name} requires LOW < HIGH, received {lo}, {hi}.")
        out[name] = (lo, hi)
    return out


def build_common_parameter_plan(
    args: argparse.Namespace,
    *,
    param_names: Sequence[str],
    default_modes: dict[str, str],
    default_inits: dict[str, float],
    default_bounds: dict[str, tuple[float, float]],
    log_params: Sequence[str] = (),
) -> CommonParameterPlan:
    modes = _parse_assignments(getattr(args, "param_mode", None), option="--param-mode")
    inits_raw = _parse_assignments(getattr(args, "param_init", None), option="--param-init")
    fixed_raw = _parse_assignments(getattr(args, "param_fixed", None), option="--param-fixed")
    bounds = _parse_bounds(getattr(args, "param_bounds", None))
    log_set = {canonical_param_name(name) for name in log_params}

    configs: dict[str, FitParameterConfig] = {}
    for raw_name in param_names:
        name = canonical_param_name(raw_name)
        mode = normalize_parameter_mode(modes.get(name, "fixed" if name in fixed_raw else default_modes.get(name, "free")), param_name=name)
        init = float(inits_raw.get(name, default_inits.get(name, np.nan)))
        fixed = float(fixed_raw[name]) if name in fixed_raw else None
        configs[name] = make_parameter_config(
            name=name,
            mode=mode,
            init=init,
            fixed=fixed,
            bounds=bounds.get(name, default_bounds.get(name, (-np.inf, np.inf))),
            log=name in log_set,
        )
    return CommonParameterPlan(configs=configs)


def master_selectors_from_args(args: argparse.Namespace) -> dict[str, object]:
    selectors: dict[str, object] = {}
    for arg_name, col_name in MASTER_SELECTOR_ARGS:
        values = split_cli_values(getattr(args, arg_name, None))
        if values is not None:
            selectors[col_name] = values
    for col_name in ("td_ms", "N", "Hz"):
        value = getattr(args, col_name, None)
        if value is not None:
            selectors[col_name] = float(value)
    row_kind = getattr(args, "row_kind", None)
    if row_kind:
        selectors["row_kind"] = str(row_kind)
    return selectors


def load_master_input(args: argparse.Namespace, *, default_row_kind: str, signal_rotated: bool | None = None) -> pd.DataFrame:
    master_path = getattr(args, "master_parquet", None)
    if master_path is None:
        raise ValueError("load_master_input requires --master-parquet.")
    master = load_master_table(master_path)
    selectors = master_selectors_from_args(args)
    row_kind = selectors.pop("row_kind", None) or default_row_kind
    if row_kind in {"signal", "signal_rotated"}:
        rotated = signal_rotated
        if rotated is None:
            rotated = row_kind == "signal_rotated"
        out = select_signal(master, rotated=rotated, **selectors)
    else:
        out = filter_master_rows(master, row_kind=row_kind, **selectors)
    if out.empty:
        raise ValueError(f"No master rows matched row_kind={row_kind!r} and selectors={selectors}.")
    return out


def _selector_mask(df: pd.DataFrame, selectors: dict[str, object]) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    for col, value in selectors.items():
        values = split_cli_values(value if isinstance(value, list) else [str(value)])
        if values is None or (len(values) == 1 and values[0].upper() == "ALL"):
            continue
        if len(values) == 1:
            try:
                numeric = float(values[0])
            except ValueError:
                numeric = None
            if numeric is not None and np.isfinite(numeric):
                mask &= np.isclose(pd.to_numeric(df[col], errors="coerce"), numeric, atol=1e-6)
                continue
        mask &= df[col].astype(str).isin({str(v) for v in values})
    return mask


def update_master_signal_correction_factors(args: argparse.Namespace, f_by_direction: dict[str, float] | None) -> None:
    if not f_by_direction:
        return
    master_path = getattr(args, "master_parquet", None)
    if master_path is None:
        return

    master_path = Path(master_path)
    master = load_master_table(master_path)
    selectors = master_selectors_from_args(args)
    mask = _selector_mask(master, selectors)
    if not mask.any():
        raise ValueError(f"No master rows matched selectors for grad_correction_factor update: {selectors}")

    direction = master.loc[mask, "direction"].astype(str)
    master.loc[mask, "grad_correction_factor"] = direction.map(lambda d: float(f_by_direction.get(str(d), 1.0))).to_numpy()
    write_master_table(master, master_path)
    print("Updated master grad_correction_factor:", master_path)


def signal_correction_factors_from_rows(df: pd.DataFrame) -> dict[str, float]:
    if "grad_correction_factor" not in df.columns:
        raise ValueError("Selected master rows do not contain grad_correction_factor. Run grad_correction first.")
    rows = df[["direction", "grad_correction_factor"]].copy()
    rows["grad_correction_factor"] = pd.to_numeric(rows["grad_correction_factor"], errors="coerce")
    missing = rows.loc[~np.isfinite(rows["grad_correction_factor"]), "direction"].dropna().astype(str).unique().tolist()
    if missing:
        raise ValueError(f"Missing grad_correction_factor in master rows for directions: {sorted(missing)}")
    grouped = rows.groupby("direction", as_index=False)["grad_correction_factor"].mean()
    return {str(d): float(f) for d, f in zip(grouped["direction"], grouped["grad_correction_factor"])}


def require_contrast_correction_factor_columns(df: pd.DataFrame) -> None:
    missing_cols = {"grad_correction_factor_1", "grad_correction_factor_2"} - set(df.columns)
    if missing_cols:
        raise ValueError(
            f"Selected contrast rows do not contain {sorted(missing_cols)}. Run grad_correction first."
        )
    missing = df[
        ~np.isfinite(pd.to_numeric(df["grad_correction_factor_1"], errors="coerce"))
        | ~np.isfinite(pd.to_numeric(df["grad_correction_factor_2"], errors="coerce"))
    ]
    if not missing.empty:
        cols = [c for c in ["subj", "sheet", "direction", "td_ms", "N_1", "N_2"] if c in missing.columns]
        preview = missing[cols].drop_duplicates().head(10).to_string(index=False)
        raise ValueError(f"Missing contrast gradient correction factors in master rows:\n{preview}")


def contrast_correction_factors_from_rows(df: pd.DataFrame, *, side: int | None = None) -> dict[str, float | tuple[float, float]]:
    require_contrast_correction_factor_columns(df)
    rows = df[["direction", "grad_correction_factor_1", "grad_correction_factor_2"]].copy()
    rows["grad_correction_factor_1"] = pd.to_numeric(rows["grad_correction_factor_1"], errors="coerce")
    rows["grad_correction_factor_2"] = pd.to_numeric(rows["grad_correction_factor_2"], errors="coerce")
    grouped = rows.groupby("direction", as_index=False)[["grad_correction_factor_1", "grad_correction_factor_2"]].mean()
    if side == 1:
        return {str(d): float(f) for d, f in zip(grouped["direction"], grouped["grad_correction_factor_1"])}
    if side == 2:
        return {str(d): float(f) for d, f in zip(grouped["direction"], grouped["grad_correction_factor_2"])}
    return {
        str(row["direction"]): (float(row["grad_correction_factor_1"]), float(row["grad_correction_factor_2"]))
        for _, row in grouped.iterrows()
    }


def append_fit_params_outputs(
    fit_params: pd.DataFrame,
    args: argparse.Namespace,
    *,
    fit_kind: str,
    model: str,
    source: str,
    out_basename: str | None = None,
) -> None:
    if fit_params is None or fit_params.empty:
        return

    rows = fit_params.copy()
    rows["fit_kind"] = str(fit_kind)
    rows["fit_model"] = str(model)
    rows["fit_source"] = str(source)
    for attr, col in [
        ("analysis_id", "analysis_id"),
        ("subj", "subj"),
        ("sheet", "sheet"),
        ("roi", "roi"),
        ("direction", "direction"),
        ("td_ms", "td_ms"),
        ("N", "N"),
        ("Hz", "Hz"),
    ]:
        if not hasattr(args, attr):
            continue
        value = getattr(args, attr)
        if value is None:
            continue
        if isinstance(value, list):
            values = split_cli_values(value)
            if values is None or len(values) != 1:
                continue
            value = values[0]
        if col not in rows.columns or rows[col].isna().all():
            rows[col] = value

    out_root = getattr(args, "out_root", None)
    if out_root is not None:
        if out_basename is None:
            out_basename = f"fit_params.{model}"
        out_path = Path(out_root) / f"{out_basename}.parquet"
        if out_path.exists():
            existing = pd.read_parquet(out_path)
            combined = pd.concat([existing, rows], ignore_index=True, sort=False).drop_duplicates().reset_index(drop=True)
        else:
            combined = rows.drop_duplicates().reset_index(drop=True)
        write_table_outputs(combined, out_path, xlsx_path=out_path.with_suffix(".xlsx"))
        print("Saved fit_params:", out_path)

    if bool(getattr(args, "append_fit_params_to_master", False)):
        master_path = getattr(args, "master_parquet", None)
        if master_path is None:
            raise ValueError("--append-fit-params-to-master requires --master-parquet.")
        append_master_rows(
            master_path,
            rows,
            row_kind="fit_params",
            analysis_id=None,
            out_path=master_path,
        )
        print("Appended fit_params to master:", master_path)
