from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


VALID_PARAMETER_MODES = ("fixed", "free", "global_td", "global_contrast")


@dataclass(frozen=True)
class FitParameterConfig:
    name: str
    mode: str
    init: float
    fixed: float | None
    bounds: tuple[float, float]
    log: bool = False


def finite_or_none(value: float | None) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except Exception:
        return None
    return out if np.isfinite(out) else None


def normalize_parameter_mode(mode: str, *, param_name: str) -> str:
    out = str(mode).strip()
    if out not in VALID_PARAMETER_MODES:
        raise ValueError(
            f"Invalid mode for {param_name}: {mode!r}. "
            f"Expected one of {list(VALID_PARAMETER_MODES)}."
        )
    return out


def fixed_parameter_value(config: FitParameterConfig) -> float:
    fixed = finite_or_none(config.fixed)
    if fixed is not None:
        return float(fixed)
    init = finite_or_none(config.init)
    if init is not None:
        return float(init)
    raise ValueError(f"{config.name} mode is fixed but no finite fixed/init value was provided.")


def make_parameter_config(
    *,
    name: str,
    mode: str,
    init: float,
    fixed: float | None,
    bounds: tuple[float, float],
    log: bool = False,
) -> FitParameterConfig:
    normalized_mode = normalize_parameter_mode(mode, param_name=name)
    if log and float(bounds[0]) <= 0:
        raise ValueError(f"{name} uses log-optimization and needs a positive lower bound. Received {bounds}.")
    return FitParameterConfig(
        name=str(name),
        mode=normalized_mode,
        init=float(init),
        fixed=finite_or_none(fixed),
        bounds=(float(bounds[0]), float(bounds[1])),
        log=bool(log),
    )


def to_optimization_values(values: np.ndarray, names: Sequence[str], log_params: set[str]) -> np.ndarray:
    out = []
    for name, value in zip(names, values):
        out.append(float(np.log(value)) if name in log_params else float(value))
    return np.asarray(out, dtype=float)


def from_optimization_values(values: np.ndarray, names: Sequence[str], log_params: set[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, value in zip(names, values):
        out[name] = float(np.exp(value)) if name in log_params else float(value)
    return out


def parameter_mode_summary(
    param_configs: dict[str, FitParameterConfig],
    mode: str,
    *,
    ordered_names: Sequence[str],
) -> str:
    names = [name for name in ordered_names if param_configs[name].mode == mode]
    return ",".join(names)
