from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


VALID_G_TYPES = {"g", "g_max", "g_lin_max", "g_thorsten"}
VALID_AXIS_BASES = {
    "bvalue",
    "bvalue_g",
    "bvalue_g_lin_max",
    "bvalue_thorsten",
    *VALID_G_TYPES,
}

_AXIS_TO_GRADIENT_BASE = {
    "bvalue": "g",
    "g": "g",
    "bvalue_g": "g",
    "g_max": "g_max",
    "g_lin_max": "g_lin_max",
    "bvalue_g_lin_max": "g_lin_max",
    "g_thorsten": "g_thorsten",
    "bvalue_thorsten": "g_thorsten",
}

_AXIS_TO_BVALUE_COLUMN = {
    "bvalue": "bvalue",
    "bvalue_g": "bvalue_g",
    "bvalue_g_lin_max": "bvalue_g_lin_max",
    "bvalue_thorsten": "bvalue_thorsten",
}


@dataclass(frozen=True)
class AxisBundle:
    requested_axis: str
    axis_base: str
    gradient_base: str
    side: int | None
    gradient_raw: np.ndarray
    gradient_corr: np.ndarray
    bvalue_raw: np.ndarray | None
    bvalue_corr: np.ndarray | None
    axis_raw: np.ndarray
    axis_corr: np.ndarray


def b_from_g(
    g_mTm: np.ndarray,
    *,
    N: float,
    gamma: float,
    delta_ms: float,
    delta_app_ms: float,
    g_type: str,
) -> np.ndarray:
    """
    Compute b = N * gamma^2 * delta^2 * delta_app * g^2 / 1e9.

    Expected units:
      gamma: 1/(ms*mT)
      delta_ms, delta_app_ms: ms
      g_mTm: mT/m
      b: s/mm^2
    """
    g = np.asarray(g_mTm, dtype=float)
    if g_type not in VALID_G_TYPES:
        raise ValueError(f"b_from_g: unrecognized g_type {g_type!r}. Allowed values: {sorted(VALID_G_TYPES)}.")
    return N * (gamma**2) * (delta_ms**2) * delta_app_ms * (g**2) / 1e9


# NOGSE: b = b_per_g2 * g^2, with b_per_g2 in s/mm^2 per (mT/m)^2. The OGSE formula above
# (N, delta, Delta_app) does not apply to NOGSE. Pulseq NOGSE uses gamma/2pi = 42.576 MHz/T.
NOGSE_B_PER_G2_KEY = "b_per_g2_s_per_mm2_per_mt2_m2"
NOGSE_GAMMA_RAD_PER_MS_MT = 2.0 * np.pi * 42.576


def nogse_trapezoid_b_per_g2(
    *,
    N: float,
    TN_ms: float,
    x_ms: float,
    ramp_time_ms: float,
    gamma: float = NOGSE_GAMMA_RAD_PER_MS_MT,
) -> float:
    """Exact ``b / g^2`` of the played trapezoidal NOGSE waveform.

    The waveform is the one used by the Pulseq NOGSE generator: ``N - 1``
    alternating CPMG lobes of total duration ``x``, followed by a Hahn lobe of
    duration ``y = TN - (N - 1) x``. Both blocks are split around the refocusing
    pulses, and every lobe has the same linear rise/fall time.

    Units are ms, mT/m and rad/(ms*mT); the result is s/mm^2 per (mT/m)^2.
    """
    n_float = float(N)
    n = int(round(n_float))
    TN_ms, x_ms, ramp_time_ms, gamma = map(float, (TN_ms, x_ms, ramp_time_ms, gamma))
    if not np.isfinite(n_float) or not np.isclose(n_float, n) or n < 2:
        raise ValueError(f"NOGSE N must be an integer >= 2, got {N!r}.")
    if not all(np.isfinite(v) and v > 0.0 for v in (TN_ms, x_ms, ramp_time_ms, gamma)):
        raise ValueError(
            f"NOGSE TN, x, ramp time and gamma must be positive; got "
            f"TN={TN_ms}, x={x_ms}, ramp={ramp_time_ms}, gamma={gamma}."
        )

    y_ms = TN_ms - (n - 1) * x_ms
    half_x = x_ms / 2.0
    half_y = y_ms / 2.0
    if y_ms <= 0.0 or min(half_x, half_y) < 2.0 * ramp_time_ms:
        raise ValueError(
            f"Invalid trapezoidal NOGSE timing N={n}, TN={TN_ms} ms, x={x_ms} ms, "
            f"y={y_ms} ms, ramp={ramp_time_ms} ms: each half-lobe must be at least 2*ramp."
        )

    # (duration, sign) pairs, matching the Pulseq generator. A unit gradient is
    # sufficient because b scales quadratically with gradient amplitude.
    lobes: list[tuple[float, float]] = [(half_x, 1.0)]
    sign = -1.0
    for _ in range(n - 2):
        lobes.append((x_ms, sign))
        sign *= -1.0
    lobes.extend([(half_x, sign), (half_y, sign), (half_y, -sign)])

    # Integrate q(t)^2 exactly on each linear-gradient segment. If
    # g(t)=g0+a*t, then q(t)=q0+gamma*(g0*t+a*t^2/2).
    q0 = 0.0
    b_ms3 = 0.0
    for duration, lobe_sign in lobes:
        flat = duration - 2.0 * ramp_time_ms
        segments = [(ramp_time_ms, 0.0, lobe_sign)]
        if flat > 0.0:
            segments.append((flat, lobe_sign, lobe_sign))
        segments.append((ramp_time_ms, lobe_sign, 0.0))
        for h, g0, g1 in segments:
            a = (g1 - g0) / h
            c0 = q0
            c1 = gamma * g0
            c2 = gamma * a / 2.0
            b_ms3 += (
                c0**2 * h
                + c0 * c1 * h**2
                + (c1**2 + 2.0 * c0 * c2) * h**3 / 3.0
                + c1 * c2 * h**4 / 2.0
                + c2**2 * h**5 / 5.0
            )
            q0 += gamma * (g0 + g1) * h / 2.0

    if not np.isclose(q0, 0.0, rtol=0.0, atol=1e-9):
        raise ValueError(f"NOGSE waveform is not balanced (final q={q0:.9g}).")
    return b_ms3 / 1e9


def nogse_rect_b_per_g2(*, N: float, TN_ms: float, x_ms: float, gamma: float) -> float:
    """
    Rectangular-lobe NOGSE b / g^2 = gamma^2 [(N - 1) x^3 + y^3] / 12 / 1e9, y = TN - (N - 1) x
    (the b-value of M_nogse_free).

    Expected units:
      gamma: 1/(ms*mT)
      TN_ms, x_ms: ms
      result: s/mm^2 per (mT/m)^2
    """
    N, TN_ms, x_ms = float(N), float(TN_ms), float(x_ms)
    y_ms = TN_ms - (N - 1.0) * x_ms
    if not (np.isfinite(y_ms) and y_ms >= 0.0 and x_ms >= 0.0 and N >= 1.0):
        raise ValueError(f"Invalid NOGSE timing N={N}, TN={TN_ms} ms, x={x_ms} ms (y={y_ms} ms).")
    return (gamma**2) * ((N - 1.0) * x_ms**3 + y_ms**3) / 12.0 / 1e9


def b_from_g_nogse(g_mTm: np.ndarray, *, b_per_g2: float) -> np.ndarray:
    """NOGSE b (s/mm^2) = b_per_g2 [s/mm^2 per (mT/m)^2] * g^2, g in mT/m."""
    return float(b_per_g2) * np.asarray(g_mTm, dtype=float) ** 2


def g_from_b_nogse(bvalue: np.ndarray, *, b_per_g2: float) -> np.ndarray:
    """NOGSE g (mT/m) = sqrt(b / b_per_g2); b = 0 gives g = 0."""
    if not np.isfinite(float(b_per_g2)) or float(b_per_g2) <= 0.0:
        raise ValueError(f"NOGSE b_per_g2 must be positive, got {b_per_g2!r}.")
    b = np.asarray(bvalue, dtype=float)
    g = np.sqrt(np.clip(b, 0.0, None) / float(b_per_g2))
    g[np.isclose(b, 0.0)] = 0.0
    return g


def normalize_axis_base(axis: str) -> str:
    base = str(axis).strip()
    if base.endswith("_1") or base.endswith("_2"):
        base = base[:-2]
    if base not in VALID_AXIS_BASES:
        raise ValueError(f"Unrecognized axis {axis!r}. Allowed values: {sorted(VALID_AXIS_BASES)}.")
    return base


def split_axis_side(axis: str) -> tuple[str, int | None]:
    raw = str(axis).strip()
    if raw.endswith("_1") or raw.endswith("_2"):
        return normalize_axis_base(raw[:-2]), int(raw[-1])
    return normalize_axis_base(raw), None


def axis_uses_bvalue(axis: str) -> bool:
    return normalize_axis_base(axis).startswith("bvalue")


def gradient_base_for_axis(axis: str) -> str:
    base = normalize_axis_base(axis)
    return _AXIS_TO_GRADIENT_BASE[base]


def axes_share_gradient_family(axis_a: str, axis_b: str) -> bool:
    return gradient_base_for_axis(axis_a) == gradient_base_for_axis(axis_b)


def default_plot_axis_for_fit(axis: str, *, side: int | None = None) -> str:
    base = normalize_axis_base(axis)
    if side is None:
        return base
    return f"{base}_{int(side)}"


def coerce_positive_correction_factor(value: float | None) -> float:
    if value is None:
        return 1.0
    out = float(value)
    if not np.isfinite(out) or out <= 0.0:
        return 1.0
    return out


def gradient_column_name(axis: str, *, side: int | None = None) -> str:
    base = gradient_base_for_axis(axis)
    return base if side is None else f"{base}_{int(side)}"


def bvalue_column_name(axis: str, *, side: int | None = None) -> str | None:
    base = normalize_axis_base(axis)
    col = _AXIS_TO_BVALUE_COLUMN.get(base)
    if col is None:
        return None
    return col if side is None else f"{col}_{int(side)}"


def extract_gradient_array(
    df: pd.DataFrame,
    *,
    axis: str,
    side: int | None = None,
) -> np.ndarray:
    col = gradient_column_name(axis, side=side)
    if col not in df.columns:
        raise KeyError(f"Missing required gradient column {col!r}. Columns={list(df.columns)}")
    arr = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
    return arr


def _extract_stored_bvalue_array(
    df: pd.DataFrame,
    *,
    axis: str,
    side: int | None = None,
) -> np.ndarray | None:
    col = bvalue_column_name(axis, side=side)
    if col is None or col not in df.columns:
        return None
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def bvalue_from_gradient(
    gradient_values: np.ndarray,
    *,
    axis: str,
    N: float | None,
    gamma: float,
    delta_ms: float | None,
    Delta_app_ms: float | None,
) -> np.ndarray:
    if N is None or delta_ms is None or Delta_app_ms is None:
        raise ValueError(
            "Computing a b-value axis from corrected gradients requires N, delta_ms, and Delta_app_ms."
        )
    return b_from_g(
        np.asarray(gradient_values, dtype=float),
        N=float(N),
        gamma=float(gamma),
        delta_ms=float(delta_ms),
        delta_app_ms=float(Delta_app_ms),
        g_type=gradient_base_for_axis(axis),
    )


def axis_from_gradient(
    gradient_values: np.ndarray,
    *,
    axis: str,
    N: float | None,
    gamma: float,
    delta_ms: float | None,
    Delta_app_ms: float | None,
) -> np.ndarray:
    base = normalize_axis_base(axis)
    gradient_values = np.asarray(gradient_values, dtype=float)
    if axis_uses_bvalue(base):
        return bvalue_from_gradient(
            gradient_values,
            axis=base,
            N=N,
            gamma=gamma,
            delta_ms=delta_ms,
            Delta_app_ms=Delta_app_ms,
        )
    return gradient_values


def build_axis_bundle(
    df: pd.DataFrame,
    *,
    axis: str,
    correction_factor: float = 1.0,
    gamma: float = 267.5221900,
    N: float | None = None,
    delta_ms: float | None = None,
    Delta_app_ms: float | None = None,
    side: int | None = None,
) -> AxisBundle:
    axis_base, axis_side = split_axis_side(axis)
    resolved_side = axis_side if side is None else int(side)
    if axis_side is not None and side is not None and int(axis_side) != int(side):
        raise ValueError(f"Axis {axis!r} and explicit side={side} disagree.")

    f_corr = coerce_positive_correction_factor(correction_factor)
    gradient_raw = extract_gradient_array(df, axis=axis_base, side=resolved_side)
    gradient_corr = gradient_raw * float(f_corr)

    bvalue_raw = _extract_stored_bvalue_array(df, axis=axis_base, side=resolved_side)
    bvalue_corr: np.ndarray | None = None

    if axis_uses_bvalue(axis_base):
        if bvalue_raw is None:
            if N is None or delta_ms is None or Delta_app_ms is None:
                raise ValueError(
                    f"Axis {axis!r} needs N, delta_ms, and Delta_app_ms to derive a b-value."
                )
            bvalue_raw = bvalue_from_gradient(
                gradient_raw,
                axis=axis_base,
                N=N,
                gamma=gamma,
                delta_ms=delta_ms,
                Delta_app_ms=Delta_app_ms,
            )
        # A b-value axis is the b the scanner applied: the gradient-correction factor only
        # rescales g for the OGSE/NOGSE signal models, so the b axis is never corrected.
        bvalue_corr = np.asarray(bvalue_raw, dtype=float)

    axis_raw = np.asarray(bvalue_raw, dtype=float) if axis_uses_bvalue(axis_base) else np.asarray(gradient_raw, dtype=float)
    axis_corr = np.asarray(bvalue_corr, dtype=float) if axis_uses_bvalue(axis_base) else np.asarray(gradient_corr, dtype=float)

    return AxisBundle(
        requested_axis=str(axis),
        axis_base=axis_base,
        gradient_base=gradient_base_for_axis(axis_base),
        side=resolved_side,
        gradient_raw=np.asarray(gradient_raw, dtype=float),
        gradient_corr=np.asarray(gradient_corr, dtype=float),
        bvalue_raw=None if bvalue_raw is None else np.asarray(bvalue_raw, dtype=float),
        bvalue_corr=None if bvalue_corr is None else np.asarray(bvalue_corr, dtype=float),
        axis_raw=axis_raw,
        axis_corr=axis_corr,
    )


__all__ = [
    "AxisBundle",
    "VALID_AXIS_BASES",
    "VALID_G_TYPES",
    "axes_share_gradient_family",
    "axis_from_gradient",
    "axis_uses_bvalue",
    "b_from_g",
    "b_from_g_nogse",
    "build_axis_bundle",
    "bvalue_column_name",
    "bvalue_from_gradient",
    "coerce_positive_correction_factor",
    "default_plot_axis_for_fit",
    "extract_gradient_array",
    "gradient_base_for_axis",
    "g_from_b_nogse",
    "gradient_column_name",
    "NOGSE_B_PER_G2_KEY",
    "NOGSE_GAMMA_RAD_PER_MS_MT",
    "nogse_rect_b_per_g2",
    "nogse_trapezoid_b_per_g2",
    "normalize_axis_base",
    "split_axis_side",
]
