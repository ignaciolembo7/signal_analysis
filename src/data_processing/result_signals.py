from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from data_processing.master_table import build_analysis_id_from_columns
from data_processing.schema import finalize_clean_signal_long
from fitting.b_from_g import b_from_g


def norm_key(value: str) -> str:
    return " ".join(str(value).strip().lower().split())


def row_get(row: pd.Series, keys: list[str], default=None):
    """Look up keys in row case/space-insensitively and return the first non-empty value."""
    norm_map = {norm_key(c): c for c in list(row.index)}
    for key in keys:
        col = norm_map.get(norm_key(key))
        if col is None:
            continue
        value = row[col]
        if value is not None and not (isinstance(value, float) and np.isnan(value)):
            return value
    return default


def to_float(value):
    return float(pd.to_numeric(value, errors="coerce"))


def ensure_numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def first_present(columns: list[str], candidates: list[str]) -> str | None:
    for col in candidates:
        if col in columns:
            return col
    return None


def sanitize_token(value) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "group"


def split_strip_tokens(values: list[str] | None) -> list[str]:
    tokens: list[str] = []
    for value in values or []:
        tokens.extend(t for t in re.split(r"[,\s]+", str(value)) if t)
    return tokens


def output_strip_tokens(cli_tokens: list[str] | None) -> list[str]:
    env_tokens = os.environ.get("NOGSE_OUTPUT_STEM_STRIP", "")
    return split_strip_tokens([env_tokens]) + split_strip_tokens(cli_tokens)


def output_stem_from_results_file(
    results_file: Path,
    *,
    grouped_g_curve: bool,
    strip_tokens: list[str],
) -> str:
    stem = results_file.stem
    if grouped_g_curve:
        stem = re.sub(r"_G\d+(?:p\d+|\.\d+)?(?=_)", "", stem, flags=re.IGNORECASE)
        stem = re.sub(r"_\d+(?=_results$)", "", stem, flags=re.IGNORECASE)
    for token in strip_tokens:
        stem = stem.replace(token, "")
    stem = re.sub(r"_+", "_", stem).strip("_")
    return sanitize_token(stem)


def truthy_env(name: str) -> bool:
    raw = os.environ.get(name, "")
    return str(raw).strip().lower() in {"1", "true", "yes", "y", "on"}


def finite_float(value) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def normalize_protocol_label(protocol, group_value, tn_value) -> str | None:
    if protocol is None or not str(protocol).strip():
        return None

    base = str(protocol).strip()
    base = re.sub(r"(^|_)BW\d+(?=_|$)", r"\1", base, flags=re.IGNORECASE)
    base = re.sub(r"(^|_)G\d+(?:p\d+|\.\d+)?(?=_|$)", r"\1", base, flags=re.IGNORECASE)

    group_float = finite_float(group_value)
    if group_float is not None:
        group_num = int(round(group_float))
        updated, nsubs = re.subn(r"Exp\d+", f"Exp{group_num:02d}", base, count=1, flags=re.IGNORECASE)
        base = updated if nsubs else f"{base}_group{group_num:03d}"

    tn_float = finite_float(tn_value)
    if tn_float is not None:
        base = re.sub(r"TN\d+(?:p\d+|\.\d+)?", f"TN{int(round(tn_float))}", base, count=1, flags=re.IGNORECASE)

    base = re.sub(r"_+", "_", base).strip("_")
    return base or None


@dataclass(frozen=True)
class CleanSequenceParams:
    sheet: str
    subj: str
    protocol: str | None
    sequence: object
    group: object
    G: float
    TN: float
    x: float
    y: float
    type: object
    Hz: float
    bmax: float
    N: int
    delta_ms: float
    Delta_app_ms: float
    max_dur_ms: float
    tm_ms: float
    td_ms: float
    TE: float
    TR: float
    g_thorsten: float | None

    def as_columns(self) -> dict[str, object]:
        out = {
            "subj": self.subj,
            "sheet": self.sheet,
            "protocol": self.protocol,
            "sequence": self.sequence,
            "group": self.group,
            "G": self.G,
            "TN": self.TN,
            "x": self.x,
            "y": self.y,
            "type": self.type,
            "Hz": self.Hz,
            "bmax": self.bmax,
            "N": self.N,
            "delta_ms": self.delta_ms,
            "Delta_app_ms": self.Delta_app_ms,
            "max_dur_ms": self.max_dur_ms,
            "tm_ms": self.tm_ms,
            "td_ms": self.td_ms,
            "TE": self.TE,
            "TR": self.TR,
        }
        if self.g_thorsten is not None:
            out["g_thorsten"] = self.g_thorsten
        return out


def detect_gradient_input_kind(stats: dict[str, pd.DataFrame]) -> str:
    any_df = next(iter(stats.values()))
    cols = [str(c) for c in any_df.columns]

    if first_present(cols, ["gvalues", "gval", "g"]) is not None:
        return "g"
    if first_present(cols, ["bvalues", "bvalue", "bval", "b"]) is not None:
        return "b"
    raise ValueError(f"Could not detect a supported gradient column in results table. Columns={cols}")


def is_single_point_results(stats: dict[str, pd.DataFrame]) -> bool:
    return all(len(df.index) == 1 for df in stats.values())


def find_matching_gradient_vector(results_file: Path) -> Path | None:
    seq_name = results_file.name.removesuffix("_results.xlsx")
    project_root = infer_project_root(results_file)
    search_root = project_root / "Data-NIFTI"
    if not search_root.is_dir():
        return None

    for ext in (".gvec", ".bvec"):
        matches = sorted(search_root.rglob(f"{seq_name}{ext}"))
        if matches:
            return matches[0]
    return None


def infer_project_root(path: Path) -> Path:
    resolved = path.resolve()
    for parent in [resolved.parent] + list(resolved.parents):
        if (parent / "Data-NIFTI").is_dir():
            return parent
        if parent.name == "Data-signals":
            return parent.parent
    return resolved.parent


def direction_from_vector_path(vector_path: Path | None) -> str:
    if vector_path is None:
        return "1"

    arr = np.loadtxt(vector_path, dtype=float)
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1:
        if arr.size != 3:
            raise ValueError(f"Expected 3 values or 3xN values in {vector_path}, got shape={arr.shape}")
        vec = arr
    else:
        if arr.shape[0] != 3:
            raise ValueError(f"Expected 3 rows in {vector_path}, got shape={arr.shape}")
        vec = np.nanmean(arr, axis=1)

    avec = np.abs(vec)
    axis = int(np.nanargmax(avec))
    if not np.isfinite(avec[axis]):
        return "1"

    return str(axis + 1)


def aggregate_g_results(
    stats: dict[str, pd.DataFrame],
    *,
    direction_label: str,
    source_file: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    g_value_ref: float | None = None

    for stat, df in stats.items():
        gcol = first_present(list(df.columns), ["gvalues", "gval", "g"])
        if gcol is None:
            raise ValueError(f"Missing g-values column in sheet/stat {stat}. Columns={list(df.columns)}")

        rois = [c for c in df.columns if c != gcol]
        gvals = pd.to_numeric(df[gcol], errors="coerce").dropna().unique().tolist()
        if len(gvals) != 1:
            raise ValueError(f"Expected a single unique g value in sheet/stat {stat}, found {gvals}.")

        g_value = float(gvals[0])
        if g_value_ref is None:
            g_value_ref = g_value
        elif not np.isclose(g_value_ref, g_value, rtol=0.0, atol=1e-9):
            raise ValueError(f"Inconsistent g values across stat sheets: first={g_value_ref}, current={g_value} ({stat}).")

        summary = df[rois].apply(pd.to_numeric, errors="coerce").mean(axis=0)
        b_step = 0 if np.isclose(g_value, 0.0, rtol=0.0, atol=1e-9) else 1

        for roi in rois:
            rows.append(
                {
                    "stat": stat,
                    "roi": str(roi),
                    "direction": str(direction_label),
                    "b_step": int(b_step),
                    "bvalue": np.nan,
                    "g": float(g_value),
                    "value": float(summary.loc[roi]),
                    "source_file": source_file,
                    "gradient_axis_kind": "g",
                }
            )

    return pd.DataFrame(rows)


def assign_group_bsteps(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    g_numeric = pd.to_numeric(out["g"], errors="coerce")
    unique_g = np.sort(pd.unique(g_numeric.dropna()))

    step_map: dict[float, int] = {}
    next_step = 1
    for gv in unique_g:
        if np.isclose(float(gv), 0.0, rtol=0.0, atol=1e-9):
            step_map[float(gv)] = 0
        else:
            step_map[float(gv)] = next_step
            next_step += 1

    def step_for_value(value) -> int:
        fv = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        if not np.isfinite(fv):
            return pd.NA
        for key, step in step_map.items():
            if np.isclose(float(fv), float(key), rtol=0.0, atol=1e-9):
                return int(step)
        return pd.NA

    out["b_step"] = out["g"].map(step_for_value)
    return out


def add_direct_g_derivatives(df: pd.DataFrame, *, gamma: float) -> pd.DataFrame:
    out = df.copy()
    out = ensure_numeric(out, ["g", "b_step", "N", "delta_ms", "Delta_app_ms", "g_thorsten"])

    if "g" not in out.columns:
        return out

    g_numeric = pd.to_numeric(out["g"], errors="coerce")
    if g_numeric.dropna().empty:
        return out

    gmax = out.groupby("direction")["g"].transform("max")
    out["g_max"] = gmax

    b_step = pd.to_numeric(out["b_step"], errors="coerce")
    b_step_max = b_step.groupby(out["direction"]).transform("max")
    with np.errstate(divide="ignore", invalid="ignore"):
        out["g_lin_max"] = gmax * (b_step / b_step_max)
    out.loc[b_step_max <= 0, "g_lin_max"] = g_numeric

    if "g_thorsten" not in out.columns:
        out["g_thorsten"] = pd.NA

    n_value = single_numeric_value(out, "N")
    delta_ms = single_numeric_value(out, "delta_ms")
    delta_app_ms = single_numeric_value(out, "Delta_app_ms")
    if n_value is None or delta_ms is None or delta_app_ms is None:
        for col in ["bvalue_g", "bvalue_g_lin_max", "bvalue_thorsten"]:
            if col not in out.columns:
                out[col] = pd.NA
        return out

    out["bvalue_g"] = b_from_g(
        pd.to_numeric(out["g"], errors="coerce").to_numpy(float),
        N=n_value,
        gamma=gamma,
        delta_ms=delta_ms,
        delta_app_ms=delta_app_ms,
        g_type="g",
    )
    out["bvalue_g_lin_max"] = b_from_g(
        pd.to_numeric(out["g_lin_max"], errors="coerce").to_numpy(float),
        N=n_value,
        gamma=gamma,
        delta_ms=delta_ms,
        delta_app_ms=delta_app_ms,
        g_type="g_lin_max",
    )
    out["bvalue_thorsten"] = b_from_g(
        pd.to_numeric(out["g_thorsten"], errors="coerce").to_numpy(float),
        N=n_value,
        gamma=gamma,
        delta_ms=delta_ms,
        delta_app_ms=delta_app_ms,
        g_type="g_thorsten",
    )
    return out


def single_numeric_value(df: pd.DataFrame, col: str) -> float | None:
    if col not in df.columns:
        return None
    values = pd.to_numeric(df[col], errors="coerce").dropna().unique()
    if len(values) != 1:
        return None
    return float(values[0])


def single_text_value(df: pd.DataFrame, col: str) -> str | None:
    if col not in df.columns:
        return None
    values = df[col].dropna().astype(str).str.strip().unique()
    values = [v for v in values if v]
    if len(values) != 1:
        return None
    return str(values[0])


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def master_analysis_id(df: pd.DataFrame, *, row_kind: str) -> str:
    preferred = (
        "subj",
        "sheet",
        "type",
        "protocol",
        "group",
        "TN",
        "N",
        "td_ms",
        "Hz",
        "G",
        "sequence",
    )
    unique_cols = []
    for col in preferred:
        if col not in df.columns:
            continue
        values = pd.Series(df[col]).dropna().unique()
        if len(values) == 1:
            unique_cols.append(col)
    return build_analysis_id_from_columns(df, columns=unique_cols, prefix=row_kind)


def filter_existing_to_incoming_curve(existing: pd.DataFrame, incoming: pd.DataFrame) -> pd.DataFrame:
    out = existing.copy()

    for col in ["sheet", "subj", "type"]:
        value = single_text_value(incoming, col)
        if value is not None and col in out.columns:
            out = out.loc[out[col].astype(str).str.strip() == value].copy()

    for col in ["group", "TN", "N"]:
        value = single_numeric_value(incoming, col)
        if value is not None and col in out.columns:
            existing_values = pd.to_numeric(out[col], errors="coerce")
            out = out.loc[np.isclose(existing_values.to_numpy(float), value, rtol=0.0, atol=1e-6)].copy()

    return out


def merge_into_group_curve(existing: pd.DataFrame | None, incoming: pd.DataFrame, *, gamma: float) -> pd.DataFrame:
    if existing is None or existing.empty:
        merged = incoming.copy()
    else:
        existing = filter_existing_to_incoming_curve(existing, incoming)
        src_files = set(incoming["source_file"].astype(str).tolist())
        keep_existing = ~existing["source_file"].astype(str).isin(src_files)
        merged = pd.concat([existing.loc[keep_existing].copy(), incoming.copy()], ignore_index=True)

    merged = assign_group_bsteps(merged)
    merged = add_direct_g_derivatives(merged, gamma=gamma)
    drop_cols = [c for c in ["S0", "value_norm"] if c in merged.columns]
    if drop_cols:
        merged = merged.drop(columns=drop_cols)
    merged = add_s0_and_value_norm(merged)
    merged = finalize_clean_signal_long(merged)
    merged = merged.sort_values(["stat", "roi", "direction", "b_step", "g"], kind="stable").reset_index(drop=True)
    return merged


def extract_clean_sequence_params(row: pd.Series, meta) -> CleanSequenceParams:
    sheet = str(row_get(row, ["sheet"], meta.sheet))
    subj = row_get(row, ["subj"], None)
    if subj is None or not str(subj).strip():
        raise ValueError(
            "Matched sequence_parameters row is missing a required 'subj' value. "
            f"sheet={sheet!r}, seq={meta.seq}"
        )
    subj = str(subj).strip()

    group_value = row_get(row, ["group"], getattr(meta, "group", None))
    g_value = to_float(row_get(row, ["G"], getattr(meta, "G", np.nan)))
    tn_value = to_float(row_get(row, ["TN"], getattr(meta, "TN", np.nan)))
    protocol_raw = row_get(row, ["protocol", "Protocol*"], None)
    protocol = normalize_protocol_label(protocol_raw, group_value, tn_value)

    max_dur_ms = to_float(row_get(row, ["max_dur_ms", "d_ms", "Max duration d  [ms]"], None))
    tm_ms = to_float(row_get(row, ["tm_ms", "TM_ms", "mixing time TM  [ms]"], None))
    td_ms = to_float(row_get(row, ["td_ms", "TN"], np.nan))
    if not np.isfinite(td_ms) and np.isfinite(max_dur_ms) and np.isfinite(tm_ms):
        td_ms = 2.0 * max_dur_ms + tm_ms

    g_thorsten = to_float(row_get(row, ["g_thorsten", "G thorsten [mT/m]"], np.nan))
    if not np.isfinite(g_thorsten):
        g_thorsten = None
    else:
        g_thorsten = float(np.sqrt(2.0) * abs(g_thorsten))

    return CleanSequenceParams(
        sheet=sheet,
        subj=subj,
        protocol=protocol,
        sequence=row_get(row, ["sequence", "seq"], None),
        group=group_value,
        G=g_value,
        TN=tn_value,
        x=to_float(row_get(row, ["x"], np.nan)),
        y=to_float(row_get(row, ["y"], np.nan)),
        type=row_get(row, ["type", "seq_type"], None),
        Hz=to_float(row_get(row, ["Hz", "Frecuency [Hz]"], meta.Hz)),
        bmax=to_float(row_get(row, ["bmax", "bval_max [s/mm2]"], meta.bmax)),
        N=int(to_float(row_get(row, ["N"], 1))),
        delta_ms=to_float(row_get(row, ["delta_ms", "delta [ms]"], meta.delta_ms)),
        Delta_app_ms=to_float(row_get(row, ["Delta_app_ms", "delta_app_ms", "Delta_app [ms]"], meta.Delta_ms)),
        max_dur_ms=max_dur_ms,
        tm_ms=tm_ms,
        td_ms=td_ms,
        TE=to_float(row_get(row, ["TE", "TE_ms", "Echo time TE  [ms]"], np.nan)),
        TR=to_float(row_get(row, ["TR", "TR_ms", "Repetition time TR  [ms]"], np.nan)),
        g_thorsten=g_thorsten,
    )


def add_clean_sequence_params(df: pd.DataFrame, params: CleanSequenceParams) -> pd.DataFrame:
    out = df.copy()
    for col, value in params.as_columns().items():
        out[col] = value
    return out


def add_s0_and_value_norm(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out = ensure_numeric(out, ["b_step", "value"])

    s0 = (
        out.loc[out["b_step"] == 0]
        .groupby(["stat", "roi", "direction"], as_index=False)["value"]
        .mean()
        .rename(columns={"value": "S0"})
    )
    out = out.merge(s0, on=["stat", "roi", "direction"], how="left")

    with np.errstate(divide="ignore", invalid="ignore"):
        out["value_norm"] = out["value"] / out["S0"]

    m0 = (out["b_step"] == 0) & out["S0"].notna()
    out.loc[m0, "value_norm"] = 1.0
    return out


def add_g_and_b_derivatives(
    df: pd.DataFrame,
    *,
    gamma: float,
    N: int,
    delta_ms: float,
    Delta_app_ms: float,
    g_thorsten: float | None,
) -> pd.DataFrame:
    out = df.copy()
    out = ensure_numeric(out, ["bvalue", "b_step"])

    b = out["bvalue"].to_numpy(dtype=float)
    denom = N * (gamma**2) * (delta_ms**2) * Delta_app_ms
    g = np.sqrt((b * 1e9) / denom)
    g[np.isclose(b, 0.0)] = 0.0
    out["g"] = g

    b_step_max = int(np.nanmax(out["b_step"].to_numpy(dtype=float)))
    if not np.isfinite(b_step_max) or b_step_max <= 0:
        b_step_max = 1

    gmax_by_dir = out.loc[out["b_step"] == b_step_max].groupby("direction")["g"].max()
    out["g_max"] = out["direction"].map(gmax_by_dir)
    out["g_max"] = out["g_max"].fillna(np.nanmax(out["g"].to_numpy(dtype=float)))

    out["g_lin_max"] = out["g_max"] * (out["b_step"] / float(b_step_max))

    if g_thorsten is None or not np.isfinite(g_thorsten):
        out["g_thorsten"] = np.nan
    else:
        out["g_thorsten"] = float(g_thorsten) * (out["b_step"] / float(b_step_max))

    out["bvalue_g"] = b_from_g(
        pd.to_numeric(out["g"], errors="coerce").to_numpy(float),
        N=N,
        gamma=gamma,
        delta_ms=delta_ms,
        delta_app_ms=Delta_app_ms,
        g_type="g",
    )
    out["bvalue_g_lin_max"] = b_from_g(
        pd.to_numeric(out["g_lin_max"], errors="coerce").to_numpy(float),
        N=N,
        gamma=gamma,
        delta_ms=delta_ms,
        delta_app_ms=Delta_app_ms,
        g_type="g_lin_max",
    )
    out["bvalue_thorsten"] = b_from_g(
        pd.to_numeric(out["g_thorsten"], errors="coerce").to_numpy(float),
        N=N,
        gamma=gamma,
        delta_ms=delta_ms,
        delta_app_ms=Delta_app_ms,
        g_type="g_thorsten",
    )
    return out
