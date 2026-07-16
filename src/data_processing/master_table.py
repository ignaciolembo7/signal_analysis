from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
import re

import numpy as np
import pandas as pd

from data_processing.io import write_table_outputs
from data_processing.schema import ensure_dproj_long_schema, ensure_signal_long_schema
from tools.strict_columns import raise_on_unrecognized_column_names


MASTER_ROW_KINDS = ("signal", "signal_rotated", "contrast", "dproj", "fit_params", "fit_points")

SIGNAL_ROW_KINDS = {"signal", "signal_rotated"}
CONTRAST_REQUIRED_COLUMNS = {"stat", "roi", "direction", "b_step", "value", "value_norm"}
MASTER_BASE_COLUMNS = [
    "row_kind",
    "analysis_id",
    "subj",
    "sheet",
    "roi",
    "direction",
    "stat",
    "b_step",
    "value",
    "value_norm",
    "D_proj",
    "grad_correction_factor",
    "alpha_macro",
    "td_ms",
    "N",
    "Hz",
    "source_file",
    "source_path",
    "source_hash",
]


def _as_list(value: object | Sequence[object] | None) -> list[object] | None:
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def split_selector_values(values: object | Sequence[object] | None) -> list[str] | None:
    """Normalize CLI selector values, accepting repeated or comma-separated args."""
    raw_values = _as_list(values)
    if raw_values is None:
        return None
    out: list[str] = []
    for raw in raw_values:
        for token in str(raw).replace(",", " ").split():
            token = token.strip()
            if token:
                out.append(token)
    return out or None


def _is_all(values: Sequence[object] | None) -> bool:
    return values is None or (len(values) == 1 and str(values[0]).upper() == "ALL")


def _normalize_token(value: object) -> str:
    text = str(value).strip()
    if text == "":
        return "NA"
    try:
        num = float(text)
    except ValueError:
        pass
    else:
        if np.isfinite(num):
            if abs(num - round(num)) < 1e-9:
                text = str(int(round(num)))
            else:
                text = f"{num:.6g}"
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", text.replace(".", "p"))
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "NA"


def _coerce_row_kind(row_kind: str | None) -> str:
    kind = "signal" if row_kind is None else str(row_kind).strip()
    if kind not in MASTER_ROW_KINDS:
        raise ValueError(f"Unsupported row_kind={row_kind!r}. Expected one of {list(MASTER_ROW_KINDS)}.")
    return kind


def _with_missing_master_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in MASTER_BASE_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA
    front = [col for col in MASTER_BASE_COLUMNS if col in out.columns]
    rest = [col for col in out.columns if col not in front]
    return out[front + rest]


def _coerce_label_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "direction" in out.columns:
        out["direction"] = out["direction"].map(lambda x: pd.NA if pd.isna(x) else str(x))
    return out


def build_analysis_id_from_columns(
    row_or_df: pd.Series | pd.DataFrame | Mapping[str, object],
    *,
    columns: Sequence[str] = ("subj", "sheet", "type", "td_ms", "N", "Hz", "direction"),
    prefix: str | None = None,
) -> str:
    """
    Build a stable, filename-free analysis id from explicit table columns.

    When a DataFrame is passed, each selected column must have one non-null
    value. This keeps accidental multi-analysis ids from being silently made.
    """
    if isinstance(row_or_df, pd.DataFrame):
        values: dict[str, object] = {}
        for col in columns:
            if col not in row_or_df.columns:
                continue
            unique = pd.Series(row_or_df[col]).dropna().unique()
            if len(unique) == 1:
                values[col] = unique[0]
            elif len(unique) > 1:
                raise ValueError(f"Column {col!r} is not unique; cannot build analysis_id.")
    else:
        values = dict(row_or_df)

    tokens = []
    if prefix:
        tokens.append(_normalize_token(prefix))
    for col in columns:
        value = values.get(col)
        if value is None or pd.isna(value):
            continue
        tokens.append(f"{_normalize_token(col)}-{_normalize_token(value)}")
    if not tokens:
        raise ValueError("Could not build analysis_id because none of the requested columns had values.")
    return "_".join(tokens)


def validate_master_table(df: pd.DataFrame) -> pd.DataFrame:
    """Validate light master-table invariants without changing scientific values."""
    raise_on_unrecognized_column_names(df.columns, context="validate_master_table")
    if "row_kind" not in df.columns:
        raise ValueError("Master table is missing required column 'row_kind'.")

    out = df.copy()
    invalid = sorted(set(out["row_kind"].dropna().astype(str)) - set(MASTER_ROW_KINDS))
    if invalid:
        raise ValueError(f"Unsupported row_kind values in master table: {invalid}.")

    signal_mask = out["row_kind"].astype(str).isin(SIGNAL_ROW_KINDS)
    if signal_mask.any():
        ensure_signal_long_schema(out.loc[signal_mask].drop(columns=["row_kind"], errors="ignore"))

    dproj_mask = out["row_kind"].astype(str).eq("dproj")
    if dproj_mask.any():
        ensure_dproj_long_schema(out.loc[dproj_mask].drop(columns=["row_kind"], errors="ignore"))

    contrast_mask = out["row_kind"].astype(str).eq("contrast")
    if contrast_mask.any():
        missing = CONTRAST_REQUIRED_COLUMNS - set(out.columns)
        if missing:
            raise ValueError(f"Master contrast rows are missing required columns: {sorted(missing)}")

    return _with_missing_master_columns(_coerce_label_columns(out))


def load_master_table(path: str | Path, *, validate: bool = True) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    suffix = p.suffix.lower()
    if suffix == ".parquet":
        df = pd.read_parquet(p)
    elif suffix == ".csv":
        df = pd.read_csv(p)
    else:
        raise ValueError(f"Unsupported master table format: {p}. Use parquet or csv.")
    return validate_master_table(df) if validate else df


def write_master_table(
    df: pd.DataFrame,
    path: str | Path,
    *,
    validate: bool = True,
) -> Path:
    out = validate_master_table(df) if validate else df
    return write_table_outputs(out, path)


def normalize_master_rows(
    rows: pd.DataFrame,
    *,
    row_kind: str | None = None,
    analysis_id: str | None = None,
) -> pd.DataFrame:
    kind = _coerce_row_kind(row_kind)
    out = rows.copy()
    if "row_kind" not in out.columns:
        out["row_kind"] = kind
    elif row_kind is not None:
        out["row_kind"] = out["row_kind"].fillna(kind).replace("", kind)

    kinds = set(out["row_kind"].dropna().astype(str))
    if kinds & SIGNAL_ROW_KINDS:
        mask = out["row_kind"].astype(str).isin(SIGNAL_ROW_KINDS)
        normalized = ensure_signal_long_schema(out.loc[mask].drop(columns=["row_kind"], errors="ignore"))
        normalized["row_kind"] = out.loc[mask, "row_kind"].to_numpy()
        out = pd.concat([normalized, out.loc[~mask]], ignore_index=True, sort=False)
    if "analysis_id" not in out.columns:
        out["analysis_id"] = analysis_id if analysis_id is not None else pd.NA
    elif analysis_id is not None:
        out["analysis_id"] = out["analysis_id"].fillna(analysis_id).replace("", analysis_id)
    return validate_master_table(out)


def append_master_rows(
    master: pd.DataFrame | str | Path | None,
    rows: pd.DataFrame,
    *,
    row_kind: str | None = None,
    analysis_id: str | None = None,
    out_path: str | Path | None = None,
    drop_duplicates: bool = True,
) -> pd.DataFrame:
    new_rows = normalize_master_rows(rows, row_kind=row_kind, analysis_id=analysis_id)
    if master is None:
        combined = new_rows
    elif isinstance(master, (str, Path)):
        combined = pd.concat([load_master_table(master), new_rows], ignore_index=True, sort=False)
    else:
        combined = pd.concat([validate_master_table(master), new_rows], ignore_index=True, sort=False)

    combined = validate_master_table(combined)
    if drop_duplicates:
        combined = combined.drop_duplicates().reset_index(drop=True)
    if out_path is not None:
        write_master_table(combined, out_path)
    return combined


def filter_master_rows(df: pd.DataFrame, **selectors: object) -> pd.DataFrame:
    out = validate_master_table(df)
    for col, raw_values in selectors.items():
        values = _as_list(raw_values)
        if _is_all(values):
            continue
        if col not in out.columns:
            raise KeyError(f"Cannot filter master table by missing column {col!r}.")
        if len(values) == 1 and isinstance(values[0], float):
            out = out[np.isclose(pd.to_numeric(out[col], errors="coerce"), float(values[0]), atol=1e-6)]
        else:
            wanted = {str(v) for v in values}
            out = out[out[col].astype(str).isin(wanted)]
    return out.reset_index(drop=True)


def filter_table_rows(df: pd.DataFrame, **selectors: object) -> pd.DataFrame:
    """Filter a generic table by explicit column selectors.

    Unlike ``filter_master_rows``, this helper does not validate row_kind/schema,
    so it is suitable for cumulative auxiliary tables such as
    master_fit_params.parquet.
    """
    out = df.copy()
    for col, raw_values in selectors.items():
        values = _as_list(raw_values)
        if _is_all(values):
            continue
        if col not in out.columns:
            raise KeyError(f"Cannot filter table by missing column {col!r}.")
        if len(values) == 1:
            value = values[0]
            try:
                num = float(value)
            except (TypeError, ValueError):
                num = None
            if num is not None and np.isfinite(num):
                numeric = pd.to_numeric(out[col], errors="coerce")
                out = out[np.isclose(numeric, num, atol=1e-6, equal_nan=False)]
                continue
        wanted = {str(v) for v in values}
        out = out[out[col].astype(str).isin(wanted)]
    return out.reset_index(drop=True)


def select_signal(
    master: pd.DataFrame,
    *,
    rotated: bool | None = None,
    **selectors: object,
) -> pd.DataFrame:
    kinds = list(SIGNAL_ROW_KINDS)
    if rotated is True:
        kinds = ["signal_rotated"]
    elif rotated is False:
        kinds = ["signal"]
    return filter_master_rows(master, row_kind=kinds, **selectors)


def select_plot_signal(
    master: pd.DataFrame,
    *,
    rotated: bool | None = True,
    **selectors: object,
) -> pd.DataFrame:
    return select_signal(master, rotated=rotated, **selectors)


def select_contrast_pair(
    master: pd.DataFrame,
    *,
    N_1: float | int | None = None,
    N_2: float | int | None = None,
    **selectors: object,
) -> pd.DataFrame:
    out = filter_master_rows(master, row_kind="contrast", **selectors)
    if N_1 is not None:
        out = out[np.isclose(pd.to_numeric(out["N_1"], errors="coerce"), float(N_1), atol=1e-6)]
    if N_2 is not None:
        out = out[np.isclose(pd.to_numeric(out["N_2"], errors="coerce"), float(N_2), atol=1e-6)]
    return out.reset_index(drop=True)


def select_contrast(master: pd.DataFrame, **selectors: object) -> pd.DataFrame:
    return filter_master_rows(master, row_kind="contrast", **selectors)


def load_fit_params_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    suffix = p.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(p)
    if suffix == ".csv":
        return pd.read_csv(p)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(p, sheet_name=0)
    raise ValueError(f"Unsupported fit params table format: {p}")


def select_fit_params(
    table: pd.DataFrame,
    *,
    fit_kind: object | Sequence[object] | None = None,
    model: object | Sequence[object] | None = None,
    ok_only: bool = False,
    **selectors: object,
) -> pd.DataFrame:
    out = table.copy()
    if fit_kind is not None and "fit_kind" in out.columns:
        selectors["fit_kind"] = fit_kind
    if model is not None:
        if "fit_model" in out.columns and "model" not in out.columns:
            selectors["fit_model"] = model
        elif "model" in out.columns:
            selectors["model"] = model
    out = filter_table_rows(out, **selectors)
    if ok_only and "ok" in out.columns:
        out = out[out["ok"].fillna(False).astype(bool)].copy()
    return out.reset_index(drop=True)


def select_alpha_macro(table: pd.DataFrame, **selectors: object) -> pd.DataFrame:
    out = table.copy()
    if "alpha_macro" not in out.columns and "alpha_inf" in out.columns:
        out = out.rename(columns={"alpha_inf": "alpha_macro"})
    if "alpha_macro_error" not in out.columns and "alpha_macro_se" in out.columns:
        out = out.rename(columns={"alpha_macro_se": "alpha_macro_error"})
    if "alpha_macro_error" not in out.columns and "alpha_inf_se" in out.columns:
        out = out.rename(columns={"alpha_inf_se": "alpha_macro_error"})

    possible_fit_kind = {"alpha_macro", "alpha_macro_summary", "dproj_alpha_macro"}
    if "fit_kind" in out.columns:
        mask = out["fit_kind"].astype(str).isin(possible_fit_kind)
        if mask.any():
            out = out[mask].copy()
    elif "row_kind" in out.columns:
        mask = out["row_kind"].astype(str).isin(possible_fit_kind | {"fit_params"})
        if mask.any():
            out = out[mask].copy()

    out = filter_table_rows(out, **selectors)
    if "alpha_macro" not in out.columns:
        raise KeyError("Selected table does not contain alpha_macro.")
    out["alpha_macro"] = pd.to_numeric(out["alpha_macro"], errors="coerce")
    if "alpha_macro_error" in out.columns:
        out["alpha_macro_error"] = pd.to_numeric(out["alpha_macro_error"], errors="coerce")
    return out.dropna(subset=["alpha_macro"]).reset_index(drop=True)
