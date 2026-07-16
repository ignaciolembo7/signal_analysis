from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re
import pandas as pd

DEFAULT_STAT_SHEETS = {0: "avg", 1: "std", 2: "med", 3: "mad", 4: "mode"}
LEGACY_XLSX_STAT_SHEETS = {
    "Sheet1": "avg",
    "Sheet2": "std",
    "Sheet3": "med",
    "Sheet4": "mad",
    "Sheet5": "mode",
}
NAMED_STAT_SHEETS = {
    "mean": "avg",
    "avg": "avg",
    "std": "std",
    "median": "med",
    "med": "med",
    "mad": "mad",
    "mode": "mode",
}
@dataclass(frozen=True)
class Layout:
    nbvals: int | None
    ndirs: int | None

def infer_layout_from_filename(path: str | Path) -> Layout:
    p = Path(path)
    m_bval = re.search(r"(\d+)bval", p.name)
    m_dir  = re.search(r"(\d+)(?:ortho)?dir", p.name, flags=re.IGNORECASE)
    nbvals = int(m_bval.group(1)) if m_bval else None
    ndirs  = int(m_dir.group(1)) if m_dir else None
    return Layout(nbvals=nbvals, ndirs=ndirs)

NEW_STAT_SHEETS = {"Mean": "avg", "Min": "min", "Max": "max", "Area": "area"}


def _read_xlsx_stats(xls: pd.ExcelFile) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    sheet_names = list(xls.sheet_names)
    exact_map = {str(name): str(name) for name in sheet_names}
    norm_map = {str(name).strip().lower(): str(name) for name in sheet_names}

    # New contrast/area workbook format.
    for sheet_name, stat_name in NEW_STAT_SHEETS.items():
        if sheet_name in exact_map:
            out[stat_name] = pd.read_excel(xls, sheet_name=exact_map[sheet_name], engine="openpyxl")
    if out:
        return out

    # MATLAB-style output written as Sheet1..Sheet5.
    for sheet_name, stat_name in LEGACY_XLSX_STAT_SHEETS.items():
        if sheet_name in exact_map:
            out[stat_name] = pd.read_excel(xls, sheet_name=exact_map[sheet_name], engine="openpyxl")
    if out:
        return out

    # Additional compatibility when sheets are already named by stat.
    for sheet_name, stat_name in NAMED_STAT_SHEETS.items():
        if sheet_name in norm_map:
            out[stat_name] = pd.read_excel(xls, sheet_name=norm_map[sheet_name], engine="openpyxl")
    if out:
        return out

    # Fallback: if nothing matches, read the first sheet as avg.
    return {"avg": pd.read_excel(xls, sheet_name=0, engine="openpyxl")}

def read_result_xls(path: str | Path, stat_sheets=DEFAULT_STAT_SHEETS) -> dict[str, pd.DataFrame]:
    """
    Read legacy (.xls) and new (.xlsx) results into dict[stat -> DataFrame].
    - Legacy .xls: sheets by index 0..4 -> avg/std/med/mad/mode
    - Legacy .xlsx: Sheet1..Sheet5 or mean/std/median/mad/mode sheets
    - New .xlsx: sheets named Mean/Min/Max/Area -> avg/min/max/area
    """
    p = Path(path)
    suf = p.suffix.lower()

    out: dict[str, pd.DataFrame] = {}

    if suf == ".xls":
        # Legacy .xls requires xlrd.
        xls = pd.ExcelFile(p, engine="xlrd")
        for sheet_idx, stat_name in stat_sheets.items():
            out[stat_name] = pd.read_excel(xls, sheet_name=sheet_idx, engine="xlrd")
        return out

    # New workbook format: .xlsx.
    xls = pd.ExcelFile(p, engine="openpyxl")
    return _read_xlsx_stats(xls)


def read_table_file(path: str | Path) -> pd.DataFrame:
    """Read one physical table file as a DataFrame."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(p)
    if suffix == ".csv":
        return pd.read_csv(p)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(p, sheet_name=0)
    raise ValueError(f"Unsupported table format for {p}")


def write_table_outputs(
    df: pd.DataFrame,
    parquet_path: str | Path,
    *,
    xlsx_path: str | Path | None = None,
    csv_path: str | Path | None = None,
) -> Path:
    """
    Write a canonical table output bundle, using parquet as the primary artifact.
    Optional xlsx/csv siblings preserve the repository's existing output contracts.
    """
    out_parquet = Path(parquet_path)
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_parquet, index=False)

    if xlsx_path is not None:
        out_xlsx = Path(xlsx_path)
        out_xlsx.parent.mkdir(parents=True, exist_ok=True)
        _EXCEL_MAX_COLS = 16384
        if df.shape[1] > _EXCEL_MAX_COLS:
            import warnings
            warnings.warn(
                f"Skipping {out_xlsx.name}: {df.shape[1]} columns exceeds Excel limit "
                f"({_EXCEL_MAX_COLS}). Parquet output is complete.",
                stacklevel=2,
            )
        else:
            df.to_excel(out_xlsx, index=False)

    if csv_path is not None:
        out_csv = Path(csv_path)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)

    return out_parquet


def sanitize_output_token(value: object) -> str:
    text = str(value).strip()
    if not text:
        return "NA"
    token = re.sub(r"[^A-Za-z0-9._-]+", "-", text)
    token = token.strip("-")
    return token or "NA"


def _direction_token(values: list[str] | tuple[str, ...] | None) -> str:
    if values is None:
        return "direction_ALL"
    normalized = sorted({sanitize_output_token(v) for v in values if str(v).strip()})
    if not normalized:
        return "direction_ALL"
    if len(normalized) == 1:
        return f"direction_{normalized[0]}"
    return f"direction_{'-'.join(normalized)}"


def fit_params_output_basename(
    *,
    model: str,
    axis: str,
    ycol: str,
    directions: list[str] | tuple[str, ...] | None,
) -> str:
    return ".".join(
        [
            "fit_params",
            sanitize_output_token(model),
            sanitize_output_token(axis),
            sanitize_output_token(ycol),
            _direction_token(directions),
        ]
    )


def write_xlsx_csv_outputs(
    df: pd.DataFrame,
    xlsx_path: str | Path,
    *,
    csv_path: str | Path | None = None,
    sheet_name: str | None = None,
) -> Path:
    """Write an xlsx table and its optional csv sibling with consistent directory handling."""
    out_xlsx = Path(xlsx_path)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    if sheet_name is None:
        df.to_excel(out_xlsx, index=False)
    else:
        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    if csv_path is not None:
        out_csv = Path(csv_path)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)

    return out_xlsx


def write_xlsx_sheets(
    tables: dict[str, pd.DataFrame],
    xlsx_path: str | Path,
) -> Path:
    """Write a multi-sheet xlsx workbook with consistent directory handling."""
    out_xlsx = Path(xlsx_path)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        for sheet_name, df in tables.items():
            df.to_excel(writer, sheet_name=str(sheet_name)[:31], index=False)
    return out_xlsx
