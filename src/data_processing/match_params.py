from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re
import pandas as pd
import numpy as np

@dataclass(frozen=True)
class ResultMeta:
    sheet: str | None
    seq: int | None
    Hz: float | None
    bmax: float | None
    group: int | None
    G: float | None
    TN: float | None
    N: float | None

    # Legacy: _d55 is treated as max_dur_ms.
    d_ms: float | None

    # New style: _d1.5 (delta) and _Delta11.6 (Delta_app).
    delta_ms: float | None
    Delta_ms: float | None

    ndirs: int | None
    nbvals: int | None
    encoding: str | None


def parse_results_filename(path: str | Path) -> ResultMeta:
    p = Path(path)
    name = p.name

    parts_lower = [part.lower() for part in p.parts]
    if "results" in parts_lower and parts_lower.index("results") + 1 < len(p.parts) - 1:
        sheet = p.parts[parts_lower.index("results") + 1].strip()
    else:
        parent_name = p.parent.name.strip()
        if parent_name and parent_name.lower() != "results":
            sheet = parent_name
        elif "_ep2d" in name:
            sheet = name.split("_ep2d")[0]
        elif "_d" in name:
            sheet = name.split("_d")[0]
        elif "_Delta" in name:
            sheet = name.split("_Delta")[0]
        else:
            sheet = p.stem.split("_")[0]
    sheet = re.sub(r"_(\d+)bval_(\d+)dir.*$", "", sheet, flags=re.IGNORECASE)

    def _int(rx: str) -> int | None:
        m = re.search(rx, name, re.IGNORECASE)
        return int(m.group(1)) if m else None

    def _float(rx: str) -> float | None:
        m = re.search(rx, name, re.IGNORECASE)
        return float(m.group(1)) if m else None
    
    def _float_with_p(rx: str) -> float | None:
        m = re.search(rx, name, re.IGNORECASE)
        if not m:
            return None
        return float(m.group(1).replace("p", "."))

    nbvals = _int(r"(\d+)bval")
    ndirs  = _int(r"(\d+)(?:ortho)?dir")
    group = _int(r"_(\d{3})_(?:NOGSE|PGSE)")
    G = _float_with_p(r"_G(\d+(?:p\d+|\.\d+)?)")
    TN = _float_with_p(r"_TN(\d+(?:p\d+|\.\d+)?)")
    N = _float_with_p(r"_N(\d+(?:p\d+|\.\d+)?)")

    # New style: interpret _d as delta_ms only when _Delta is present.
    delta_ms = _float(r"_d(\d+(?:\.\d+)?)") if re.search(r"_Delta", name, re.IGNORECASE) else None
    Delta_ms = _float(r"_Delta(\d+(?:\.\d+)?)")

    # Legacy style: _d55 is used only when there is no _Delta token.
    d_ms = _float_with_p(r"_d(\d+(?:p\d+|\.\d+)?)") if delta_ms is None else None

    Hz   = _float(r"Hz(\d+(?:\.\d+)?)")
    bmax = _float(r"_b(\d+(?:\.\d+)?)")

    seq  = _int(r"_(\d+)_results")
    encoding = "OGSE" if re.search(r"OGSE", name, re.IGNORECASE) else ("PGSE" if re.search(r"PGSE", name, re.IGNORECASE) else None)

    return ResultMeta(
        sheet=sheet, seq=seq, Hz=Hz, bmax=bmax, group=group, G=G, TN=TN, N=N,
        d_ms=(float(d_ms) if d_ms is not None else None),
        delta_ms=delta_ms, Delta_ms=Delta_ms,
        ndirs=ndirs, nbvals=nbvals, encoding=encoding
    )


def _norm_sheet(x: str) -> str:
    s = str(x).strip().upper()
    s = s.replace("-", "_")
    s = re.sub(r"_(\d+)BVAL_(\d+)DIR.*$", "", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def _filter_close(df: pd.DataFrame, col: str, target: float, *, atol: float = 1e-6) -> pd.DataFrame:
    if col not in df.columns:
        return df
    v = pd.to_numeric(df[col], errors="coerce")
    return df[np.isclose(v.to_numpy(float), float(target), rtol=0.0, atol=float(atol))]


def _sequence_match_column(df: pd.DataFrame) -> str | None:
    for col in ("sequence", "seq"):
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col
    return None


def select_params_row(params: pd.DataFrame, meta: ResultMeta) -> pd.Series | None:
    """
    Strict match policy:
    1) The sheet name must match exactly.
    2) If a group number is available, match it inside that sheet.
    3) If a sequence number is available, match the Excel sequence column inside that sheet/group.
    4) If more than one row still remains, use gradient, Hz, and timing fields to disambiguate.
    5) Do not silently fall back to another sheet.
    
    Returns:
        pd.Series with the matched row, or None if no match is found.
    """
    df = params.copy()

    # Strict sheet match
    if meta.sheet is None or "sheet" not in df.columns:
        return None

    target_sheet = str(meta.sheet).strip()
    sheet_values = df["sheet"].astype(str).str.strip()
    df = df[sheet_values == target_sheet]

    if df.empty:
        return None

    # Group disambiguates direct-g acquisitions whose sequence numbers belong to separate curves.
    if meta.group is not None and "group" in df.columns:
        df = _filter_close(df, "group", float(meta.group), atol=1e-6)

    if meta.seq is not None:
        seq_col = _sequence_match_column(df)
        if seq_col is not None:
            seq_values = pd.to_numeric(df[seq_col], errors="coerce")
            df = df[seq_values == int(meta.seq)]

            if df.empty:
                return None

    if meta.G is not None and "G" in df.columns:
        df = _filter_close(df, "G", float(meta.G), atol=1e-6)

    if meta.TN is not None and "TN" in df.columns:
        df = _filter_close(df, "TN", float(meta.TN), atol=1e-3)

    if meta.N is not None and "N" in df.columns:
        df = _filter_close(df, "N", float(meta.N), atol=1e-6)

    # Hz
    Hz = 0 if (meta.Hz is None and meta.encoding == "PGSE") else meta.Hz
    if Hz is not None and "Hz" in df.columns:
        hzcol = pd.to_numeric(df["Hz"], errors="coerce")
        if float(Hz) == 0.0:
            df = df[(hzcol.isna()) | (np.isclose(hzcol.to_numpy(float), 0.0, rtol=0.0, atol=1e-9))]
        else:
            df = _filter_close(df, "Hz", float(Hz), atol=1e-6)

    # Timing fields
    if meta.delta_ms is not None and "delta_ms" in df.columns:
        df = _filter_close(df, "delta_ms", float(meta.delta_ms), atol=1e-3)

    if meta.Delta_ms is not None and "delta_app_ms" in df.columns:
        df = _filter_close(df, "delta_app_ms", float(meta.Delta_ms), atol=1e-3)

    if meta.d_ms is not None:
        if "max_dur_ms" in df.columns:
            df = _filter_close(df, "max_dur_ms", float(meta.d_ms), atol=1e-3)
        elif "d_ms" in df.columns:
            df = _filter_close(df, "d_ms", float(meta.d_ms), atol=1e-3)

    if df.empty:
        return None

    if len(df) > 1:
        raise ValueError(
            "More than one parameter row matched inside the exact sheet. "
            "Refine the metadata or make the Excel table more specific.\n"
            f"sheet={target_sheet!r}, seq={meta.seq}, Hz={meta.Hz}, d_ms={meta.d_ms}, "
            f"delta_ms={meta.delta_ms}, Delta_ms={meta.Delta_ms}"
        )

    return df.iloc[0]
