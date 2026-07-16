from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

def _normalize_header_token(value) -> str:
    token = str(value).strip().lower()
    token = token.replace("*", "")
    token = " ".join(token.split())
    return token


def _find_protocol_header_row(df: pd.DataFrame) -> int:
    for idx in range(len(df)):
        row = {_normalize_header_token(v) for v in df.iloc[idx].tolist() if pd.notna(v)}
        if "protocol" in row and ("seq" in row or "sequence" in row):
            return int(idx)
    raise ValueError("Could not find a header row with protocol/sequence columns in this sheet.")


def _already_has_header(df: pd.DataFrame) -> bool:
    cols = {_normalize_header_token(c) for c in df.columns if pd.notna(c)}
    return "protocol" in cols and ("seq" in cols or "sequence" in cols)

def read_sequence_params_xlsx(path: str | Path) -> pd.DataFrame:
    """
    Read the sequence-parameters workbook into one tidy row per (sheet, protocol, seq, Hz).
    Extra columns such as a manually curated 'subj' tag are preserved.
    """
    xls = pd.ExcelFile(path)
    out = []

    for sheet in xls.sheet_names:
        sub = pd.read_excel(xls, sheet_name=sheet)

        if not _already_has_header(sub):
            raw = pd.read_excel(xls, sheet_name=sheet, header=None)
            hdr = _find_protocol_header_row(raw)
            sub = raw.iloc[hdr:, :].copy()
            sub.columns = sub.iloc[0].tolist()
            sub = sub.iloc[1:].reset_index(drop=True)

        # Drop empty columns.
        sub = sub.loc[:, [c for c in sub.columns if str(c) != "nan"]]

        # Normalize names required by the pipeline.
        rename = {
            "Protocol*": "protocol",
            "Protocol": "protocol",
            "seq": "seq",
            "sequence": "seq",
            "Frecuency [Hz]": "Hz",
            "bval_max [s/mm2]": "bmax",
            "delta [ms]": "delta_ms",
            "Delta_app [ms]": "delta_app_ms",
            "N": "N",
            "G thorsten [mT/m]": "g_thorsten",
            "type of seq": "seq_type",
            # Canonical name: max_dur_ms (formerly d_ms).
            "max duration d  [ms]": "max_dur_ms",
            "Echo time TE  [ms]": "TE_ms",
            "Repetition time TR  [ms]": "TR_ms",
            "mixing time TM  [ms]": "TM_ms",
        }
        sub = sub.rename(columns={k: v for k, v in rename.items() if k in sub.columns})

        sub["sheet"] = sheet

        # Clean numeric types; many cells arrive as '-' or strings.
        for c in [
            "seq", "Hz", "bmax", "delta_ms", "delta_app_ms", "N",
            "group", "G", "TN", "x", "y",
            "g_thorsten", "max_dur_ms", "TE_ms", "TR_ms", "TM_ms"
        ]:
            if c in sub.columns:
                sub[c] = pd.to_numeric(sub[c], errors="coerce")

        # Notebook convention: PGSE uses N=1 when the value is empty.
        if "seq_type" in sub.columns and "N" in sub.columns:
            is_pgse = sub["seq_type"].astype(str).str.contains("PGSE", case=False, na=False)
            sub.loc[is_pgse & sub["N"].isna(), "N"] = 1

        out.append(sub)

    params = pd.concat(out, ignore_index=True)

    # Keep relevant columns first when they exist.
    col_order = [
        c for c in [
            "sheet","subj","protocol","seq","group","G","TN","Hz","bmax","max_dur_ms",
            "delta_ms","delta_app_ms","N","x","y","g_thorsten",
            "type","seq_type","TE_ms","TR_ms","TM_ms"
        ] if c in params.columns
    ]
    rest = [c for c in params.columns if c not in col_order]
    return params[col_order + rest]
