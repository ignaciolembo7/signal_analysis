from __future__ import annotations

from pathlib import Path
import re


def canonical_sheet_name(name: str | None) -> str | None:
    if name is None:
        return None
    s = str(name).strip()
    if not s:
        return None

    m = re.match(r"^(\d{8}[-_][^_]+)", s)
    if m:
        return m.group(1)

    m = re.match(r"^(.+?)_(?:N\d|td)", s)
    if m:
        return m.group(1)

    return s


def infer_brain_group(sheet: str | None, source_name: str | None = None) -> str:
    raw = str(sheet or source_name or "").strip()
    if not raw:
        return "UNKNOWN"

    stem = Path(raw).stem
    match = re.match(r"^\d{8}[-_](.+)$", stem)
    tail = match.group(1) if match else stem
    token = tail.split("_")[0]
    token = re.sub(r"-\d+$", "", token)
    return token or stem


def infer_subj_label(sheet: str | None, source_name: str | None = None) -> str:
    token = infer_brain_group(sheet, source_name=source_name)
    token = re.sub(r"(?<=[A-Za-z])\d+$", "", token)
    token = str(token).strip()
    return token or "UNKNOWN"
