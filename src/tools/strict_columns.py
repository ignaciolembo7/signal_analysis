from __future__ import annotations

from collections.abc import Iterable


UNRECOGNIZED_COLUMN_BASES = {
    "axis",
    "signal_norm",
    "gthorsten",
    "gthorsten_mTm",
    "bvalue_gthorsten",
}

SUFFIXED_UNRECOGNIZED_COLUMN_BASES = {
    "axis",
    "signal_norm",
    "gthorsten",
    "bvalue_gthorsten",
}

CANONICAL_COLUMN_HINT = "Use canonical names such as 'direction', 'value_norm', and 'g_thorsten'."


def find_unrecognized_column_names(columns: Iterable[str]) -> list[str]:
    found: list[str] = []
    for column in columns:
        if column in UNRECOGNIZED_COLUMN_BASES:
            found.append(str(column))
            continue
        if any(
            column == f"{base}_1" or column == f"{base}_2"
            for base in SUFFIXED_UNRECOGNIZED_COLUMN_BASES
        ):
            found.append(str(column))
    return sorted(found)


def raise_on_unrecognized_column_names(columns: Iterable[str], *, context: str) -> None:
    found = find_unrecognized_column_names(columns)
    if found:
        raise ValueError(f"{context}: unrecognized column names: {found}. {CANONICAL_COLUMN_HINT}")
