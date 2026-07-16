from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile
from typing import Sequence

import pandas as pd

from fitting.cli_common import load_master_input


@dataclass
class RecipeTableSelection:
    """Rows selected for one pipeline step.

    New code should use ``df`` directly. ``paths`` exists only for older
    functions that still expect parquet files.
    """

    df: pd.DataFrame | None
    paths: list[Path]
    source: str
    _temp_dir: tempfile.TemporaryDirectory[str] | None = None

    def cleanup(self) -> None:
        if self._temp_dir is not None:
            self._temp_dir.cleanup()
            self._temp_dir = None


def selected_rows_or_legacy_paths(
    args,
    *,
    legacy_paths: Sequence[Path] | None,
    default_row_kind: str,
    temp_prefix: str,
    signal_rotated: bool | None = None,
    temp_name: str = "master_selection.long.parquet",
) -> RecipeTableSelection:
    """Select rows from master, or keep explicit legacy paths.

    This is a transition helper: it keeps scripts master-table-first while old
    backend functions are gradually converted from parquet paths to DataFrames.
    """
    master_path = getattr(args, "master_parquet", None)
    if master_path is not None:
        rows = load_master_input(args, default_row_kind=default_row_kind, signal_rotated=signal_rotated)
        temp_dir = tempfile.TemporaryDirectory(prefix=temp_prefix)
        path = Path(temp_dir.name) / temp_name
        rows.to_parquet(path, index=False)
        return RecipeTableSelection(df=rows, paths=[path], source="master", _temp_dir=temp_dir)

    paths = [Path(p) for p in (legacy_paths or []) if p is not None]
    if not paths:
        raise ValueError("Pass --master-parquet with selectors, or explicit legacy input path(s).")
    return RecipeTableSelection(df=None, paths=paths, source="legacy")


def selected_rows_or_legacy_table(
    args,
    *,
    legacy_path: Path | None,
    default_row_kind: str,
    temp_prefix: str,
    signal_rotated: bool | None = None,
    temp_name: str = "master_selection.long.parquet",
) -> RecipeTableSelection:
    paths = [] if legacy_path is None else [Path(legacy_path)]
    return selected_rows_or_legacy_paths(
        args,
        legacy_paths=paths,
        default_row_kind=default_row_kind,
        temp_prefix=temp_prefix,
        signal_rotated=signal_rotated,
        temp_name=temp_name,
    )


def selected_rows_or_legacy_dataframe(
    args,
    *,
    legacy_path: Path | None,
    default_row_kind: str,
    signal_rotated: bool | None = None,
) -> RecipeTableSelection:
    """Return selected rows as a DataFrame for backends that no longer need files."""
    master_path = getattr(args, "master_parquet", None)
    if master_path is not None and legacy_path is None:
        rows = load_master_input(args, default_row_kind=default_row_kind, signal_rotated=signal_rotated)
        return RecipeTableSelection(df=rows, paths=[], source="master")

    if legacy_path is None:
        raise ValueError("Pass --master-parquet with selectors, or an explicit legacy input path.")
    path = Path(legacy_path)
    return RecipeTableSelection(df=pd.read_parquet(path), paths=[path], source="legacy")
